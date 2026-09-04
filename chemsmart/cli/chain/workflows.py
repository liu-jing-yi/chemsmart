"""Predefined chain workflow subcommands (pKa, Fukui, redox, reaction).

These are Click groups under ``chain``, not ``-s/--steps`` pipeline jobs.
Submit requires ``--program {gaussian,orca}`` and hydrates ``ctx.obj`` from
the chain YAML alias so existing program builders can run unchanged.
Analyze is attached for pKa, Fukui, and redox only.
"""

import functools
import importlib
import logging

import click
from click.core import ParameterSource

from chemsmart.analysis.fukui import FUKUI_MODES, ORCA_FUKUI_MODES
from chemsmart.cli.fukui import fukui as fukui_analyze
from chemsmart.cli.job import click_job_options
from chemsmart.cli.pka import analyze as pka_analyze
from chemsmart.cli.pka import batch_analyze as pka_batch_analyze
from chemsmart.cli.pka import (
    click_pka_proton_options,
    click_pka_shared_options,
    resolve_pka_entropy_cutoff,
)
from chemsmart.cli.reaction import (
    click_reaction_shared_options,
    merge_reaction_options,
    store_reaction_shared,
)
from chemsmart.cli.redox import analyze as redox_analyze
from chemsmart.cli.redox import click_redox_shared_options, store_redox_shared
from chemsmart.cli.utils import (
    CHAIN_CLI_DEFAULTS_KEY,
    CHAIN_PROJECT_SETTINGS_KEY,
)
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import return_objects_and_indices_from_string_index

logger = logging.getLogger(__name__)

CHAIN_WORKFLOW_COMMANDS = ("pka", "fukui", "redox", "reaction")
CHAIN_WORKFLOW_PROGRAMS = ("gaussian", "orca")
_ANALYZE_SUBCOMMANDS = frozenset({"analyze", "batch-analyze"})
_TABLE_SUFFIXES = (".csv", ".tsv", ".txt")
_WORKFLOW_PROGRAM_KEY = "chain_workflow_program"


def click_workflow_program_option(f):
    """``--program`` for chain workflow submit (not required for analyze)."""

    @click.option(
        "--program",
        type=click.Choice(CHAIN_WORKFLOW_PROGRAMS, case_sensitive=False),
        required=False,
        help=(
            "QC program for submit. Theory and solvent come from that "
            "program's alias in the chain project YAML."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_fukui_submit_options(f):
    """Fukui submit options shared by Gaussian and ORCA program commands."""

    @click.option(
        "--mode",
        default="mulliken",
        show_default=True,
        type=click.Choice(list(FUKUI_MODES), case_sensitive=False),
        help="Charges to be used for Fukui Indices calculations.",
    )
    @click.option(
        "-rcc",
        "--radical-cation-charge",
        type=int,
        default=None,
        help=(
            "Override charge for the radical-cation job. "
            "Default is derived from the neutral charge."
        ),
    )
    @click.option(
        "-rcm",
        "--radical-cation-multiplicity",
        type=int,
        default=None,
        help=(
            "Override multiplicity for the radical-cation job. "
            "Default is derived from the neutral multiplicity."
        ),
    )
    @click.option(
        "-rac",
        "--radical-anion-charge",
        type=int,
        default=None,
        help=(
            "Override charge for the radical-anion job. "
            "Default is derived from the neutral charge."
        ),
    )
    @click.option(
        "-ram",
        "--radical-anion-multiplicity",
        type=int,
        default=None,
        help=(
            "Override multiplicity for the radical-anion job. "
            "Default is derived from the neutral multiplicity."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def require_workflow_program(program):
    """Require ``--program`` for workflow submit."""
    if program is None:
        raise click.UsageError(
            "--program is required for submit (gaussian or orca)."
        )
    return program.lower()


def require_chain_project(ctx, action="workflow submit"):
    """Require chain ``-p/--project`` settings on ``ctx.obj``."""
    settings = None
    if ctx.obj is not None:
        settings = ctx.obj.get(CHAIN_PROJECT_SETTINGS_KEY)
    if settings is None:
        raise click.UsageError(f"-p/--project is required for {action}.")
    return settings


def hydrate_program_ctx_from_chain(ctx, program, job_name=None):
    """Fill gaussian/orca builder keys from the chain YAML alias.

    Loads ``GaussianProjectSettings`` / ``ORCAProjectSettings`` via
    ``project_for(program)``, builds default job settings, and sets
    ``project_settings``, ``job_settings``, ``keywords``, ``molecules``,
    ``molecule_indices``, ``filename``, and ``label`` on ``ctx.obj``.
    Charge and multiplicity come from chain ``-c/-m``. Theory and solvent
    come from the program project, not from extra chain-group flags.
    """
    from chemsmart.cli.chain.chain import _chain_label
    from chemsmart.io.molecules.structure import Molecule

    chain_settings = require_chain_project(ctx)
    defaults = ctx.obj.get(CHAIN_CLI_DEFAULTS_KEY) or {}
    filename = defaults.get("filename")
    label = defaults.get("label")
    append_label = defaults.get("append_label")
    index = defaults.get("index")
    charge = defaults.get("charge")
    multiplicity = defaults.get("multiplicity")

    try:
        project_name = chain_settings.project_for(program)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc
    project_settings, job_settings = _default_program_settings(
        program, project_name
    )

    keywords = ("charge", "multiplicity")
    if charge is not None:
        job_settings.charge = charge
    if multiplicity is not None:
        job_settings.multiplicity = multiplicity

    skip_structure_load = bool(filename) and str(filename).lower().endswith(
        _TABLE_SUFFIXES
    )
    molecules = None
    molecule_indices = None
    if filename and not skip_structure_load:
        molecules = Molecule.from_filepath(
            filepath=filename, index=":", return_list=True
        )
        if molecules is None:
            raise click.UsageError(
                f"Could not obtain molecule from {filename}."
            )
        if index is not None:
            molecules, molecule_indices = (
                return_objects_and_indices_from_string_index(
                    list_of_objects=molecules, index=index
                )
            )
        if not isinstance(molecules, list):
            molecules = [molecules]
            if molecule_indices is not None and not isinstance(
                molecule_indices, list
            ):
                molecule_indices = [molecule_indices]

    resolved_label = _workflow_label(
        filename, label, append_label, job_name, _chain_label
    )

    logger.info(
        "Hydrating %s context from chain project %s (alias %s)",
        program,
        chain_settings.PROJECT_NAME,
        project_name,
    )
    ctx.obj["project_settings"] = project_settings
    ctx.obj["job_settings"] = job_settings
    ctx.obj["keywords"] = keywords
    ctx.obj["molecules"] = molecules
    ctx.obj["molecule_indices"] = molecule_indices
    ctx.obj["label"] = resolved_label
    ctx.obj["filename"] = filename
    ctx.obj[_WORKFLOW_PROGRAM_KEY] = program


def _default_program_settings(program, project_name):
    if program == "gaussian":
        from chemsmart.jobs.gaussian.settings import GaussianJobSettings
        from chemsmart.settings.gaussian import GaussianProjectSettings

        return (
            GaussianProjectSettings.from_project(project_name),
            GaussianJobSettings.default(),
        )
    if program == "orca":
        from chemsmart.jobs.orca.settings import ORCAJobSettings
        from chemsmart.settings.orca import ORCAProjectSettings

        return (
            ORCAProjectSettings.from_project(project_name),
            ORCAJobSettings.default(),
        )
    allowed = ", ".join(CHAIN_WORKFLOW_PROGRAMS)
    raise click.UsageError(
        f"Unsupported workflow program {program!r}. "
        f"Allowed programs: {allowed}."
    )


def _workflow_label(filename, label, append_label, job_name, chain_label):
    if filename:
        resolved = chain_label(filename, label, append_label)
        if label is None and append_label is None and job_name:
            return f"{resolved}_{job_name}"
        return resolved
    if append_label is not None:
        raise click.UsageError(
            "-a/--append-label requires a structure file (-f/--filename)."
        )
    if label is not None:
        return clean_label(label)
    return job_name or "output"


def _params_for(command, params):
    names = {param.name for param in command.params}
    return {name: value for name, value in params.items() if name in names}


def _import_program_workflow(program, workflow):
    module_name = f"chemsmart.cli.{program}.{workflow}"
    try:
        module = importlib.import_module(module_name)
        if workflow == "pka":
            return module.pka
        if workflow == "fukui":
            return module.fukui
        if workflow == "redox":
            return module.redox
        if workflow == "reaction":
            return module.reaction
    except (ImportError, AttributeError) as exc:
        raise click.UsageError(
            f"No {program} {workflow} command is registered."
        ) from exc
    raise click.UsageError(f"No {program} {workflow} command is registered.")


def _resolved_skip_completed(ctx):
    """Prefer COMMANDLINE ``-R``/``-S``, else the chain-group stored value."""
    current = ctx
    while current is not None:
        source = current.get_parameter_source("skip_completed")
        if source is ParameterSource.COMMANDLINE:
            return current.params["skip_completed"]
        current = current.parent
    if ctx.obj is not None:
        from chemsmart.cli.chain.chain import _SKIP_COMPLETED_KEY

        if _SKIP_COMPLETED_KEY in ctx.obj:
            return ctx.obj[_SKIP_COMPLETED_KEY]
    return True


def _params_with_resolved_skip_completed(ctx, command, extra_params=None):
    params = dict(ctx.params)
    if extra_params:
        params.update(extra_params)
    names = {param.name for param in command.params}
    if "skip_completed" in names:
        params["skip_completed"] = _resolved_skip_completed(ctx)
    return params


def _prepare_workflow_submit(ctx, program, job_name):
    program = require_workflow_program(program)
    require_chain_project(ctx)
    hydrate_program_ctx_from_chain(ctx, program, job_name=job_name)
    return program


def _invoke_program_workflow(ctx, program, workflow, extra_params=None):
    command = _import_program_workflow(program, workflow)
    params = _params_with_resolved_skip_completed(ctx, command, extra_params)
    return ctx.invoke(command, **_params_for(command, params))


def _invoke_program_subcommand(
    ctx, program, workflow, subcommand_name, extra_params=None
):
    group = _import_program_workflow(program, workflow)
    command = group.commands.get(subcommand_name)
    if command is None:
        raise click.UsageError(
            f"No {program} {workflow} {subcommand_name} command is registered."
        )
    params = _params_with_resolved_skip_completed(ctx, command, extra_params)
    return ctx.invoke(command, **_params_for(command, params))


def _store_pka_shared(ctx, kwargs):
    s_freq_cutoff, entropy_method = resolve_pka_entropy_cutoff(
        kwargs.get("cutoff_entropy_grimme"),
        kwargs.get("cutoff_entropy_truhlar"),
    )
    ctx.ensure_object(dict)
    ctx.obj["pka_shared"] = dict(
        scheme=kwargs["scheme"],
        reference=kwargs["reference"],
        reference_proton_index=kwargs["reference_proton_index"],
        reference_color_code=kwargs["reference_color_code"],
        reference_charge=kwargs["reference_charge"],
        reference_multiplicity=kwargs["reference_multiplicity"],
        reference_conjugate_base_charge=kwargs[
            "reference_conjugate_base_charge"
        ],
        reference_conjugate_base_multiplicity=kwargs[
            "reference_conjugate_base_multiplicity"
        ],
        delta_g_proton=kwargs["delta_g_proton"],
        conjugate_base_charge=kwargs["conjugate_base_charge"],
        conjugate_base_multiplicity=kwargs["conjugate_base_multiplicity"],
        solvent_model=kwargs["solvent_model"],
        solvent_id=kwargs["solvent_id"],
        temperature=kwargs["temperature"],
        concentration=kwargs["concentration"],
        pressure=kwargs["pressure"],
        cutoff_entropy_grimme=s_freq_cutoff,
        cutoff_enthalpy=kwargs["cutoff_enthalpy"],
        entropy_method=entropy_method,
        skip_completed=kwargs["skip_completed"],
    )
    ctx.obj["pka_proton_index"] = kwargs.get("proton_index")
    ctx.obj["pka_color_code"] = kwargs.get("color_code")


def _register_pka(chain_group):
    @chain_group.group("pka", cls=MyGroup, invoke_without_command=True)
    @click_job_options
    @click_pka_shared_options
    @click_pka_proton_options
    @click_workflow_program_option
    @click.pass_context
    def pka(ctx, program, **kwargs):
        """pKa workflow submit (``--program``) and output analysis.

        Submit uses the chain YAML alias for ``--program``. Analyze is
        backend-independent and does not need ``-p``, ``-f``, or ``--program``.
        """
        kwargs["skip_completed"] = _resolved_skip_completed(ctx)
        _store_pka_shared(ctx, kwargs)
        if ctx.invoked_subcommand in _ANALYZE_SUBCOMMANDS:
            return None
        program = _prepare_workflow_submit(ctx, program, job_name="pka")
        if ctx.invoked_subcommand is None:
            return _invoke_program_workflow(ctx, program, "pka")
        return None

    @pka.command("submit", cls=MyCommand)
    @click_job_options
    @click_pka_proton_options
    @click.pass_context
    def submit(ctx, skip_completed, proton_index, color_code, **kwargs):
        """Submit a single-molecule pKa calculation."""
        program = ctx.obj[_WORKFLOW_PROGRAM_KEY]
        return _invoke_program_subcommand(
            ctx,
            program,
            "pka",
            "submit",
            extra_params={
                "skip_completed": skip_completed,
                "proton_index": proton_index,
                "color_code": color_code,
                **kwargs,
            },
        )

    @pka.command("batch", cls=MyCommand)
    @click_job_options
    @click_pka_proton_options
    @click.pass_context
    def batch(ctx, skip_completed, proton_index, color_code, **kwargs):
        """Batch pKa submission from a CSV table or multi-molecule CDXML."""
        program = ctx.obj[_WORKFLOW_PROGRAM_KEY]
        return _invoke_program_subcommand(
            ctx,
            program,
            "pka",
            "batch",
            extra_params={
                "skip_completed": skip_completed,
                "proton_index": proton_index,
                "color_code": color_code,
                **kwargs,
            },
        )

    pka.add_command(pka_analyze)
    pka.add_command(pka_batch_analyze)
    return pka


def _register_fukui(chain_group):
    @chain_group.group("fukui", cls=MyGroup, invoke_without_command=True)
    @click_job_options
    @click_fukui_submit_options
    @click_workflow_program_option
    @click.pass_context
    def fukui(
        ctx,
        program,
        skip_completed,
        mode,
        radical_cation_charge,
        radical_cation_multiplicity,
        radical_anion_charge,
        radical_anion_multiplicity,
        **kwargs,
    ):
        """Fukui workflow submit (``--program``) and output analysis.

        Submit uses the chain YAML alias for ``--program``. Analyze is
        backend-independent and does not need ``-p``, ``-f``, or ``--program``.
        """
        if ctx.invoked_subcommand is not None:
            return None
        program = _prepare_workflow_submit(ctx, program, job_name="fukui")
        if program == "orca" and mode.lower() not in ORCA_FUKUI_MODES:
            allowed = ", ".join(ORCA_FUKUI_MODES)
            raise click.UsageError(
                f"ORCA Fukui does not support mode {mode!r}. "
                f"Allowed modes: {allowed}."
            )
        return _invoke_program_workflow(
            ctx,
            program,
            "fukui",
            extra_params={
                "skip_completed": skip_completed,
                "mode": mode,
                "radical_cation_charge": radical_cation_charge,
                "radical_cation_multiplicity": radical_cation_multiplicity,
                "radical_anion_charge": radical_anion_charge,
                "radical_anion_multiplicity": radical_anion_multiplicity,
                **kwargs,
            },
        )

    fukui.add_command(fukui_analyze, name="analyze")
    return fukui


def _register_redox(chain_group):
    @chain_group.group("redox", cls=MyGroup, invoke_without_command=True)
    @click_job_options
    @click_redox_shared_options
    @click_workflow_program_option
    @click.pass_context
    def redox(ctx, program, skip_completed, **kwargs):
        """Redox workflow submit (``--program``) and output analysis.

        Submit uses the chain YAML alias for ``--program``. Analyze is
        backend-independent and does not need ``-p``, ``-f``, or ``--program``.
        """
        store_redox_shared(ctx, kwargs)
        if ctx.invoked_subcommand is not None:
            return None
        program = _prepare_workflow_submit(ctx, program, job_name="redox")
        return _invoke_program_workflow(
            ctx,
            program,
            "redox",
            extra_params={"skip_completed": skip_completed, **kwargs},
        )

    redox.add_command(redox_analyze)
    return redox


def _register_reaction(chain_group):
    @chain_group.group("reaction", cls=MyGroup, invoke_without_command=True)
    @click_job_options
    @click_reaction_shared_options
    @click_workflow_program_option
    @click.pass_context
    def reaction(
        ctx,
        program,
        skip_completed,
        reactants,
        products,
        ts_guess,
        **kwargs,
    ):
        """Reaction workflow submit (``--program``). No analyze command.

        Uses ``GaussianReactionJob`` or ``ORCAReactionJob`` (QST vs NEB)
        according to ``--program``. Project settings come from the chain
        YAML alias.
        """
        store_reaction_shared(
            ctx,
            reactants=reactants,
            products=products,
            ts_guess=ts_guess,
        )
        program = _prepare_workflow_submit(ctx, program, job_name="reaction")
        if ctx.invoked_subcommand is None:
            return _invoke_program_workflow(
                ctx,
                program,
                "reaction",
                extra_params={"skip_completed": skip_completed, **kwargs},
            )
        return None

    @reaction.command("submit", cls=MyCommand)
    @click_job_options
    @click_reaction_shared_options
    @click.pass_context
    def submit(ctx, skip_completed, reactants, products, ts_guess, **kwargs):
        """Submit a single reaction workflow."""
        reactants, products, ts_guess = merge_reaction_options(
            ctx,
            reactants=reactants,
            products=products,
            ts_guess=ts_guess,
        )
        program = ctx.obj[_WORKFLOW_PROGRAM_KEY]
        return _invoke_program_subcommand(
            ctx,
            program,
            "reaction",
            "submit",
            extra_params={
                "skip_completed": skip_completed,
                "reactants": reactants,
                "products": products,
                "ts_guess": ts_guess,
                **kwargs,
            },
        )

    @reaction.command("batch", cls=MyCommand)
    @click_job_options
    @click.pass_context
    def batch(ctx, skip_completed, **kwargs):
        """Batch reaction submission from a CSV table grouped by reaction_id."""
        program = ctx.obj[_WORKFLOW_PROGRAM_KEY]
        return _invoke_program_subcommand(
            ctx,
            program,
            "reaction",
            "batch",
            extra_params={"skip_completed": skip_completed, **kwargs},
        )

    return reaction


def register_chain_workflows(chain_group):
    """Attach pKa, Fukui, redox, and reaction groups to ``chain``."""
    _register_pka(chain_group)
    _register_fukui(chain_group)
    _register_redox(chain_group)
    _register_reaction(chain_group)
