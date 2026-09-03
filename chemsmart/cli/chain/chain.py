"""CLI group for YAML chain pipelines and nested program slices."""

import functools
import logging
import os

import click

from chemsmart.cli.crest import crest
from chemsmart.cli.gaussian import gaussian
from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filename_options,
    click_job_options,
)
from chemsmart.cli.orca import orca
from chemsmart.cli.utils import CHAIN_PROJECT_SETTINGS_KEY
from chemsmart.cli.xtb import xtb
from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.chain_steps import build_chain_job, parse_chain_step
from chemsmart.settings.chain import ChainProjectSettings
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import return_objects_and_indices_from_string_index

logger = logging.getLogger(__name__)


def click_chain_options(f):
    """Project, charge, and multiplicity options for chain jobs."""

    @click.option(
        "--project",
        "-p",
        type=str,
        required=True,
        help="Chain project settings.",
    )
    @click.option(
        "-c",
        "--charge",
        type=int,
        default=None,
        help="Charge of the molecule.",
    )
    @click.option(
        "-m",
        "--multiplicity",
        type=int,
        default=None,
        help="Multiplicity of the molecule.",
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def _chain_label(filename, label, append_label):
    if label is not None and append_label is not None:
        raise click.UsageError(
            "Give only one of -l/--label and -a/--append-label."
        )
    if append_label is not None:
        stem = os.path.splitext(os.path.basename(filename))[0]
        label = f"{stem}_{append_label}"
    elif label is None:
        label = os.path.splitext(os.path.basename(filename))[0]
    return clean_label(label)


def _molecule_from_filename(filename, index):
    molecules = Molecule.from_filepath(
        filepath=filename, index=":", return_list=True
    )
    if molecules is None:
        raise click.UsageError(f"Could not obtain molecule from {filename}.")
    if index is not None:
        molecules, _ = return_objects_and_indices_from_string_index(
            list_of_objects=molecules, index=index
        )
    if not isinstance(molecules, list):
        molecules = [molecules]
    if not molecules:
        raise click.UsageError(f"No molecule selected from {filename}.")
    return molecules[-1]


def _parse_cli_steps(steps, chain_settings):
    parsed = []
    for token in steps:
        try:
            step = parse_chain_step(token)
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc
        if step.program not in ChainProjectSettings.PROGRAMS:
            allowed = ", ".join(ChainProjectSettings.PROGRAMS)
            raise click.UsageError(
                f"Unknown chain program {step.program!r}. "
                f"Allowed programs: {allowed}."
            )
        try:
            chain_settings.project_for(step.program)
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc
        parsed.append(step)
    return parsed


@click.group(cls=MyGroup, invoke_without_command=True)
@click_chain_options
@click_filename_options
@click_file_label_and_index_options
@click_job_options
@click.option(
    "-s",
    "--steps",
    "steps",
    multiple=True,
    metavar="PROGRAM:JOB",
    help="Pipeline step (repeatable), e.g. gaussian:opt.",
)
@click.pass_context
def chain(
    ctx,
    project,
    filename,
    label,
    append_label,
    index,
    charge,
    multiplicity,
    skip_completed,
    steps,
):
    """Run a YAML chain pipeline or a nested program slice.

    ``chemsmart run/sub chain -p NAME -f mol.xyz -s gaussian:opt`` runs a
    pipeline from repeatable ``-s/--steps``. Nested ``gaussian``,
    ``orca``, ``xtb``, and ``crest`` commands reuse those programs' CLIs
    with this file's project aliases, e.g.
    ``chain -p NAME gaussian pka ...``.
    """
    ctx.ensure_object(dict)
    try:
        chain_settings = ChainProjectSettings.from_project(project)
    except FileNotFoundError as exc:
        raise click.UsageError(str(exc)) from exc
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc
    ctx.obj[CHAIN_PROJECT_SETTINGS_KEY] = chain_settings

    if ctx.invoked_subcommand is not None:
        if steps:
            raise click.UsageError(
                "Cannot combine -s/--steps with a nested program command."
            )
        return None

    if not steps:
        raise click.UsageError(
            "Chain pipeline requires at least one -s/--steps."
        )
    if filename is None:
        raise click.UsageError(
            "Chain pipeline requires a structure file (-f/--filename)."
        )

    parsed_steps = _parse_cli_steps(steps, chain_settings)
    molecule = _molecule_from_filename(filename, index)
    job_label = _chain_label(filename, label, append_label)
    logger.info(
        "Building chain pipeline %s from project %s", job_label, project
    )
    try:
        return build_chain_job(
            chain_settings,
            steps=parsed_steps,
            molecule=molecule,
            label=job_label,
            charge=charge,
            multiplicity=multiplicity,
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=skip_completed,
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


chain.add_command(gaussian)
chain.add_command(orca)
chain.add_command(xtb)
chain.add_command(crest)
