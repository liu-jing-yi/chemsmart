"""Shared redox CLI options, job builder, and backend-independent analysis.

Job submission lives under Gaussian, ORCA, or chain and uses
``register_redox_cli`` / ``build_redox_job``. Output analysis is
``chemsmart run redox analyze`` (also ``chemsmart run chain redox analyze``).
"""

from __future__ import annotations

import functools
import inspect
import logging
from dataclasses import dataclass

import click

from chemsmart.cli.job import click_job_options
from chemsmart.cli.pka import (
    click_pka_thermochemistry_options,
    pka_gas_phase_data,
    pka_solvent_scf_energy,
    resolve_pka_entropy_cutoff,
)
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.constants import FARADAY, energy_conversion
from chemsmart.utils.utils import check_charge_and_multiplicity

logger = logging.getLogger(__name__)

_REGISTRY: dict[str, RedoxReference] = {}


@dataclass(frozen=True)
class RedoxReference:
    """One experimental or computational redox reference couple.

    Attributes:
        name: Registry key, e.g. ``fc_fc+``.
        E_ref_V: Reference reduction potential in volts.
        n_electrons: Electrons transferred in the couple.
        scale: Reporting scale, e.g. ``Fc/Fc+`` or ``SHE``.
        couple_label: Chemical couple, e.g. ``Fc/Fc+``.
        ox_file: Optional path to the oxidized geometry.
        red_file: Optional path to the reduced geometry.
        ox_charge: Optional charge of the oxidized species.
        ox_multiplicity: Optional multiplicity of the oxidized species.
        red_charge: Optional charge of the reduced species.
        red_multiplicity: Optional multiplicity of the reduced species.
        solvent: Optional solvent label associated with ``E_ref_V``.
    """

    name: str
    E_ref_V: float
    n_electrons: int
    scale: str
    couple_label: str
    ox_file: str | None = None
    red_file: str | None = None
    ox_charge: int | None = None
    ox_multiplicity: int | None = None
    red_charge: int | None = None
    red_multiplicity: int | None = None
    solvent: str | None = None

    def __post_init__(self):
        if not self.name:
            raise ValueError("RedoxReference.name is required.")
        if self.n_electrons < 1:
            raise ValueError("n_electrons must be a positive integer.")


def register_redox_reference(reference):
    """Register a :class:`RedoxReference` under ``reference.name``.

    Args:
        reference (RedoxReference): Couple to add to the registry.

    Raises:
        TypeError: If ``reference`` is not a :class:`RedoxReference`.
        ValueError: If ``reference.name`` is already registered.
    """
    if not isinstance(reference, RedoxReference):
        raise TypeError(
            "reference must be a RedoxReference, "
            f"got {type(reference).__name__}."
        )
    if reference.name in _REGISTRY:
        raise ValueError(
            f"Redox reference {reference.name!r} is already registered."
        )
    _REGISTRY[reference.name] = reference


def get_redox_reference(name):
    """Return the registered :class:`RedoxReference` for ``name``.

    Args:
        name (str): Registry key.

    Returns:
        RedoxReference: The registered couple.

    Raises:
        ValueError: If ``name`` is not registered.
    """
    reference = _REGISTRY.get(name)
    if reference is None:
        available = ", ".join(sorted(_REGISTRY)) or "(none)"
        raise ValueError(
            f"Unknown redox reference {name!r}. Available: {available}."
        )
    return reference


def list_redox_references():
    """Return registered couples sorted by name."""
    return tuple(_REGISTRY[name] for name in sorted(_REGISTRY))


def resolve_redox_reference(reference=None):
    """Return a :class:`RedoxReference` from an object, name, or default.

    Args:
        reference: A :class:`RedoxReference`, registry name, or ``None``
            (defaults to ``fc_fc+``).

    Returns:
        RedoxReference: Resolved couple.

    Raises:
        TypeError: If ``reference`` is not a supported type.
    """
    if reference is None:
        return get_redox_reference("fc_fc+")
    if isinstance(reference, RedoxReference):
        return reference
    if isinstance(reference, str):
        return get_redox_reference(reference)
    raise TypeError(
        "reference must be a RedoxReference, registry name, or None, "
        f"got {type(reference).__name__}."
    )


def _species_terms(gas_file, solv_file, thermo_kwargs):
    """Return gas, correction, solvent, and solution G terms in Hartree."""
    e_gas_au, g_corr_au = pka_gas_phase_data(gas_file, **thermo_kwargs)
    e_solv_au = pka_solvent_scf_energy(solv_file)
    return {
        "E_gas_au": e_gas_au,
        "G_corr_au": g_corr_au,
        "E_solv_au": e_solv_au,
        "G_soln_au": e_solv_au + g_corr_au,
    }


def _require_files(**files):
    missing = [name for name, path in files.items() if not path]
    if missing:
        raise ValueError(
            "Missing required files for redox exchange analysis: "
            + ", ".join(missing)
        )


def compute_redox_potential(
    ox_gas_file,
    red_gas_file,
    ref_ox_gas_file,
    ref_red_gas_file,
    ox_solv_file,
    red_solv_file,
    ref_ox_solv_file,
    ref_red_solv_file,
    reference=None,
    n_electrons=None,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    entropy_method="grimme",
):
    """Compute the exchange redox potential from output files.

    Dual-level solution free energies follow the pKa helpers:
    ``G_soln = E_solv + (qh-G − E_gas)``.

    Reaction: ``Ox + Ref_red → Red + Ref_ox``

    ``ΔG_exchange = G(Red) + G(Ref_ox) − G(Ox) − G(Ref_red)``

    ``E_target = E_ref − ΔG_exchange / (n F)``

    Args:
        ox_gas_file: Gas-phase opt+freq output for the oxidized target.
        red_gas_file: Gas-phase opt+freq output for the reduced target.
        ref_ox_gas_file: Gas-phase opt+freq output for the oxidized reference.
        ref_red_gas_file: Gas-phase opt+freq output for the reduced reference.
        ox_solv_file: Solution-phase SP output for the oxidized target.
        red_solv_file: Solution-phase SP output for the reduced target.
        ref_ox_solv_file: Solution-phase SP output for the oxidized reference.
        ref_red_solv_file: Solution-phase SP output for the reduced reference.
        reference: :class:`RedoxReference`, registry name, or ``None``
            (defaults to ``fc_fc+``).
        n_electrons: Electrons transferred. Defaults to the reference couple.
            Must match the reference couple when given.
        temperature: Temperature in Kelvin for quasi-harmonic G.
        concentration: Concentration in mol/L.
        pressure: Pressure in atm.
        cutoff_entropy_grimme: Entropy cutoff in cm⁻¹ (Grimme qRRHO).
        cutoff_enthalpy: Enthalpy cutoff in cm⁻¹.
        entropy_method: Quasi-RRHO entropy method (``grimme`` or ``truhlar``).

    Returns:
        dict: Exchange ΔG (au, kcal/mol, eV, J/mol), ``E_ref_V``,
        ``E_target_V``, reference metadata, and per-species G terms.
    """
    _require_files(
        ox_gas_file=ox_gas_file,
        red_gas_file=red_gas_file,
        ref_ox_gas_file=ref_ox_gas_file,
        ref_red_gas_file=ref_red_gas_file,
        ox_solv_file=ox_solv_file,
        red_solv_file=red_solv_file,
        ref_ox_solv_file=ref_ox_solv_file,
        ref_red_solv_file=ref_red_solv_file,
    )
    reference = resolve_redox_reference(reference)
    if n_electrons is None:
        n_electrons = reference.n_electrons
    elif n_electrons != reference.n_electrons:
        raise ValueError(
            f"n_electrons ({n_electrons}) must match the reference couple "
            f"{reference.name!r} (n={reference.n_electrons})."
        )
    if n_electrons < 1:
        raise ValueError("n_electrons must be a positive integer.")

    thermo_kwargs = dict(
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
        entropy_method=entropy_method,
    )
    ox = _species_terms(ox_gas_file, ox_solv_file, thermo_kwargs)
    red = _species_terms(red_gas_file, red_solv_file, thermo_kwargs)
    ref_ox = _species_terms(ref_ox_gas_file, ref_ox_solv_file, thermo_kwargs)
    ref_red = _species_terms(
        ref_red_gas_file, ref_red_solv_file, thermo_kwargs
    )

    delta_g_exchange_au = (
        red["G_soln_au"]
        + ref_ox["G_soln_au"]
        - ox["G_soln_au"]
        - ref_red["G_soln_au"]
    )
    delta_g_exchange_j_mol = energy_conversion(
        "hartree", "j/mol", delta_g_exchange_au
    )
    delta_g_exchange_kcal_mol = energy_conversion(
        "hartree", "kcal/mol", delta_g_exchange_au
    )
    delta_g_exchange_ev = energy_conversion(
        "hartree", "ev", delta_g_exchange_au
    )
    e_target_v = reference.E_ref_V - (
        delta_g_exchange_j_mol / (n_electrons * FARADAY)
    )
    return {
        "delta_G_exchange_au": delta_g_exchange_au,
        "delta_G_exchange_kcal_mol": delta_g_exchange_kcal_mol,
        "delta_G_exchange_eV": delta_g_exchange_ev,
        "delta_G_exchange_J_mol": delta_g_exchange_j_mol,
        "E_ref_V": reference.E_ref_V,
        "E_target_V": e_target_v,
        "scale": reference.scale,
        "couple_label": reference.couple_label,
        "reference_name": reference.name,
        "n_electrons": n_electrons,
        "temperature": temperature,
        "G_soln_ox_au": ox["G_soln_au"],
        "G_soln_red_au": red["G_soln_au"],
        "G_soln_ref_ox_au": ref_ox["G_soln_au"],
        "G_soln_ref_red_au": ref_red["G_soln_au"],
        "E_solv_ox_au": ox["E_solv_au"],
        "E_solv_red_au": red["E_solv_au"],
        "E_solv_ref_ox_au": ref_ox["E_solv_au"],
        "E_solv_ref_red_au": ref_red["E_solv_au"],
        "G_corr_ox_au": ox["G_corr_au"],
        "G_corr_red_au": red["G_corr_au"],
        "G_corr_ref_ox_au": ref_ox["G_corr_au"],
        "G_corr_ref_red_au": ref_red["G_corr_au"],
        "E_gas_ox_au": ox["E_gas_au"],
        "E_gas_red_au": red["E_gas_au"],
        "E_gas_ref_ox_au": ref_ox["E_gas_au"],
        "E_gas_ref_red_au": ref_red["E_gas_au"],
    }


def format_redox_summary(result):
    """Return a formatted summary of an exchange redox calculation."""
    lines = [
        "=" * 78,
        "Redox Potential - Dual-level Exchange Scheme",
        "=" * 78,
        "Reaction: Ox + Ref_red → Red + Ref_ox",
        f"Reference: {result['reference_name']} "
        f"({result['couple_label']}, {result['scale']})",
        f"n = {result['n_electrons']}",
        f"Temperature: {result['temperature']} K",
        "",
        "Method:",
        "  G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)",
        "  G_soln = E_solv + G_corr  (solution free energy)",
        "  ΔG_exchange = G(Red) + G(Ref_ox) − G(Ox) − G(Ref_red)",
        "  E_target = E_ref − ΔG_exchange / (n F)",
        "-" * 78,
        "",
        "Gas-Phase Electronic Energies (E_gas, au):",
        f"  Ox:      {result['E_gas_ox_au']:.10f}",
        f"  Red:     {result['E_gas_red_au']:.10f}",
        f"  Ref_ox:  {result['E_gas_ref_ox_au']:.10f}",
        f"  Ref_red: {result['E_gas_ref_red_au']:.10f}",
        "",
        "Thermal Corrections (G_corr = qh-G - E_gas, au):",
        f"  Ox:      {result['G_corr_ox_au']:.10f}",
        f"  Red:     {result['G_corr_red_au']:.10f}",
        f"  Ref_ox:  {result['G_corr_ref_ox_au']:.10f}",
        f"  Ref_red: {result['G_corr_ref_red_au']:.10f}",
        "",
        "Solvent Single-Point Energies (E_solv, au):",
        f"  Ox:      {result['E_solv_ox_au']:.10f}",
        f"  Red:     {result['E_solv_red_au']:.10f}",
        f"  Ref_ox:  {result['E_solv_ref_ox_au']:.10f}",
        f"  Ref_red: {result['E_solv_ref_red_au']:.10f}",
        "",
        "Solution Free Energies (G_soln = E_solv + G_corr, au):",
        f"  Ox:      {result['G_soln_ox_au']:.10f}",
        f"  Red:     {result['G_soln_red_au']:.10f}",
        f"  Ref_ox:  {result['G_soln_ref_ox_au']:.10f}",
        f"  Ref_red: {result['G_soln_ref_red_au']:.10f}",
        "-" * 78,
        "",
        "Redox Potential:",
        f"  ΔG_exchange = {result['delta_G_exchange_au']:.10f} au",
        f"              = {result['delta_G_exchange_kcal_mol']:.4f} kcal/mol",
        f"              = {result['delta_G_exchange_eV']:.4f} eV",
        f"              = {result['delta_G_exchange_J_mol']:.4f} J/mol",
        f"  E_ref = {result['E_ref_V']:.4f} V ({result['scale']})",
        "",
        f"  *** E_target = {result['E_target_V']:.4f} V "
        f"({result['scale']}) ***",
        "=" * 78,
    ]
    return "\n".join(lines)


register_redox_reference(
    RedoxReference(
        name="fc_fc+",
        E_ref_V=0.0,
        n_electrons=1,
        scale="Fc/Fc+",
        couple_label="Fc/Fc+",
    )
)


def click_redox_reference_options(f):
    """Reference couple and electron-count options for submit and analyze."""

    @click.option(
        "-n",
        "--n-electrons",
        type=int,
        default=None,
        help=(
            "Electrons transferred. Defaults to the reference couple. "
            "Must match the reference couple when given."
        ),
    )
    @click.option(
        "-r",
        "--reference",
        type=str,
        default="fc_fc+",
        show_default=True,
        help="Registry name of the reference couple.",
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_redox_submit_structure_options(f):
    """Target/reference geometry and charge options for redox submit."""

    @click.option(
        "-rrm",
        "--ref-red-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the reduced reference. Defaults to the "
            "registry value or the reduced-form multiplicity of Ref_ox."
        ),
    )
    @click.option(
        "-rrc",
        "--ref-red-charge",
        type=int,
        default=None,
        help=(
            "Charge of the reduced reference. Defaults to the registry "
            "value or ox − n of the oxidized reference."
        ),
    )
    @click.option(
        "-rom",
        "--ref-ox-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the oxidized reference. Defaults to the "
            "registry value or the loaded geometry."
        ),
    )
    @click.option(
        "-roc",
        "--ref-ox-charge",
        type=int,
        default=None,
        help=(
            "Charge of the oxidized reference. Defaults to the registry "
            "value or the loaded geometry."
        ),
    )
    @click.option(
        "-rr",
        "--ref-red",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help=(
            "Reduced reference geometry. Required unless the registered "
            "couple provides red_file."
        ),
    )
    @click.option(
        "-ro",
        "--ref-ox",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help=(
            "Oxidized reference geometry. Required unless the registered "
            "couple provides ox_file."
        ),
    )
    @click.option(
        "-rdm",
        "--red-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the reduced target. Defaults to the "
            "reduced-form multiplicity of the oxidized target."
        ),
    )
    @click.option(
        "-rdc",
        "--red-charge",
        type=int,
        default=None,
        help=(
            "Charge of the reduced target. Defaults to ox − n from the "
            "parent oxidized charge."
        ),
    )
    @click.option(
        "-rd",
        "--red",
        type=click.Path(exists=True, dir_okay=False),
        default=None,
        help=(
            "Reduced target geometry. Defaults to the oxidized structure "
            "from parent -f with charge ox − n."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_redox_shared_options(f):
    """Submit options shared by Gaussian, ORCA, and chain redox commands."""
    f = click_redox_submit_structure_options(f)
    f = click_redox_reference_options(f)
    return click_pka_thermochemistry_options(f)


def click_redox_analyze_options(f):
    """Output-file options for exchange redox analysis."""
    f = click.option(
        "-rrs",
        "--ref-red-solv",
        "ref_red_solv_file",
        type=click.Path(exists=True, dir_okay=False),
        required=True,
        help="Solution-phase SP output for the reduced reference.",
    )(f)
    f = click.option(
        "-ros",
        "--ref-ox-solv",
        "ref_ox_solv_file",
        type=click.Path(exists=True, dir_okay=False),
        required=True,
        help="Solution-phase SP output for the oxidized reference.",
    )(f)
    f = click.option(
        "-rds",
        "--red-solv",
        "red_solv_file",
        type=click.Path(exists=True, dir_okay=False),
        required=True,
        help="Solution-phase SP output for the reduced target.",
    )(f)
    f = click.option(
        "-oxs",
        "--ox-solv",
        "ox_solv_file",
        type=click.Path(exists=True, dir_okay=False),
        required=True,
        help="Solution-phase SP output for the oxidized target.",
    )(f)
    f = click.option(
        "-rrg",
        "--ref-red-gas",
        "ref_red_gas_file",
        type=click.Path(exists=True, dir_okay=False),
        required=True,
        help="Gas-phase opt+freq output for the reduced reference.",
    )(f)
    f = click.option(
        "-rog",
        "--ref-ox-gas",
        "ref_ox_gas_file",
        type=click.Path(exists=True, dir_okay=False),
        required=True,
        help="Gas-phase opt+freq output for the oxidized reference.",
    )(f)
    f = click.option(
        "-rdg",
        "--red-gas",
        "red_gas_file",
        type=click.Path(exists=True, dir_okay=False),
        required=True,
        help="Gas-phase opt+freq output for the reduced target.",
    )(f)
    f = click.option(
        "-oxg",
        "--ox-gas",
        "ox_gas_file",
        type=click.Path(exists=True, dir_okay=False),
        required=True,
        help="Gas-phase opt+freq output for the oxidized target.",
    )(f)
    return f


def store_redox_shared(ctx, kwargs):
    """Record redox CLI options on ``ctx.obj`` for submit and analyze."""
    s_freq_cutoff, entropy_method = resolve_pka_entropy_cutoff(
        kwargs.get("cutoff_entropy_grimme"),
        kwargs.get("cutoff_entropy_truhlar"),
    )
    ctx.ensure_object(dict)
    ctx.obj["redox_shared"] = dict(
        reference=kwargs.get("reference") or "fc_fc+",
        n_electrons=kwargs.get("n_electrons"),
        red_file=kwargs.get("red"),
        red_charge=kwargs.get("red_charge"),
        red_multiplicity=kwargs.get("red_multiplicity"),
        ref_ox_file=kwargs.get("ref_ox"),
        ref_red_file=kwargs.get("ref_red"),
        ref_ox_charge=kwargs.get("ref_ox_charge"),
        ref_ox_multiplicity=kwargs.get("ref_ox_multiplicity"),
        ref_red_charge=kwargs.get("ref_red_charge"),
        ref_red_multiplicity=kwargs.get("ref_red_multiplicity"),
        temperature=kwargs.get("temperature", 298.15),
        concentration=kwargs.get("concentration", 1.0),
        pressure=kwargs.get("pressure", 1.0),
        cutoff_entropy_grimme=s_freq_cutoff,
        cutoff_enthalpy=kwargs.get("cutoff_enthalpy", 100.0),
        entropy_method=entropy_method,
    )


def _base_settings_kwargs(opt_settings, settings_cls):
    from chemsmart.jobs.gaussian.settings import GaussianJobSettings
    from chemsmart.jobs.orca.settings import ORCAJobSettings

    if issubclass(settings_cls, GaussianJobSettings):
        base_cls = GaussianJobSettings
    elif issubclass(settings_cls, ORCAJobSettings):
        base_cls = ORCAJobSettings
    else:
        raise TypeError(
            "settings_cls must be a Gaussian or ORCA job settings class, "
            f"got {settings_cls!r}."
        )
    params = {
        name
        for name, param in inspect.signature(
            base_cls.__init__
        ).parameters.items()
        if name != "self"
        and param.kind
        not in (
            inspect.Parameter.VAR_POSITIONAL,
            inspect.Parameter.VAR_KEYWORD,
        )
    }
    return {
        key: value
        for key, value in vars(opt_settings).items()
        if key in params and value is not None
    }


def _first_non_none(*values):
    for value in values:
        if value is not None:
            return value
    return None


def build_redox_job(ctx, job_cls, settings_cls, skip_completed, **kwargs):
    """Create one redox job from parent ``-f`` and shared redox options."""
    if kwargs:
        store_redox_shared(ctx, kwargs)
    shared = ctx.obj["redox_shared"]

    molecules = ctx.obj.get("molecules")
    if not molecules:
        raise click.UsageError(
            "Redox submit requires a structure file via parent -f/--filename."
        )

    project_settings = ctx.obj["project_settings"]
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]
    opt_settings = project_settings.opt_settings().merge(
        job_settings, keywords=keywords
    )
    sp_settings = project_settings.sp_settings()
    settings_kwargs = _base_settings_kwargs(opt_settings, settings_cls)
    settings_kwargs["solvent_model"] = _first_non_none(
        settings_kwargs.get("solvent_model"),
        opt_settings.solvent_model,
        sp_settings.solvent_model,
    )
    settings_kwargs["solvent_id"] = _first_non_none(
        settings_kwargs.get("solvent_id"),
        opt_settings.solvent_id,
        sp_settings.solvent_id,
    )
    settings_kwargs.update(
        n_electrons=shared["n_electrons"],
        red_file=shared["red_file"],
        red_charge=shared["red_charge"],
        red_multiplicity=shared["red_multiplicity"],
        reference=shared["reference"],
        ref_ox_file=shared["ref_ox_file"],
        ref_red_file=shared["ref_red_file"],
        ref_ox_charge=shared["ref_ox_charge"],
        ref_ox_multiplicity=shared["ref_ox_multiplicity"],
        ref_red_charge=shared["ref_red_charge"],
        ref_red_multiplicity=shared["ref_red_multiplicity"],
        temperature=shared["temperature"],
        concentration=shared["concentration"],
        pressure=shared["pressure"],
        cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
        cutoff_enthalpy=shared["cutoff_enthalpy"],
    )
    settings_kwargs = {
        key: value
        for key, value in settings_kwargs.items()
        if value is not None
    }

    label = ctx.obj["label"]
    logger.info("Creating %s job %s", job_cls.__name__, label)
    try:
        settings = settings_cls(**settings_kwargs)
        check_charge_and_multiplicity(settings)
        return job_cls(
            molecule=molecules[-1],
            settings=settings,
            label=label,
            jobrunner=ctx.obj["jobrunner"],
            skip_completed=skip_completed,
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


def register_redox_cli(parent_group, job_cls, settings_cls):
    """Attach ``redox`` submit to a Gaussian or ORCA Click group."""

    @parent_group.command("redox", cls=MyCommand)
    @click_job_options
    @click_redox_shared_options
    @click.pass_context
    def redox(ctx, skip_completed, **kwargs):
        """Submit a dual-level exchange redox calculation.

        Oxidized target comes from parent ``-f``; reduced target uses the
        same geometry (or ``--red``) with charge ``ox − n``. Reference
        geometries come from the registry and/or ``--ref-ox`` / ``--ref-red``.

        Analyze completed outputs with ``chemsmart run redox analyze``.
        """
        return build_redox_job(
            ctx, job_cls, settings_cls, skip_completed, **kwargs
        )

    return redox


@click.group(name="redox", cls=MyGroup)
@click_pka_thermochemistry_options
@click_redox_reference_options
@click.pass_context
def redox(
    ctx,
    temperature,
    concentration,
    pressure,
    cutoff_entropy_grimme,
    cutoff_entropy_truhlar,
    cutoff_enthalpy,
    reference,
    n_electrons,
):
    """Backend-independent redox output analysis.

    Submit jobs with ``chemsmart run/sub gaussian|orca … redox`` or
    ``chemsmart run/sub chain … redox --program {gaussian,orca}``.
    """
    store_redox_shared(
        ctx,
        dict(
            temperature=temperature,
            concentration=concentration,
            pressure=pressure,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_entropy_truhlar=cutoff_entropy_truhlar,
            cutoff_enthalpy=cutoff_enthalpy,
            reference=reference,
            n_electrons=n_electrons,
        ),
    )


@redox.command("analyze", cls=MyCommand)
@click_redox_analyze_options
@click.pass_context
def analyze(
    ctx,
    ox_gas_file,
    red_gas_file,
    ref_ox_gas_file,
    ref_red_gas_file,
    ox_solv_file,
    red_solv_file,
    ref_ox_solv_file,
    ref_red_solv_file,
):
    """Compute the exchange redox potential from existing output files.

    Dual-level solution free energies: ``G_soln = E_solv + (qh-G − E_gas)``.

    Reaction: ``Ox + Ref_red → Red + Ref_ox``

    ``E_target = E_ref − ΔG_exchange / (n F)``
    """
    shared = {}
    if ctx.obj is not None:
        shared = ctx.obj.get("redox_shared") or {}
    try:
        result = compute_redox_potential(
            ox_gas_file=ox_gas_file,
            red_gas_file=red_gas_file,
            ref_ox_gas_file=ref_ox_gas_file,
            ref_red_gas_file=ref_red_gas_file,
            ox_solv_file=ox_solv_file,
            red_solv_file=red_solv_file,
            ref_ox_solv_file=ref_ox_solv_file,
            ref_red_solv_file=ref_red_solv_file,
            reference=shared.get("reference"),
            n_electrons=shared.get("n_electrons"),
            temperature=shared.get("temperature", 298.15),
            concentration=shared.get("concentration", 1.0),
            pressure=shared.get("pressure", 1.0),
            cutoff_entropy_grimme=shared.get("cutoff_entropy_grimme", 100.0),
            cutoff_enthalpy=shared.get("cutoff_enthalpy", 100.0),
            entropy_method=shared.get("entropy_method") or "grimme",
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc
    click.echo(format_redox_summary(result))
    return None
