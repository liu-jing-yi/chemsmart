"""
Shared pKa Click commands, submission, and site preflight.

``chemsmart run pka`` computes pKa from existing output files.
Gaussian and ORCA ``pka submit`` and ``pka batch`` use the drivers here.

Subcommands
-----------
analyze        Compute pKa from 8 existing output files.
batch-analyze  Batch pKa from a table of output file paths.
"""

import copy
import functools
import logging
import os
from pathlib import Path
from typing import NamedTuple

import click

from chemsmart.analysis.pka import compute_pka, print_pka_summary
from chemsmart.cli.thermochemistry.thermochemistry import (
    resolve_entropy_cutoff,
    thermochemistry_cutoff_options,
    thermochemistry_temp_pressure_conc_options,
)
from chemsmart.io.file import PKaCDXFile
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)

CHEMDRAW_MOLECULAR_FRAGMENT_WARNING = (
    "This file contains more than one ChemDraw molecular fragment. "
    "Salts, counterions, explicit solvent, catalyst/ligand pairs, and "
    "other disconnected components may not be independent molecules."
)


def resolve_pka_entropy_cutoff(cutoff_entropy_grimme, cutoff_entropy_truhlar):
    """Resolve pKa entropy cutoff; default to Grimme 100 cm⁻¹ when unset."""
    s_freq_cutoff, entropy_method = resolve_entropy_cutoff(
        cutoff_entropy_grimme, cutoff_entropy_truhlar
    )
    if s_freq_cutoff is None:
        return 100.0, "grimme"
    return s_freq_cutoff, entropy_method


def click_pka_pkb_options(f):
    """pKb conversion options for pKa submission and analysis."""
    f = click.option(
        "--pks",
        type=float,
        default=None,
        help=(
            "Solvent autoprotolysis constant. If --pkb is set and --pks "
            "is omitted, 14.0 is used. Supplying --pks also enables "
            "pKb reporting."
        ),
    )(f)
    f = click.option(
        "--pkb",
        is_flag=True,
        default=False,
        help=(
            "Submit: protonate the input free base (B → BH+) then run pKa "
            "jobs. Analysis: also print pKb = pKs - pKa. If --pks is "
            "omitted, 14.0 is used."
        ),
    )(f)
    return f


def click_pka_thermochemistry_options(f):
    """Thermochemistry options reused by pKa submission and analysis."""
    f = thermochemistry_temp_pressure_conc_options(
        f,
        temperature_required=False,
        temperature_default=298.15,
        concentration_default=1.0,
        pressure_default=1.0,
        concentration_short="-c",
    )
    return thermochemistry_cutoff_options(
        f,
        enthalpy_default=100.0,
    )


def resolve_pka_sampling_options(sampling=False, num_conformers=1):
    """Return ``(sampling, num_conformers)`` for pKa submission.

    ``-N/--num-conformers`` greater than 1 requires ``--sampling``.
    """
    sampling = bool(sampling)
    if num_conformers is None:
        num_conformers = 1
    if num_conformers > 1 and not sampling:
        raise click.UsageError(
            "-N/--num-conformers requires --sampling when greater than 1."
        )
    return sampling, num_conformers


def click_pka_shared_options(f):
    f = click_pka_thermochemistry_options(f)
    f = click_pka_pkb_options(f)

    @click.option(
        "-s",
        "--scheme",
        "scheme",
        type=click.Choice(["direct", "proton exchange"]),
        default="proton exchange",
        help=(
            "Thermodynamic cycle type. 'proton exchange' uses a reference acid (default). 'direct' uses G_soln(H+) in water."
        ),
    )
    @click.option(
        "-r",
        "--reference",
        type=click.Path(exists=True),
        default=None,
        help=(
            "Path to geometry file for reference acid (HRef) for proton exchange cycle."
        ),
    )
    @click.option(
        "-rpi",
        "--reference-proton-index",
        type=int,
        default=None,
        help=(
            "1-based index of the proton to remove from reference acid (HRef). Required when --reference is provided."
        ),
    )
    @click.option(
        "-rcc",
        "--reference-color-code",
        type=int,
        default=None,
        help=(
            "CDXML colour-table index identifying the proton in the reference acid file."
        ),
    )
    @click.option(
        "-rc",
        "--reference-charge",
        type=int,
        default=None,
        help="Charge of the reference acid (HRef).",
    )
    @click.option(
        "-rm",
        "--reference-multiplicity",
        type=int,
        default=None,
        help="Multiplicity of the reference acid (HRef).",
    )
    @click.option(
        "--reference-conjugate-base-charge",
        type=int,
        default=None,
        help=(
            "Charge of the reference conjugate base (Ref-). Defaults to (reference_charge - 1)."
        ),
    )
    @click.option(
        "--reference-conjugate-base-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the reference conjugate base (Ref-). Defaults to reference_multiplicity."
        ),
    )
    @click.option(
        "-dG",
        "--delta-g-proton",
        type=float,
        default=None,
        help=(
            "G_soln(H+) in kcal/mol for the direct cycle. If omitted, a "
            "T-dependent aqueous default is computed from Kelly, Cramer, "
            "and Truhlar ΔG_solv(H+) = -265.9 kcal/mol."
        ),
    )
    @click.option(
        "--conjugate-base-charge",
        type=int,
        default=None,
        help="Charge of the conjugate base (A-). Defaults to (charge - 1).",
    )
    @click.option(
        "--conjugate-base-multiplicity",
        type=int,
        default=None,
        help=(
            "Multiplicity of the conjugate base (A-). Defaults to multiplicity."
        ),
    )
    @click.option(
        "-sm",
        "--solvent-model",
        type=str,
        default=None,
        help="Solvation model for solution phase SP.",
    )
    @click.option(
        "-si",
        "--solvent-id",
        type=str,
        default=None,
        help="Solvent ID for solution phase SP (default: project setting or water).",
    )
    @click.option(
        "--sampling/--no-sampling",
        default=False,
        type=bool,
        help=(
            "Enable CREST conformational sampling before DFT for HA, A-, "
            "and, when set, the reference acid. Off by default."
        ),
    )
    @click.option(
        "-N",
        "--num-conformers",
        type=click.IntRange(min=1),
        default=1,
        show_default=True,
        help=(
            "Number of lowest-energy CREST conformers per sampled "
            "species. Distinct from -n/--num-cores. Each conformer gets "
            "one gas-phase opt+freq job and one matching solvent "
            "single-point. N = 1 keeps legacy filenames without _c1. "
            "Values greater than 1 require --sampling."
        ),
    )
    @click.option(
        "--preview",
        is_flag=True,
        default=False,
        help=(
            "Print a preflight table for the input structure and stop "
            "before creating or submitting Gaussian, ORCA, or CREST jobs."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def is_pka_cdxml_input(filename):
    """Return True when *filename* is a ChemDraw CDX/CDXML structure file."""
    return bool(filename) and str(filename).lower().endswith(
        (".cdx", ".cdxml")
    )


# Conservative unique-site SMARTS. Acid patterns map to the ionizable
# hydrogen; base patterns map to the heavy atom to protonate.
PKA_ACID_SMARTS = (
    "[CX3](=O)[OX2H1][#1]",
    "[c][OX2H1][#1]",
    "[SX2H1][#1]",
    "[NX4;+1][#1]",
    "[nH;+1][#1]",
)
PKB_BASE_SMARTS = (
    "[NX3;H2,H1;!$(NC=[O,S])]",
    "[NX3;H0;!$(NC=[O,S]);!$(N=*)]",
    "[nX2;H0]",
)
_IONIZABLE_SITE_SMARTS = {
    "acid": PKA_ACID_SMARTS,
    "base": PKB_BASE_SMARTS,
}
_IONIZABLE_SITE_ATOMIC_NUM = {
    "acid": 1,
    "base": 7,
}


@functools.lru_cache(maxsize=None)
def _ionizable_smarts_mol(smarts):
    from rdkit import Chem

    pattern = Chem.MolFromSmarts(smarts)
    if pattern is None:
        raise ValueError(f"Invalid SMARTS pattern: {smarts}")
    return pattern


def _prepare_rdkit_mol_for_ionizable_site(molecule):
    """Return an RDKit mol with aromaticity and N charges suitable for SMARTS."""
    from rdkit import Chem

    rdkit_mol = molecule.to_rdkit()
    for atom in rdkit_mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            atom.SetIsAromatic(False)
            for bond in atom.GetBonds():
                bond.SetIsAromatic(False)
                if bond.GetBondType() == Chem.BondType.AROMATIC:
                    bond.SetBondType(Chem.BondType.SINGLE)
        elif atom.GetAtomicNum() == 7 and atom.GetDegree() == 4:
            atom.SetFormalCharge(1)
    rdkit_mol.UpdatePropertyCache(strict=False)
    sanitize_ops = (
        Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_ADJUSTHS
    )
    try:
        Chem.SanitizeMol(rdkit_mol, sanitizeOps=sanitize_ops)
    except Chem.MolSanitizeException:
        Chem.SetAromaticity(rdkit_mol)

    overall_charge = 0 if molecule.charge is None else int(molecule.charge)
    remaining = overall_charge - Chem.GetFormalCharge(rdkit_mol)
    if remaining > 0:
        candidates = [
            atom
            for atom in rdkit_mol.GetAtoms()
            if atom.GetAtomicNum() == 7
            and atom.GetFormalCharge() == 0
            and atom.GetIsAromatic()
            and atom.GetDegree() == 3
            and any(n.GetAtomicNum() == 1 for n in atom.GetNeighbors())
        ]
        if len(candidates) == remaining:
            for atom in candidates:
                atom.SetFormalCharge(1)
    return rdkit_mol


def resolve_ionizable_site(molecule, mode="acid"):
    """Return the unique 1-based ionizable site from acid or base SMARTS.

    Acid mode returns the acidic hydrogen index. Base mode returns the
    heavy-atom index to protonate. Matching uses :meth:`Molecule.to_rdkit`.

    Args:
        molecule: Structure to search.
        mode: ``"acid"`` or ``"base"``.

    Returns:
        int: 1-based atom index of the unique matching site.

    Raises:
        ValueError: If *mode* is invalid, or SMARTS finds 0 or more than one
            site. Specify ``-pi/--proton-index`` or ``-cc/--color-code``.
    """
    if mode not in _IONIZABLE_SITE_SMARTS:
        raise ValueError(f"mode must be 'acid' or 'base', got {mode!r}")

    rdkit_mol = _prepare_rdkit_mol_for_ionizable_site(molecule)
    atomic_num = _IONIZABLE_SITE_ATOMIC_NUM[mode]
    sites = set()
    for smarts in _IONIZABLE_SITE_SMARTS[mode]:
        for match in rdkit_mol.GetSubstructMatches(
            _ionizable_smarts_mol(smarts)
        ):
            sites.update(
                idx
                for idx in match
                if rdkit_mol.GetAtomWithIdx(idx).GetAtomicNum() == atomic_num
            )

    if len(sites) != 1:
        kind = "acidic hydrogen" if mode == "acid" else "basic atom"
        raise ValueError(
            f"Could not uniquely identify the {kind} "
            f"({len(sites)} SMARTS matches). "
            "Specify -pi/--proton-index or -cc/--color-code."
        )
    return next(iter(sites)) + 1


def pka_submit_site_mode(pkb=False):
    """Return ``"base"`` when submitting with ``--pkb``, else ``"acid"``."""
    return "base" if pkb else "acid"


def _cdxml_has_no_colour_markup(exc):
    """Return True when CDXML colour detection failed for lack of markup."""
    return "share the same colour" in str(exc)


def prepare_pkb_submit_molecule(molecule, site_index, opt_settings):
    """Protonate the free base and increment settings charge by 1.

    *site_index* is the 1-based heavy atom to protonate. *opt_settings.charge*
    is the charge of the input free base when set. The returned molecule is
    BH+ and the returned index is the added hydrogen, which existing pKa
    jobs remove to recover B.
    """
    import copy

    from chemsmart.io.molecules.structure import PKaMolecule

    updated = copy.copy(opt_settings)
    if updated.multiplicity is None and molecule.multiplicity is not None:
        updated.multiplicity = int(molecule.multiplicity)
    if updated.charge is None and molecule.charge is not None:
        updated.charge = int(molecule.charge)

    try:
        protonated = PKaMolecule.add_proton_at_atom(molecule, site_index)
    except ValueError as exc:
        if "hydrogen" in str(exc).lower():
            raise ValueError(
                "-pi/--proton-index with --pkb must be the heavy atom to "
                "protonate, not a hydrogen. Colour the basic atom, or drop "
                "--pkb and run pKa."
            ) from exc
        raise

    if updated.charge is not None:
        updated.charge = int(updated.charge) + 1
    elif protonated.charge is not None:
        updated.charge = int(protonated.charge)
    return protonated, protonated.proton_index, updated


def prepare_pka_submit_structure(
    molecule, proton_index, opt_settings, pkb=False
):
    """Return ``(molecule, proton_index, opt_settings)`` for a pKa job.

    When *pkb* is false, charge/multiplicity are taken from *molecule* when
    unset. When *pkb* is true, the free base is protonated and settings
    charge becomes the input charge plus one.
    """
    if pkb:
        return prepare_pkb_submit_molecule(
            molecule, proton_index, opt_settings
        )
    opt_settings = apply_pka_molecule_charge_multiplicity(
        opt_settings, molecule
    )
    return molecule, proton_index, opt_settings


def prepare_pka_submit_molecules(
    molecules, proton_index, opt_settings, pkb=False
):
    """Apply :func:`prepare_pka_submit_structure` to each input molecule."""
    if not molecules:
        return molecules, proton_index, opt_settings
    prepared = []
    new_index = proton_index
    updated = opt_settings
    for mol in molecules:
        mol, new_index, updated = prepare_pka_submit_structure(
            mol, proton_index, opt_settings, pkb=pkb
        )
        prepared.append(mol)
    return prepared, new_index, updated


def apply_pka_molecule_charge_multiplicity(opt_settings, molecule):
    """Use charge/multiplicity from *molecule* when absent on *opt_settings*."""
    import copy

    updated = copy.copy(opt_settings)
    mol_charge = molecule.charge
    mol_mult = molecule.multiplicity
    if updated.charge is None and mol_charge is not None:
        updated.charge = int(mol_charge)
    if updated.multiplicity is None and mol_mult is not None:
        updated.multiplicity = int(mol_mult)
    return updated


def require_pka_charge_multiplicity(opt_settings, source_hint=""):
    """Raise when charge or multiplicity are still unset after resolution."""
    missing = []
    if opt_settings.charge is None:
        missing.append("-c/--charge")
    if opt_settings.multiplicity is None:
        missing.append("-m/--multiplicity")
    if not missing:
        return
    suffix = f" ({source_hint})" if source_hint else ""
    raise click.UsageError(
        "Charge and multiplicity are required for pKa submission. "
        f"Missing: {', '.join(missing)}. "
        "Provide them on the parent command or use a ChemDraw structure "
        f"from which they can be inferred{suffix}."
    )


def is_pka_batch_invocation(ctx):
    """Return True when the nested ``pka`` command targets ``batch`` mode."""
    if getattr(ctx, "invoked_subcommand", None) != "pka":
        return False

    tokens = []
    current = ctx
    while current is not None:
        if getattr(current, "invoked_subcommand", None) == "pka":
            tokens.extend(str(token) for token in (current.args or []))
        current = current.parent

    if "submit" in tokens:
        return False
    return "batch" in tokens


PKA_SITE_SOURCE_EXPLICIT = "explicit"
PKA_SITE_SOURCE_COLOUR = "ChemDraw colour"
PKA_SITE_SOURCE_SMARTS = "SMARTS"
PKA_PREFLIGHT_COLUMNS = (
    "fragment",
    "label",
    "mode",
    "site",
    "source",
    "charge",
    "multiplicity",
    "status",
)


class PkaPreflightRow(NamedTuple):
    fragment: int
    label: str
    mode: str
    site: str
    source: str
    charge: str
    multiplicity: str
    status: str


def _smarts_pka_sites(molecules, filename, mode):
    if not molecules:
        raise ValueError(
            f"Could not read a molecule from {filename}. "
            "Specify -pi/--proton-index or -cc/--color-code."
        )
    return [
        (
            molecule,
            resolve_ionizable_site(molecule, mode=mode),
            PKA_SITE_SOURCE_SMARTS,
        )
        for molecule in molecules
    ]


def _cdxml_pka_sites(filename, color_code, mode):
    """Return colour sites, or SMARTS sites when the file has no markup."""
    from chemsmart.io.file import PKaCDXFile

    cdx_file = PKaCDXFile(filename)
    try:
        if mode == "base":
            site, molecules = cdx_file._resolve_proton_from_cdxml(
                color_code, mode=mode
            )
        else:
            molecules = cdx_file.get_pka_molecules(
                color_code=color_code,
                index=":",
                return_list=True,
            )
            return [
                (molecule, molecule.proton_index, PKA_SITE_SOURCE_COLOUR)
                for molecule in molecules
            ]
    except ValueError as exc:
        if color_code is not None or not _cdxml_has_no_colour_markup(exc):
            raise
        return _smarts_pka_sites(list(cdx_file.molecules), filename, mode)

    if molecules is not None:
        return [
            (molecule, molecule.proton_index, PKA_SITE_SOURCE_COLOUR)
            for molecule in molecules
        ]
    parsed = list(cdx_file.molecules)
    if not parsed:
        raise ValueError(f"Could not read a molecule from {filename}.")
    return [(parsed[-1], site, PKA_SITE_SOURCE_COLOUR)]


def resolve_pka_sites(
    filename, proton_index=None, color_code=None, mode="acid"
):
    """Return ``(molecule, site, source)`` for each ionizable site.

    Precedence is an explicit index, then ChemDraw colour, then SMARTS.
    A multi-fragment CDXML file returns one tuple per fragment. An
    explicit index returns one tuple and does not run colour or SMARTS.
    """
    if mode not in _IONIZABLE_SITE_SMARTS:
        raise ValueError(f"mode must be 'acid' or 'base', got {mode!r}")

    filename = str(filename)
    if proton_index is not None:
        from chemsmart.io.molecules.structure import Molecule

        molecule = Molecule.from_filepath(filename)
        return [(molecule, int(proton_index), PKA_SITE_SOURCE_EXPLICIT)]

    if is_pka_cdxml_input(filename):
        return _cdxml_pka_sites(filename, color_code, mode)

    from chemsmart.io.molecules.structure import Molecule

    molecule = Molecule.from_filepath(filename)
    molecules = [] if molecule is None else [molecule]
    return _smarts_pka_sites(molecules, filename, mode)


def resolve_proton_index(filename, proton_index, color_code=None, mode="acid"):
    """Resolve the ionizable site for pKa or ``--pkb`` submission.

    Site precedence is ``-pi``, then ChemDraw colour, then a unique SMARTS
    match. If a proton index is provided, it is returned directly. For
    CDX/CDXML inputs, a uniquely coloured site is auto-detected; a
    multi-fragment file yields a list of per-fragment molecules, which is
    returned as the second tuple element while the proton index is set to
    ``None`` so callers can branch to per-molecule job creation. When the
    CDXML file has no colour markup, SMARTS is used. For other structure
    files, a unique site is resolved from acid or base SMARTS.

    Args:
        filename: Input structure file path, used to detect CDX/CDXML inputs.
        proton_index: 1-based index supplied by the user, if any. Acid mode:
            hydrogen to remove. Base mode (``--pkb``): heavy atom to
            protonate.
        color_code: CDXML color-table index used for auto-detection.
        mode: ``"acid"`` (default) or ``"base"``.

    Returns:
        tuple[int | None, list | None]:
            - Site index when a single molecule is resolved.
            - ``None`` for the index with a list of per-fragment molecules
              when multiple molecules are detected in CDX/CDXML.

    Raises:
        ValueError: If required inputs are missing or inconsistent with the
            file type, or SMARTS finds 0 or more than one site.
    """
    if mode not in _IONIZABLE_SITE_SMARTS:
        raise ValueError(f"mode must be 'acid' or 'base', got {mode!r}")

    if proton_index is not None:
        return proton_index, None

    filename = str(filename)
    if not is_pka_cdxml_input(filename):
        if color_code is not None:
            raise ValueError(
                "-cc/--color-code can only be used with .cdx/.cdxml files."
            )
        from chemsmart.utils.datasets import PKaTableEntry

        if PKaTableEntry.is_submission_table(filename):
            raise ValueError(
                "Table input detected for pKa job submission. "
                "Use the 'batch' subcommand to process each table row, e.g. "
                "'pka -s direct batch' (proton_index is read from the table)."
            )

    sites = resolve_pka_sites(filename, color_code=color_code, mode=mode)
    if len(sites) > 1:
        molecules = []
        for molecule, site, _source in sites:
            molecule.proton_index = site
            molecules.append(molecule)
        return None, molecules
    return sites[0][1], None


def resolve_pka_batch_row(
    filepath, proton_index=None, color_code=None, pkb=False
):
    """Resolve site index and molecule for one pKa submission-table row.

    Each table row maps to a single job. An explicit ``proton_index`` always
    takes precedence (acidic H, or with ``pkb`` the heavy atom to protonate).
    Otherwise :func:`resolve_pka_sites` supplies the site. Multi-molecule
    CDXML files are rejected here; pass them directly as ``-f`` with
    ``pka batch`` instead.

    Returns:
        tuple[int, Molecule | PKaMolecule]: Resolved site index and structure.
    """
    mode = pka_submit_site_mode(pkb)
    sites = resolve_pka_sites(
        filepath,
        proton_index=proton_index,
        color_code=color_code,
        mode=mode,
    )
    if proton_index is None and len(sites) != 1:
        raise ValueError(
            f"CDXML file {filepath} contains {len(sites)} molecules. "
            "Submission-table rows support single-molecule CDXML files only. "
            "Pass a multi-molecule CDXML file directly as -f with pka batch."
        )
    molecule, site, _source = sites[0]
    return site, molecule


def list_pka_preflight_sites(filename, color_code, mode):
    """Return ``(molecule, site_index, source)`` in ChemDraw fragment order.

    An explicit ``--proton-index`` is not applied here. CDXML colour is
    tried first. A file with no colour markup falls through to SMARTS.
    Other structure files use SMARTS.
    """
    return resolve_pka_sites(filename, color_code=color_code, mode=mode)


def _pka_submission_backend(ctx):
    current = ctx
    while current is not None:
        if current.info_name in ("gaussian", "orca"):
            return current.info_name
        current = current.parent
    return "gaussian"


def _single_pka_job_label(ctx, filename):
    label = ctx.obj.get("label")
    if label:
        return label
    return f"{Path(str(filename)).stem}_pka"


def _fragment_pka_job_label(filename, fragment):
    return f"{Path(str(filename)).stem}_frag{fragment}_pka"


def _table_pka_job_label(filepath, backend):
    label = Path(str(filepath)).stem
    if backend == "orca":
        return ensure_pka_label_suffix(label)
    return label


def _preview_opt_settings(ctx):
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    job_settings = ctx.obj.get("job_settings")
    keywords = ctx.obj.get("keywords") or {}
    if job_settings is not None:
        opt_settings = opt_settings.merge(job_settings, keywords=keywords)
    return opt_settings


def _resolved_charge_multiplicity(opt_settings, molecule):
    charge = opt_settings.charge
    multiplicity = opt_settings.multiplicity
    if charge is None and molecule.charge is not None:
        charge = int(molecule.charge)
    if multiplicity is None and molecule.multiplicity is not None:
        multiplicity = int(molecule.multiplicity)
    return charge, multiplicity


def _fragment_charge_prefix(fragment, filename):
    if is_pka_cdxml_input(filename):
        return f"ChemDraw molecular fragment {fragment}"
    return f"Fragment {fragment}"


def _require_preview_charge_multiplicity(
    fragment, charge, multiplicity, filename
):
    missing = []
    if charge is None:
        missing.append("-c/--charge")
    if multiplicity is None:
        missing.append("-m/--multiplicity")
    if not missing:
        return
    prefix = _fragment_charge_prefix(fragment, filename)
    raise click.UsageError(
        f"{prefix}: charge and multiplicity are required before pKa "
        f"preview. Missing: {', '.join(missing)}. "
        "Provide them on the parent command or use a structure from "
        f"which they can be inferred ({filename})."
    )


def _preflight_row_for_site(
    ctx,
    molecule,
    site_index,
    source,
    fragment,
    label,
    filename,
    opt_settings=None,
    charge=None,
    multiplicity=None,
):
    shared = ctx.obj["pka_shared"]
    pkb = bool(shared.get("pkb", False))
    mode_name = "pKb" if pkb else "pKa"
    site_index = int(site_index)
    n_atoms = molecule.num_atoms
    if site_index < 1 or site_index > n_atoms:
        raise click.UsageError(
            f"Fragment {fragment}: atom index {site_index} is outside "
            f"1..{n_atoms} ({filename})."
        )
    element = molecule.chemical_symbols[site_index - 1]
    if opt_settings is None:
        opt_settings = _preview_opt_settings(ctx)
    if charge is None and multiplicity is None:
        charge, multiplicity = _resolved_charge_multiplicity(
            opt_settings, molecule
        )
    _require_preview_charge_multiplicity(
        fragment, charge, multiplicity, filename
    )

    status = "ok"
    if pkb:
        try:
            protonated, _hydrogen_index, _updated = (
                prepare_pkb_submit_molecule(molecule, site_index, opt_settings)
            )
        except ValueError as exc:
            raise click.UsageError(
                f"{_fragment_charge_prefix(fragment, filename)}: {exc}"
            ) from exc
        status = f"ok; added H {protonated.proton_index}"
    elif element != "H":
        raise click.UsageError(
            f"{_fragment_charge_prefix(fragment, filename)}: atom "
            f"{site_index} is {element}, not the hydrogen that would "
            f"be removed ({filename})."
        )

    return PkaPreflightRow(
        fragment=fragment,
        label=label,
        mode=mode_name,
        site=f"{site_index} {element}",
        source=source,
        charge=str(int(charge)),
        multiplicity=str(int(multiplicity)),
        status=status,
    )


def _molecules_for_explicit_preview(filename):
    if is_pka_cdxml_input(filename):
        from chemsmart.io.file import PKaCDXFile

        return list(PKaCDXFile(filename).molecules)
    from chemsmart.io.molecules.structure import Molecule

    molecule = Molecule.from_filepath(filename)
    if molecule is None:
        return []
    if isinstance(molecule, list):
        return list(molecule)
    return [molecule]


def _explicit_preflight_targets(ctx, filename):
    molecules = list(ctx.obj.get("molecules") or [])
    if not molecules:
        molecules = _molecules_for_explicit_preview(filename)
    if not molecules:
        raise click.UsageError(f"Could not read a molecule from {filename}.")
    indices = ctx.obj.get("molecule_indices")
    if indices and len(molecules) > 1:
        return [
            (int(index), molecule)
            for index, molecule in zip(indices, molecules)
        ]
    fragment = len(molecules)
    return [(fragment, molecules[-1])]


def _rows_for_explicit_index(ctx, filename, proton_index, opt_settings):
    targets = _explicit_preflight_targets(ctx, filename)
    base_label = _single_pka_job_label(ctx, filename)
    multiple = len(targets) > 1
    rows = []
    for fragment, molecule in targets:
        label = f"{base_label}_idx{fragment}" if multiple else base_label
        rows.append(
            _preflight_row_for_site(
                ctx,
                molecule,
                proton_index,
                PKA_SITE_SOURCE_EXPLICIT,
                fragment,
                label,
                filename,
                opt_settings=opt_settings,
            )
        )
    return rows


def _rows_for_submission_table(ctx, filename, color_code, mode, opt_settings):
    import copy

    from chemsmart.utils.datasets import PKaOutputTable, PKaTableEntry

    try:
        entries = PKaTableEntry.parse_pka_table(filename)
        PKaOutputTable.validate_pka_table_entries(
            entries, check_file_exists=True
        )
    except (FileNotFoundError, ValueError) as exc:
        raise click.UsageError(str(exc)) from exc

    backend = _pka_submission_backend(ctx)
    rows = []
    for number, entry in enumerate(entries, start=1):
        filepath = entry.filepath
        label = _table_pka_job_label(filepath, backend)
        row_settings = copy.copy(opt_settings)
        row_settings.charge = int(entry.charge)
        row_settings.multiplicity = int(entry.multiplicity)
        sites = resolve_pka_sites(
            filepath,
            proton_index=entry.proton_index,
            color_code=color_code,
            mode=mode,
        )
        if entry.proton_index is None and len(sites) != 1:
            raise click.UsageError(
                f"Row {number}: {filepath} contains {len(sites)} "
                "ChemDraw fragments. Pass a multi-fragment file as "
                "-f with pka batch, not inside a submission table."
            )
        molecule, site, source = sites[0]
        rows.append(
            _preflight_row_for_site(
                ctx,
                molecule,
                site,
                source,
                number,
                label,
                filepath,
                opt_settings=row_settings,
                charge=int(entry.charge),
                multiplicity=int(entry.multiplicity),
            )
        )
    return rows


def build_pka_preflight_rows(ctx):
    """Return preflight rows for the current pKa submission context."""
    filename = ctx.obj.get("filename")
    if not filename:
        raise click.UsageError(
            "--preview requires -f/--filename pointing at a structure "
            "file or a pKa submission table."
        )
    proton_index, color_code = resolve_pka_submit_proton_options(ctx)
    mode = pka_submit_site_mode(ctx.obj["pka_shared"].get("pkb", False))
    opt_settings = _preview_opt_settings(ctx)

    from chemsmart.utils.datasets import PKaTableEntry

    if PKaTableEntry.is_submission_table(filename):
        return _rows_for_submission_table(
            ctx, filename, color_code, mode, opt_settings
        )
    if proton_index is not None:
        return _rows_for_explicit_index(
            ctx, filename, proton_index, opt_settings
        )

    sites = list_pka_preflight_sites(filename, color_code, mode)
    multiple = len(sites) > 1
    rows = []
    for number, (molecule, site, source) in enumerate(sites, start=1):
        label = (
            _fragment_pka_job_label(filename, number)
            if multiple
            else _single_pka_job_label(ctx, filename)
        )
        rows.append(
            _preflight_row_for_site(
                ctx,
                molecule,
                site,
                source,
                number,
                label,
                filename,
                opt_settings=opt_settings,
            )
        )
    return rows


def format_pka_preflight_table(rows):
    """Return a deterministic plain-text table for *rows*."""
    records = [list(row) for row in rows]
    headers = list(PKA_PREFLIGHT_COLUMNS)
    widths = [
        max([len(header)] + [len(str(record[index])) for record in records])
        for index, header in enumerate(headers)
    ]

    def _format_line(cells):
        return "  ".join(
            str(cell).ljust(width) for cell, width in zip(cells, widths)
        ).rstrip()

    lines = [
        _format_line(headers),
        _format_line("-" * width for width in widths),
    ]
    lines.extend(_format_line(record) for record in records)
    return "\n".join(lines)


def _parent_charge_multiplicity(ctx):
    job_settings = ctx.obj.get("job_settings")
    if job_settings is None:
        return None, None
    return job_settings.charge, job_settings.multiplicity


def pka_preflight_notices(ctx):
    """Return notices for a multi-fragment ChemDraw file.

    Parent ``-c`` / ``-m`` are reported only when those command values
    are set, because they then apply to every fragment. Per-fragment
    charges parsed from the drawing are left to the table.
    """
    filename = ctx.obj.get("filename")
    if not filename or not is_pka_cdxml_input(filename):
        return []
    if len(PKaCDXFile(filename).molecules) < 2:
        return []
    notices = [CHEMDRAW_MOLECULAR_FRAGMENT_WARNING]
    charge, multiplicity = _parent_charge_multiplicity(ctx)
    parts = []
    if charge is not None:
        parts.append(f"-c/--charge {int(charge)}")
    if multiplicity is not None:
        parts.append(f"-m/--multiplicity {int(multiplicity)}")
    if not parts:
        return notices
    verb = "apply" if len(parts) > 1 else "applies"
    notices.append(
        f"Parent {' and '.join(parts)} {verb} to every "
        "ChemDraw molecular fragment."
    )
    return notices


def print_pka_preview(ctx):
    """Print the preflight table when ``--preview`` is set.

    Returns True after printing. The caller must not create jobs.
    """
    if not ctx.obj["pka_shared"].get("preview"):
        return False
    try:
        rows = build_pka_preflight_rows(ctx)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc
    for notice in pka_preflight_notices(ctx):
        click.echo(notice)
    click.echo(format_pka_preflight_table(rows))
    return True


def ensure_pka_label_suffix(raw):
    """Append ``_pka`` when *raw* does not already end with that suffix."""
    text = "" if raw is None else str(raw)
    if text.endswith("_pka"):
        return text
    return f"{text}_pka"


def log_pka_job_settings(pka_settings, proton_index, shared):
    """Log the proton site, thermodynamic cycle, and species charges."""
    logger.info(f"Proton index to remove: {proton_index}")
    logger.info(f"Thermodynamic cycle: {shared['scheme']}")
    logger.info(
        f"Protonated form (HA): charge={pka_settings.charge}, "
        f"mult={pka_settings.multiplicity}"
    )
    if shared["conjugate_base_charge"] is not None:
        cb_charge = shared["conjugate_base_charge"]
    else:
        cb_charge = pka_settings.charge - 1
    if shared["conjugate_base_multiplicity"] is not None:
        cb_mult = shared["conjugate_base_multiplicity"]
    else:
        cb_mult = pka_settings.multiplicity
    logger.info(f"Conjugate base (A-): charge={cb_charge}, mult={cb_mult}")


_PKA_REFERENCE_SHARED_KEYS = (
    "reference",
    "reference_proton_index",
    "reference_charge",
    "reference_multiplicity",
    "reference_conjugate_base_charge",
    "reference_conjugate_base_multiplicity",
)


def configure_pka_submission(ctx, *, submit_command, batch_command, **options):
    """Store shared pKa options and dispatch to submit or batch.

    *options* are the Click values from ``gaussian pka`` or ``orca pka``.
    When the group is invoked without a subcommand, a submission table
    selects *batch_command* and every other input selects *submit_command*.
    """
    s_freq_cutoff, entropy_method = resolve_pka_entropy_cutoff(
        options["cutoff_entropy_grimme"],
        options["cutoff_entropy_truhlar"],
    )
    sampling, num_conformers = resolve_pka_sampling_options(
        options["sampling"], options["num_conformers"]
    )
    shared = dict(
        scheme=options["scheme"],
        reference=options["reference"],
        reference_proton_index=options["reference_proton_index"],
        reference_color_code=options["reference_color_code"],
        reference_charge=options["reference_charge"],
        reference_multiplicity=options["reference_multiplicity"],
        reference_conjugate_base_charge=options[
            "reference_conjugate_base_charge"
        ],
        reference_conjugate_base_multiplicity=options[
            "reference_conjugate_base_multiplicity"
        ],
        delta_g_proton=options["delta_g_proton"],
        conjugate_base_charge=options["conjugate_base_charge"],
        conjugate_base_multiplicity=options["conjugate_base_multiplicity"],
        solvent_model=options["solvent_model"],
        solvent_id=options["solvent_id"],
        sampling=sampling,
        num_conformers=num_conformers,
        pkb=options["pkb"],
        pks=options["pks"],
        temperature=options["temperature"],
        concentration=options["concentration"],
        pressure=options["pressure"],
        cutoff_entropy_grimme=s_freq_cutoff,
        cutoff_enthalpy=options["cutoff_enthalpy"],
        entropy_method=entropy_method,
        skip_completed=options["skip_completed"],
        preview=options["preview"],
    )
    ctx.ensure_object(dict)
    ctx.obj["pka_shared"] = shared
    ctx.obj["pka_proton_index"] = options["proton_index"]
    ctx.obj["pka_color_code"] = options["color_code"]

    if ctx.invoked_subcommand is not None:
        return None

    from chemsmart.utils.datasets import PKaTableEntry

    filename = ctx.obj.get("filename")
    skip_completed = options["skip_completed"]
    if PKaTableEntry.is_submission_table(filename):
        return ctx.invoke(batch_command, skip_completed=skip_completed)
    return ctx.invoke(submit_command, skip_completed=skip_completed)


def _pka_opt_settings(ctx, *, mode):
    """Return project settings and merged optimization settings.

    ``submit`` always merges job settings. ``fragment`` merges when job
    settings are present and requires the context keys. ``batch`` merges
    when job settings are present and tolerates missing context keys.
    """
    project_settings = ctx.obj["project_settings"]
    opt_settings = project_settings.opt_settings()
    if mode == "batch":
        job_settings = ctx.obj.get("job_settings")
        keywords = ctx.obj.get("keywords", {})
    else:
        job_settings = ctx.obj["job_settings"]
        keywords = ctx.obj["keywords"]
    if mode == "submit" or job_settings:
        opt_settings = opt_settings.merge(job_settings, keywords=keywords)
    return project_settings, opt_settings


def _require_batch_reference_options(shared):
    """Require reference acid options for a proton-exchange batch."""
    if shared["scheme"] != "proton exchange":
        return

    missing = []
    if shared["reference"] is None:
        missing.append("-r/--reference")
    elif shared["reference_proton_index"] is None:
        ref = shared["reference"]
        if str(ref).endswith((".cdx", ".cdxml")):
            try:
                shared["reference_proton_index"] = (
                    PKaCDXFile.resolve_reference_proton(
                        ref,
                        None,
                        shared["reference_color_code"],
                    )
                )
            except click.UsageError:
                missing.append("-rpi/--reference-proton-index")
        else:
            missing.append("-rpi/--reference-proton-index")
    if shared["reference_charge"] is None:
        missing.append("-rc/--reference-charge")
    if shared["reference_multiplicity"] is None:
        missing.append("-rm/--reference-multiplicity")
    if missing:
        raise click.UsageError(
            "For proton exchange cycle with batch input, these "
            "reference acid options are required:\n  " + "\n  ".join(missing)
        )


def _batch_row_shared(shared, index):
    """Keep the reference acid on the first proton-exchange row only."""
    row_shared = copy.copy(shared)
    original_scheme = shared["scheme"]
    if index == 0 or original_scheme != "proton exchange":
        row_shared["scheme"] = original_scheme
    else:
        row_shared["scheme"] = "direct"
    if row_shared["scheme"] != "proton exchange":
        for key in _PKA_REFERENCE_SHARED_KEYS:
            row_shared[key] = None
    return row_shared


def _make_pka_job(
    job_class,
    molecule,
    settings,
    label,
    jobrunner,
    skip_completed,
    **kwargs,
):
    return job_class(
        molecule=molecule,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )


def _create_pka_jobs_from_molecules(
    ctx,
    pka_molecules,
    shared,
    skip_completed,
    *,
    job_class,
    settings_builder,
    label_for,
    program,
    **kwargs,
):
    """Create one pKa job per ChemDraw fragment."""
    validate_reference_options(shared)
    project_settings, opt_settings = _pka_opt_settings(ctx, mode="fragment")
    jobrunner = ctx.obj["jobrunner"]
    filename = ctx.obj.get("filename", "")
    basename = Path(str(filename)).stem or "pka"

    jobs = []
    for idx, pka_mol in enumerate(pka_molecules, start=1):
        label = label_for(f"{basename}_frag{idx}_pka")
        try:
            molecule, proton_index, row_opt_settings = (
                prepare_pka_submit_structure(
                    pka_mol,
                    pka_mol.proton_index,
                    opt_settings,
                    pkb=shared.get("pkb", False),
                )
            )
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc
        require_pka_charge_multiplicity(
            row_opt_settings,
            source_hint=(f"ChemDraw molecular fragment {idx} in {filename}"),
        )
        pka_settings = settings_builder(
            proton_index, shared, row_opt_settings, project_settings
        )
        logger.info(
            f"Creating {program} pKa job for fragment {idx}: "
            f"proton_index={proton_index}, label={label}"
        )
        job = _make_pka_job(
            job_class,
            molecule,
            pka_settings,
            label,
            jobrunner,
            skip_completed,
            **kwargs,
        )
        charge = pka_mol.charge
        if charge is None:
            charge = pka_settings.charge
        job._batch_entry = {
            "filepath": str(filename),
            "proton_index": pka_mol.proton_index,
            "charge": int(charge),
            "multiplicity": int(pka_settings.multiplicity),
            "scheme": shared["scheme"],
            "fragment_index": idx,
            "label": label,
        }
        jobs.append(job)

    logger.info(
        f"Created {len(jobs)} {program} pKa jobs from multi-fragment CDXML"
    )
    return jobs


def submit_pka_jobs(
    ctx,
    skip_completed,
    proton_index,
    color_code,
    *,
    job_class,
    settings_builder,
    label_for,
    batch_command,
    program,
    **kwargs,
):
    """Create pKa jobs for one structure, or expand a multi-fragment CDXML.

    *settings_builder* is called as
    ``settings_builder(proton_index, shared, opt_settings, project_settings)``.
    *label_for* maps a raw label to the program's job label.
    A submission table is handed to *batch_command*.
    """
    shared = ctx.obj["pka_shared"]
    if print_pka_preview(ctx):
        return None
    filename = ctx.obj.get("filename")

    from chemsmart.utils.datasets import PKaTableEntry

    if PKaTableEntry.is_submission_table(filename):
        return ctx.invoke(batch_command, skip_completed=skip_completed)

    proton_index, color_code = resolve_pka_submit_proton_options(
        ctx, proton_index=proton_index, color_code=color_code
    )
    try:
        proton_index, pka_molecules = resolve_proton_index(
            filename,
            proton_index,
            color_code,
            mode=pka_submit_site_mode(shared.get("pkb", False)),
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc

    if pka_molecules is not None:
        return _create_pka_jobs_from_molecules(
            ctx,
            pka_molecules,
            shared,
            skip_completed,
            job_class=job_class,
            settings_builder=settings_builder,
            label_for=label_for,
            program=program,
            **kwargs,
        )

    validate_reference_options(shared)
    project_settings, opt_settings = _pka_opt_settings(ctx, mode="submit")
    molecules = ctx.obj["molecules"]
    try:
        molecules, proton_index, opt_settings = prepare_pka_submit_molecules(
            molecules,
            proton_index,
            opt_settings,
            pkb=shared.get("pkb", False),
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc

    pka_settings = settings_builder(
        proton_index, shared, opt_settings, project_settings
    )
    require_pka_charge_multiplicity(
        pka_settings, source_hint=f"input file {filename}"
    )
    log_pka_job_settings(pka_settings, proton_index, shared)

    jobrunner = ctx.obj["jobrunner"]
    label = label_for(ctx.obj["label"])
    molecule_indices = ctx.obj.get("molecule_indices")
    if len(molecules) > 1 and molecule_indices:
        logger.info(f"Creating {len(molecules)} {program} pKa jobs")
        return [
            _make_pka_job(
                job_class,
                mol,
                pka_settings,
                f"{label}_idx{idx}",
                jobrunner,
                skip_completed,
                **kwargs,
            )
            for mol, idx in zip(molecules, molecule_indices)
        ]

    return _make_pka_job(
        job_class,
        molecules[-1],
        pka_settings,
        label,
        jobrunner,
        skip_completed,
        **kwargs,
    )


def batch_pka_jobs(
    ctx,
    skip_completed,
    proton_index,
    color_code,
    *,
    job_class,
    settings_builder,
    label_for,
    submit_command,
    program,
    **kwargs,
):
    """Create pKa jobs from a submission table or a multi-molecule CDXML.

    Only the first proton-exchange row keeps the reference acid. Later
    rows use the direct cycle. *settings_builder* and *label_for* match
    :func:`submit_pka_jobs`.
    """
    shared = ctx.obj["pka_shared"]
    if print_pka_preview(ctx):
        return None

    input_table_path = ctx.obj.get("filename")
    if not input_table_path:
        raise click.UsageError(
            "Batch mode requires the parent "
            f"{program} -f/--filename to specify the table file path."
        )

    def create_fragment_jobs(
        ctx, pka_molecules, shared, skip_completed, **fragment_kwargs
    ):
        return _create_pka_jobs_from_molecules(
            ctx,
            pka_molecules,
            shared,
            skip_completed,
            job_class=job_class,
            settings_builder=settings_builder,
            label_for=label_for,
            program=program,
            **fragment_kwargs,
        )

    def invoke_submit(ctx, **invoke_kwargs):
        return ctx.invoke(submit_command, **invoke_kwargs)

    if is_pka_cdxml_input(input_table_path):
        return batch_pka_jobs_from_cdxml(
            ctx,
            skip_completed,
            create_fragment_jobs,
            invoke_submit,
            **kwargs,
        )

    from chemsmart.utils.datasets import PKaOutputTable, PKaTableEntry

    logger.info(f"Reading {program} pKa jobs from table: {input_table_path}")
    try:
        entries = PKaTableEntry.parse_pka_table(input_table_path)
        PKaOutputTable.validate_pka_table_entries(
            entries, check_file_exists=True
        )
    except (FileNotFoundError, ValueError) as exc:
        raise click.UsageError(str(exc)) from exc

    logger.info(f"Found {len(entries)} entries in table")
    _require_batch_reference_options(shared)

    project_settings, opt_settings = _pka_opt_settings(ctx, mode="batch")
    _, color_code = resolve_pka_submit_proton_options(
        ctx, proton_index=proton_index, color_code=color_code
    )
    jobrunner = ctx.obj["jobrunner"]
    jobs = []
    for index, entry in enumerate(entries):
        filepath = entry.get("filepath") or entry.get("path") or entry.filepath
        try:
            row_proton_index, molecule = resolve_pka_batch_row(
                filepath,
                proton_index=entry.proton_index,
                color_code=color_code,
                pkb=shared.get("pkb", False),
            )
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc
        label = label_for(Path(filepath).stem)
        input_proton_index = row_proton_index
        input_charge = int(entry.charge)

        row_opt_settings = copy.copy(opt_settings)
        row_opt_settings.charge = input_charge
        row_opt_settings.multiplicity = int(entry.multiplicity)
        try:
            molecule, row_proton_index, row_opt_settings = (
                prepare_pka_submit_structure(
                    molecule,
                    row_proton_index,
                    row_opt_settings,
                    pkb=shared.get("pkb", False),
                )
            )
        except ValueError as exc:
            raise click.UsageError(str(exc)) from exc

        row_shared = _batch_row_shared(shared, index)
        pka_settings = settings_builder(
            row_proton_index,
            row_shared,
            row_opt_settings,
            project_settings,
        )
        job = _make_pka_job(
            job_class,
            molecule,
            pka_settings,
            label,
            jobrunner,
            skip_completed,
            **kwargs,
        )
        job._batch_entry = {
            "filepath": str(filepath),
            "proton_index": input_proton_index,
            "charge": input_charge,
            "multiplicity": int(entry.multiplicity),
            "scheme": row_shared["scheme"],
            "label": label,
        }
        jobs.append(job)

    logger.info(f"Created {len(jobs)} {program} pKa jobs from table")
    return jobs


def batch_pka_jobs_from_cdxml(
    ctx,
    skip_completed,
    create_jobs_fn,
    invoke_submit_fn,
    **kwargs,
):
    """Create pKa jobs from a CDXML batch input via coloured-proton detection."""
    filename = ctx.obj.get("filename")
    shared = ctx.obj["pka_shared"]
    proton_index, color_code = resolve_pka_submit_proton_options(ctx)
    try:
        proton_index, pka_molecules = resolve_proton_index(
            filename,
            proton_index,
            color_code,
            mode=pka_submit_site_mode(shared.get("pkb", False)),
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc

    if pka_molecules is not None:
        if len(pka_molecules) > 1:
            click.echo(CHEMDRAW_MOLECULAR_FRAGMENT_WARNING)
        return create_jobs_fn(
            ctx, pka_molecules, shared, skip_completed, **kwargs
        )

    return invoke_submit_fn(
        ctx,
        skip_completed=skip_completed,
        proton_index=proton_index,
        color_code=color_code,
        **kwargs,
    )


def resolve_pka_submit_proton_options(ctx, proton_index=None, color_code=None):
    """Resolve proton options for ``pka submit`` from multiple Click scopes.

    The ``pka`` group captures ``-pi/--proton-index``, but nested invocation
    via ``chemsmart run`` does not always populate ``ctx.obj`` from group-level
    options.  ``submit`` therefore re-declares the same options and this helper
    merges submit args, parent group params, and ``ctx.obj``.
    """
    parent = ctx.parent
    if proton_index is None and parent is not None:
        proton_index = parent.params.get("proton_index")
    if proton_index is None:
        proton_index = ctx.obj.get("pka_proton_index")

    if color_code is None and parent is not None:
        color_code = parent.params.get("color_code")
    if color_code is None:
        color_code = ctx.obj.get("pka_color_code")

    ctx.obj["pka_proton_index"] = proton_index
    ctx.obj["pka_color_code"] = color_code
    return proton_index, color_code


def click_pka_proton_options(f):
    """Options that identify the proton to remove (-pi / -cc).

    Applied to the ``pka`` group and ``submit`` so values are available whether
    options appear before or after the ``submit`` token (required for per-row
    ``chemsmart run`` scripts generated by ``chemsmart sub`` batch mode).
    """

    @click.option(
        "-pi",
        "--proton-index",
        type=int,
        required=False,
        help=(
            "1-based index of the proton to remove, or with --pkb the heavy "
            "atom to protonate. If omitted, a uniquely coloured ChemDraw "
            "site or a unique SMARTS match is used."
        ),
    )
    @click.option(
        "-cc",
        "--color-code",
        type=int,
        default=None,
        help=(
            "CDXML colour-table index identifying the proton to remove. "
            "If omitted, the uniquely coloured hydrogen is auto-detected."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def click_pka_analysis_scheme_options(f):
    """Scheme and proton solvation options for pKa output analysis."""
    f = click_pka_pkb_options(f)

    @click.option(
        "-s",
        "--scheme",
        "scheme",
        type=click.Choice(["direct", "proton exchange"]),
        default=None,
        help=(
            "Thermodynamic cycle for analysis. "
            "Default: proton exchange when omitted."
        ),
    )
    @click.option(
        "-dG",
        "--delta-g-proton",
        "delta_g_proton",
        type=float,
        default=None,
        help=(
            "G_soln(H+) in kcal/mol for the direct cycle. If omitted, a "
            "T-dependent aqueous default is computed. Used only with "
            "--scheme direct."
        ),
    )
    @functools.wraps(f)
    def wrapper(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper


def _resolve_pka_analysis_scheme(scheme, delta_g_proton):
    """Validate CLI options and return the active analysis scheme."""
    if scheme is None:
        scheme = "proton exchange"

    if delta_g_proton is not None and scheme != "direct":
        logger.info(
            "Ignoring -dG/--delta-g-proton because --scheme direct was not "
            "specified; using proton exchange analysis."
        )
    return scheme


def _scheme_display_name(scheme):
    names = {
        "direct": "Direct Dissociation",
        "proton exchange": "Proton Exchange",
    }
    return names.get(scheme, scheme)


def _maybe_discover_ha_ensemble(ha_gas_path, program=None):
    """Return the HA ``_c*`` ensemble list, or ``None`` if none exist."""
    from chemsmart.utils.datasets import (
        discover_pka_output_path,
        pka_output_basename_from_path,
    )

    directory = os.path.dirname(str(ha_gas_path)) or "."
    basename = pka_output_basename_from_path(ha_gas_path, "ha_gas")
    ha_discovered = discover_pka_output_path(
        basename,
        directory,
        "ha_gas",
        program=program,
        filepath_hint=ha_gas_path,
    )
    if isinstance(ha_discovered, (list, tuple)):
        return list(ha_discovered)
    return None


def validate_direct_analyze_files(ha, a, ha_solv, a_solv):
    """Validate required files for direct-cycle pKa analysis."""
    required = [
        ("ha", ha),
        ("a", a),
        ("ha-solv", ha_solv),
        ("a-solv", a_solv),
    ]
    missing = [name for name, path in required if path is None]
    if missing:
        raise click.UsageError(
            "For direct-cycle pKa analysis all four output files are required.\n"
            f"Missing: {', '.join(f'--{name}' for name in missing)}"
        )
    _require_pka_analysis_files(required)


def _auto_discover_direct_pka_files(ha_gas_path, program=None):
    """Infer A- and solvent SP paths from the HA gas-phase output path."""
    from chemsmart.utils.datasets import (
        PKA_TARGET_SUFFIX_HELP,
        pka_field_paths_exist,
    )
    from chemsmart.utils.io import discover_pka_target_companion_outputs

    results = discover_pka_target_companion_outputs(
        ha_gas_path, program=program
    )
    ha_ensemble = _maybe_discover_ha_ensemble(ha_gas_path, program=program)
    if ha_ensemble is not None:
        results["ha"] = ha_ensemble

    missing = [
        f"  {k}: {v}"
        for k, v in results.items()
        if not pka_field_paths_exist(v)
    ]
    if missing:
        raise click.UsageError(
            "Auto-discovery could not find some companion output files.\n"
            "Missing files:\n" + "\n".join(missing) + "\n\n"
            "Provide them explicitly or ensure output files follow:\n"
            + PKA_TARGET_SUFFIX_HELP
        )
    return results


def click_pka_analyze_options(f):
    f = click.option(
        "-rp",
        "--reference-pka",
        type=float,
        default=None,
        help=(
            "Experimental pKa of reference acid HRef. Required for proton exchange analysis."
        ),
    )(f)
    f = click.option(
        "-rs",
        "--ref-solv",
        "ref_solv",
        type=click.Path(exists=True),
        default=None,
        help="Ref- solvent SP output file.",
    )(f)
    f = click.option(
        "-hrs",
        "--href-solv",
        "href_solv",
        type=click.Path(exists=True),
        default=None,
        help="HRef solvent SP output file.",
    )(f)
    f = click.option(
        "-as",
        "--a-solv",
        "a_solv",
        type=click.Path(exists=True),
        default=None,
        help="A- solvent SP output file.",
    )(f)
    f = click.option(
        "-has",
        "--ha-solv",
        "ha_solv",
        type=click.Path(exists=True),
        default=None,
        help="HA solvent SP output file.",
    )(f)
    f = click.option(
        "-r",
        "--ref",
        "ref",
        type=click.Path(exists=True),
        default=None,
        help="Ref- gas-phase opt+freq output file.",
    )(f)
    f = click.option(
        "-hr",
        "--href",
        "href",
        type=click.Path(exists=True),
        default=None,
        help="HRef gas-phase opt+freq output file.",
    )(f)
    f = click.option(
        "-a",
        "--a",
        "a",
        type=click.Path(exists=True),
        default=None,
        help="A- gas-phase opt+freq output file.",
    )(f)
    f = click.option(
        "-ha",
        "--ha",
        "ha",
        type=click.Path(exists=True),
        default=None,
        help="HA gas-phase opt+freq output file.",
    )(f)
    return f


def validate_reference_options(shared):
    reference = shared["reference"]
    if reference is None:
        if shared["scheme"] == "proton exchange":
            raise click.UsageError(
                "Proton exchange cycle requires a reference acid. "
                "Use -r/--reference."
            )
        return

    if shared["scheme"] != "proton exchange":
        raise click.UsageError(
            "Reference acid file can only be used with 'proton exchange' "
            "cycle. Use -s 'proton exchange' or remove the -r option."
        )

    shared["reference_proton_index"] = PKaCDXFile.resolve_reference_proton(
        reference,
        shared["reference_proton_index"],
        shared["reference_color_code"],
    )

    missing = []
    if shared["reference_proton_index"] is None:
        missing.append(
            "-rpi/--reference-proton-index (or use a .cdxml reference with a coloured proton)"
        )
    if shared["reference_charge"] is None:
        missing.append("-rc/--reference-charge")
    if shared["reference_multiplicity"] is None:
        missing.append("-rm/--reference-multiplicity")
    if missing:
        raise click.UsageError(
            "When --reference is provided, the following options are "
            "required: " + ", ".join(missing)
        )


def _validate_pka_table_program(pka_table, program):
    """Ensure explicit -p matches every output file in the table."""
    from chemsmart.utils.datasets import pka_field_paths

    output_fields = (
        "ha_gas",
        "a_gas",
        "ha_sp",
        "a_sp",
        "href_gas",
        "ref_gas",
        "href_sp",
        "ref_sp",
    )
    for entry in pka_table.entries:
        for field in output_fields:
            for path in pka_field_paths(entry.get(field)):
                if not os.path.isfile(path):
                    continue
                detected = get_program_type_from_file(path)
                if detected != program:
                    raise click.UsageError(
                        f"File '{path}' was detected as {detected!r}, but "
                        f"batch-analyze was run with -p {program}."
                    )


def _require_pka_analysis_files(named_paths):
    """Raise UsageError if any named analysis path is missing on disk."""
    from chemsmart.utils.datasets import pka_field_paths

    missing_files = []
    for name, path in named_paths:
        paths = pka_field_paths(path)
        if not paths:
            missing_files.append(f"  --{name}: {path}")
            continue
        for filepath in paths:
            if not os.path.isfile(filepath):
                missing_files.append(f"  --{name}: {filepath}")
    if missing_files:
        raise click.UsageError(
            "One or more pKa analysis files do not exist:\n"
            + "\n".join(missing_files)
        )


def validate_analyze_files(
    ha, a, href, ref, ha_solv, a_solv, href_solv, ref_solv, reference_pka
):
    required_gas = [ha, a, href, ref]
    required_solv = [ha_solv, a_solv, href_solv, ref_solv]
    file_names = ["ha", "a", "href", "ref"]

    missing_gas = []
    missing_solv = []
    for gas, solv, name in zip(required_gas, required_solv, file_names):
        if gas is None:
            missing_gas.append(f"--{name}")
        if solv is None:
            missing_solv.append(f"--{name}-solv")

    if missing_gas or missing_solv:
        raise click.UsageError(
            "For pKa analysis all 8 output files are required.\n"
            f"Missing gas-phase: "
            f"{', '.join(missing_gas) if missing_gas else 'none'}\n"
            f"Missing solvent SP: "
            f"{', '.join(missing_solv) if missing_solv else 'none'}"
        )

    if reference_pka is None:
        raise click.UsageError(
            "-rp/--reference-pka is required for output-file analysis."
        )

    solv_names = [f"{name}-solv" for name in file_names]
    named_paths = list(zip(file_names, required_gas)) + list(
        zip(solv_names, required_solv)
    )
    _require_pka_analysis_files(named_paths)


def _auto_discover_pka_files(ha_gas_path, href_gas_path, program=None):
    """Infer companion output paths from HA and HRef gas-phase paths."""
    from chemsmart.utils.datasets import (
        PKA_REFERENCE_SUFFIX_HELP,
        PKA_TARGET_SUFFIX_HELP,
        discover_pka_reference_companion_outputs,
        pka_field_paths_exist,
    )
    from chemsmart.utils.io import discover_pka_target_companion_outputs

    if program is None:
        program = get_program_type_from_file(ha_gas_path)

    results = discover_pka_target_companion_outputs(
        ha_gas_path, program=program
    )
    results.update(discover_pka_reference_companion_outputs(href_gas_path))
    ha_ensemble = _maybe_discover_ha_ensemble(ha_gas_path, program=program)
    if ha_ensemble is not None:
        results["ha"] = ha_ensemble

    missing = [
        f"  {k}: {v}"
        for k, v in results.items()
        if not pka_field_paths_exist(v)
    ]
    if missing:
        raise click.UsageError(
            "Auto-discovery could not find some companion output files.\n"
            "Missing files:\n" + "\n".join(missing) + "\n\n"
            "Provide them explicitly or ensure output files follow:\n"
            + PKA_TARGET_SUFFIX_HELP
            + "\n"
            + PKA_REFERENCE_SUFFIX_HELP
        )
    return results


@click.group(name="pka", cls=MyGroup)
@click_pka_thermochemistry_options
@click_pka_analysis_scheme_options
@click.pass_context
def pka(
    ctx,
    temperature,
    concentration,
    pressure,
    cutoff_entropy_grimme,
    cutoff_entropy_truhlar,
    cutoff_enthalpy,
    scheme,
    delta_g_proton,
    pkb,
    pks,
):
    """Backend-independent pKa output analysis."""
    s_freq_cutoff, entropy_method = resolve_pka_entropy_cutoff(
        cutoff_entropy_grimme, cutoff_entropy_truhlar
    )

    ctx.ensure_object(dict)
    ctx.obj["pka_shared"] = dict(
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        cutoff_entropy_grimme=s_freq_cutoff,
        cutoff_enthalpy=cutoff_enthalpy,
        entropy_method=entropy_method,
        scheme=scheme,
        delta_g_proton=delta_g_proton,
        pkb=pkb,
        pks=pks,
    )


@pka.command("analyze", cls=MyCommand)
@click_pka_analyze_options
@click.pass_context
def analyze(
    ctx,
    ha,
    a,
    href,
    ref,
    ha_solv,
    a_solv,
    href_solv,
    ref_solv,
    reference_pka,
    **kwargs,
):
    """Compute pKa from existing output files (auto-detects backend).

    Uses the Dual-level Proton Exchange scheme.

    \b
    Species labels:
        HA   - target acid           HRef - reference acid
        A-   - target conjugate base Ref- - reference conjugate base

    \b
    Reaction:  HA + Ref-  ->  A- + HRef

    \b
    Only -ha and -hr are required.  The remaining six files are
    auto-discovered from the same naming convention as batch-analyze:
      <basename>_pka_A_opt.<ext>   conjugate base
      <basename>_pka_HA_sp.<ext>   HA solvent single-point
      <basename>_pka_A_sp.<ext>    conjugate base solvent SP
      (and the corresponding _pka_Ref_* / _pka_HRef_sp files for HRef)
    If <basename>_pka_HA_opt_c*.<ext> ensemble files exist, those (and
    matching A/SP _c* files) are used instead of the single-file suffixes.
    Gas-phase and solvent outputs for each species must form a matching
    pair. N = 1 keeps the legacy names above, without a _c1 suffix.
    Override any auto-discovered path with the corresponding flag.

    \b
    Examples:
      chemsmart run pka analyze \\
          -ha acid_opt.log -hr collidine-H_opt.log -rp 6.75

      chemsmart run pka analyze \\
          -ha acid.log -a base.log \\
          -hr ref.log -r ref_base.log \\
          -has acid_sp.log -as base_sp.log \\
          -hrs ref_sp.log -rs ref_base_sp.log \\
          -rp 6.75 -T 298.15
    """
    shared = ctx.obj["pka_shared"]
    scheme = _resolve_pka_analysis_scheme(
        shared.get("scheme"), shared.get("delta_g_proton")
    )

    if scheme == "direct":
        if ha is None:
            raise click.UsageError(
                "-ha/--ha is required for direct-cycle pKa analysis."
            )
        optional = {"a": a, "ha_solv": ha_solv, "a_solv": a_solv}
        if any(v is None for v in optional.values()):
            discovered = _auto_discover_direct_pka_files(ha)
            for key in optional:
                if optional[key] is None:
                    optional[key] = discovered[key]
                    logger.info(f"Auto-discovered {key}: {discovered[key]}")
            a = optional["a"]
            ha_solv = optional["ha_solv"]
            a_solv = optional["a_solv"]
            if "ha" in discovered:
                ha = discovered["ha"]
                logger.info(f"Auto-discovered HA ensemble: {ha}")

        validate_direct_analyze_files(ha, a, ha_solv, a_solv)

        logger.info("Computing pKa (Direct Dissociation)...")
        print_pka_summary(
            ha_gas_file=ha,
            a_gas_file=a,
            ha_solv_file=ha_solv,
            a_solv_file=a_solv,
            scheme="direct",
            delta_G_proton=shared["delta_g_proton"],
            temperature=shared["temperature"],
            concentration=shared["concentration"],
            pressure=shared["pressure"],
            cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
            cutoff_enthalpy=shared["cutoff_enthalpy"],
            entropy_method=shared["entropy_method"],
            pkb=shared["pkb"],
            pks=shared["pks"],
        )
        return None

    # Auto-discover missing companion files when at least -ha and -hr given
    if ha is not None and href is not None:
        optional = {
            "a": a,
            "ha_solv": ha_solv,
            "a_solv": a_solv,
            "ref": ref,
            "href_solv": href_solv,
            "ref_solv": ref_solv,
        }
        if any(v is None for v in optional.values()):
            discovered = _auto_discover_pka_files(ha, href)
            for key in optional:
                if optional[key] is None:
                    optional[key] = discovered[key]
                    logger.info(f"Auto-discovered {key}: {discovered[key]}")
            a = optional["a"]
            ha_solv = optional["ha_solv"]
            a_solv = optional["a_solv"]
            ref = optional["ref"]
            href_solv = optional["href_solv"]
            ref_solv = optional["ref_solv"]
            if "ha" in discovered:
                ha = discovered["ha"]
                logger.info(f"Auto-discovered HA ensemble: {ha}")

    validate_analyze_files(
        ha, a, href, ref, ha_solv, a_solv, href_solv, ref_solv, reference_pka
    )

    logger.info(f"Computing pKa ({_scheme_display_name(scheme)})...")
    print_pka_summary(
        ha_gas_file=ha,
        a_gas_file=a,
        href_gas_file=href,
        ref_gas_file=ref,
        ha_solv_file=ha_solv,
        a_solv_file=a_solv,
        href_solv_file=href_solv,
        ref_solv_file=ref_solv,
        pka_reference=reference_pka,
        scheme=scheme,
        temperature=shared["temperature"],
        concentration=shared["concentration"],
        pressure=shared["pressure"],
        cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
        cutoff_enthalpy=shared["cutoff_enthalpy"],
        entropy_method=shared["entropy_method"],
        pkb=shared["pkb"],
        pks=shared["pks"],
    )

    # Return None so process_pipeline skips jobrunner execution.
    return None


@pka.command("batch-analyze", cls=MyCommand)
@click.option(
    "-o",
    "--output-table",
    type=click.Path(exists=True),
    required=True,
    help=(
        "Table of precomputed output file paths.  Columns: basename, ha_gas, a_gas, href_gas, ref_gas, ha_sp, a_sp, href_sp, ref_sp, pka_ref."
    ),
)
@click.option(
    "-O",
    "--output-results",
    type=click.Path(),
    default=None,
    help="Path to write the formatted results report.  Stdout if omitted.",
)
@click.option(
    "-p",
    "--program",
    type=click.Choice(["gaussian", "orca", "auto"]),
    default="auto",
    show_default=True,
    help=(
        "Require every populated output path to match this backend.  "
        "'auto' (default) parses each file independently and supports "
        "mixed Gaussian/ORCA tables."
    ),
)
@click.pass_context
def batch_analyze(ctx, output_table, output_results, program, **kwargs):
    """Batch pKa computation from a table of precomputed output files.

    Blank reference-acid cells are filled from the previous row.

    \b
    Examples:
      chemsmart run pka batch-analyze -o outputs.csv
      chemsmart run pka batch-analyze -o outputs.csv -O results.csv
    """
    shared = ctx.obj["pka_shared"]
    scheme = _resolve_pka_analysis_scheme(
        shared.get("scheme"), shared.get("delta_g_proton")
    )

    from chemsmart.utils.datasets import PKaOutputTable

    logger.info(f"Reading pKa output table: {output_table}")
    try:
        pka_output_table = PKaOutputTable.from_file(output_table)
        pka_output_table.prepare(check_file_exists=True, scheme=scheme)
    except (FileNotFoundError, ValueError) as e:
        raise click.UsageError(str(e))

    logger.info(f"Found {len(pka_output_table)} entries in output table")

    if program != "auto":
        _validate_pka_table_program(pka_output_table, program)

    program_label = "auto" if program == "auto" else program

    logger.info(
        f"Computing pKa ({_scheme_display_name(scheme)}) "
        f"for {len(pka_output_table)} systems "
        f"(T={shared['temperature']}K, program={program_label})"
    )
    results = pka_output_table.run_pka(
        output_cls=compute_pka,
        temperature=shared["temperature"],
        concentration=shared["concentration"],
        pressure=shared["pressure"],
        cutoff_entropy_grimme=shared["cutoff_entropy_grimme"],
        cutoff_enthalpy=shared["cutoff_enthalpy"],
        entropy_method=shared["entropy_method"],
        scheme=scheme,
        delta_G_proton=shared.get("delta_g_proton"),
    )
    output_string = pka_output_table.echo_pka_output_table_results(
        results=results,
        output_results=output_results,
        temperature=shared["temperature"],
        pressure=shared["pressure"],
        scheme=scheme,
        pkb=shared["pkb"],
        pks=shared["pks"],
    )
    click.echo(output_string)

    return None
