"""
Backend-independent pKa analysis and shared pKa sampling helpers.

``chemsmart run pka`` is post-processing and does not invoke Gaussian,
ORCA, or CREST. Gaussian and ORCA pKa jobs import the sampling helpers
defined here.

Subcommands
-----------
analyze        Compute pKa from 8 existing output files.
batch-analyze  Batch pKa from a table of output file paths.
"""

import functools
import logging
import math
import os

import click

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.cli.thermochemistry.thermochemistry import (
    resolve_entropy_cutoff,
    thermochemistry_cutoff_options,
    thermochemistry_temp_pressure_conc_options,
)
from chemsmart.io.file import PKaCDXFile
from chemsmart.utils.cli import MyCommand, MyGroup
from chemsmart.utils.constants import (
    HARTREE_TO_KCAL_MOL,
    R,
    atm_to_pa,
    energy_conversion,
)
from chemsmart.utils.io import get_program_type_from_file

logger = logging.getLogger(__name__)


def resolve_pka_entropy_cutoff(cutoff_entropy_grimme, cutoff_entropy_truhlar):
    """Resolve pKa entropy cutoff; default to Grimme 100 cm⁻¹ when unset."""
    s_freq_cutoff, entropy_method = resolve_entropy_cutoff(
        cutoff_entropy_grimme, cutoff_entropy_truhlar
    )
    if s_freq_cutoff is None:
        return 100.0, "grimme"
    return s_freq_cutoff, entropy_method


def _pka_thermochemistry_kwargs(
    temperature,
    concentration,
    pressure,
    cutoff_entropy_grimme,
    cutoff_enthalpy,
    entropy_method,
):
    return {
        "temperature": temperature,
        "concentration": concentration,
        "pressure": pressure,
        "s_freq_cutoff": cutoff_entropy_grimme,
        "h_freq_cutoff": cutoff_enthalpy,
        "entropy_method": entropy_method,
        "energy_units": "hartree",
        "check_imaginary_frequencies": True,
    }


def _extract_thermochemistry_property(thermo, filepath, attr, label):
    value = getattr(thermo, attr)
    if value is None:
        raise ValueError(f"Could not extract {label} from file: {filepath}")
    return value


def pka_gas_phase_data(
    filepath,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    entropy_method="grimme",
):
    """Return gas-phase SCF energy and qh-G correction in Hartree."""
    thermo = Thermochemistry(
        filename=filepath,
        **_pka_thermochemistry_kwargs(
            temperature,
            concentration,
            pressure,
            cutoff_entropy_grimme,
            cutoff_enthalpy,
            entropy_method,
        ),
    )
    electronic_energy_j_mol = _extract_thermochemistry_property(
        thermo,
        filepath,
        "electronic_energy",
        "SCF energy",
    )
    qh_gibbs_j_mol = _extract_thermochemistry_property(
        thermo,
        filepath,
        "qrrho_gibbs_free_energy",
        "quasi-harmonic Gibbs free energy",
    )
    electronic_energy_au = energy_conversion(
        "j/mol", "hartree", electronic_energy_j_mol
    )
    qh_gibbs_au = energy_conversion("j/mol", "hartree", qh_gibbs_j_mol)
    return electronic_energy_au, qh_gibbs_au - electronic_energy_au


def pka_solvent_scf_energy(filepath):
    """Return solvent-phase SCF energy in Hartree."""
    thermo = Thermochemistry(filename=filepath)
    electronic_energy_j_mol = _extract_thermochemistry_property(
        thermo,
        filepath,
        "electronic_energy",
        "SCF energy",
    )
    return energy_conversion("j/mol", "hartree", electronic_energy_j_mol)


# Kelly, Cramer, Truhlar ΔG*_solv(H+) at 298 K in water; not T-corrected.
# J. Phys. Chem. B 2006, 110, 16066-16081.
KELLY_PROTON_SOLVATION_FREE_ENERGY_KCAL_MOL = -265.9
# Sackur–Tetrode / JANAF S°(H+, 1 atm).
PROTON_GAS_STANDARD_ENTROPY_CAL_MOL_K = 26.016
_AQUEOUS_SOLVENT_IDS = frozenset({"water", "h2o"})
DEFAULT_PKS = 14.0


def aqueous_proton_solution_free_energy_kcal_mol(
    temperature,
    delta_g_solv=KELLY_PROTON_SOLVATION_FREE_ENERGY_KCAL_MOL,
):
    """Return aqueous G*_soln(H+) in kcal/mol at *temperature* (K).

    G*_aq(H+) = G°_gas(H+, 1 atm) + RT ln(RT / P°) + ΔG*_solv(H+)

    G°_gas(H+) = 5/2 RT − T S°(H+, 1 atm), with S°(H+) = 26.016 cal mol⁻¹ K⁻¹.
    ΔG*_solv(H+) defaults to the Kelly, Cramer, and Truhlar aqueous value
    (−265.9 kcal/mol; J. Phys. Chem. B 2006, 110, 16066) and is not
    temperature-corrected. At 298.15 K the result is ≈ −270.3 kcal/mol.
    """
    r_kcal_mol_k = energy_conversion("j/mol", "kcal/mol", R)
    s_kcal_mol_k = PROTON_GAS_STANDARD_ENTROPY_CAL_MOL_K / 1000.0
    g_gas_kcal_mol = (
        2.5 * r_kcal_mol_k * temperature - temperature * s_kcal_mol_k
    )
    r_liter_atm_mol_k = R / atm_to_pa * 1000.0
    standard_molar_volume_l = r_liter_atm_mol_k * temperature
    g_standard_state_kcal_mol = (
        r_kcal_mol_k * temperature * math.log(standard_molar_volume_l)
    )
    return g_gas_kcal_mol + g_standard_state_kcal_mol + delta_g_solv


def warn_if_non_aqueous_direct_proton_default(
    scheme, delta_g_proton, solvent_id
):
    """Warn if aqueous G_soln(H+) default is used with a non-water solvent."""
    if scheme != "direct" or delta_g_proton is not None:
        return
    if solvent_id is None:
        return
    if str(solvent_id).strip().lower() in _AQUEOUS_SOLVENT_IDS:
        return
    logger.warning(
        "Computed default G_soln(H+) is for aqueous water; "
        "solvent_id=%r is not water. Pass -dG/--delta-g-proton for a "
        "literature non-aqueous G_soln(H+) value.",
        solvent_id,
    )


def pks_to_pkb(pka, pks):
    """Return pKb = pKs − pKa."""
    return pks - pka


def resolve_pkb_reporting(pkb=False, pks=None):
    """Return ``(report_pkb, pks, pks_defaulted)`` for pKb conversion.

    Reporting is enabled when ``pkb`` is true or ``pks`` is supplied.
    If reporting is enabled and ``pks`` is omitted, ``DEFAULT_PKS`` is used.
    """
    report = bool(pkb) or pks is not None
    if not report:
        return False, None, False
    if pks is None:
        return True, DEFAULT_PKS, True
    return True, pks, False


def warn_if_default_pks_non_aqueous(pks_defaulted, solvent_id):
    """Warn if default aqueous pKs=14 is used with a non-water solvent."""
    if not pks_defaulted or solvent_id is None:
        return
    if str(solvent_id).strip().lower() in _AQUEOUS_SOLVENT_IDS:
        return
    logger.warning(
        "Default pKs = 14.0 is for aqueous water; "
        "solvent_id=%r is not water. Pass --pks for a literature "
        "non-aqueous autoprotolysis constant.",
        solvent_id,
    )


def ensemble_effective_free_energy(g_values, temperature):
    """Return the ensemble effective free energy in Hartree.

    ``G_eff = -RT ln Σ exp(-G_i / RT)``, evaluated with a log-sum-exp
    shift for numerical stability. ``g_values`` and the return value are
    in Hartree; ``temperature`` is in Kelvin.
    """
    if temperature is None or temperature <= 0:
        raise ValueError("temperature must be a positive value in Kelvin.")
    energies = [float(g) for g in g_values]
    if not energies:
        raise ValueError("g_values must contain at least one free energy.")
    if len(energies) == 1:
        return energies[0]
    rt_hartree = energy_conversion("j/mol", "hartree", R * temperature)
    g_min = min(energies)
    log_sum = math.log(
        sum(math.exp(-(g - g_min) / rt_hartree) for g in energies)
    )
    return g_min - rt_hartree * log_sum


def _normalize_pka_file_arg(value):
    """Return a non-empty list of filesystem paths, or ``None``."""
    if value is None:
        return None
    if isinstance(value, (str, bytes, os.PathLike)):
        return [os.fspath(value)]
    try:
        paths = [os.fspath(path) for path in value]
    except TypeError as exc:
        raise TypeError(
            f"Expected a path or sequence of paths, got {type(value)!r}."
        ) from exc
    if not paths:
        raise ValueError("File list must be non-empty.")
    return paths


def _species_solution_free_energy(
    gas_files, solv_files, thermo_kwargs, temperature, label
):
    """Return G_soln (or G_eff), component energies, and conformer count."""
    gas_files = _normalize_pka_file_arg(gas_files)
    solv_files = _normalize_pka_file_arg(solv_files)
    if gas_files is None or solv_files is None:
        raise ValueError(f"Missing required files for {label}.")
    if len(gas_files) != len(solv_files):
        raise ValueError(
            f"{label} gas-phase and solvent file counts must match "
            f"({len(gas_files)} vs {len(solv_files)})."
        )
    e_gas_values = []
    g_corr_values = []
    e_solv_values = []
    g_soln_values = []
    for gas_file, solv_file in zip(gas_files, solv_files):
        e_gas, g_corr = pka_gas_phase_data(gas_file, **thermo_kwargs)
        e_solv = pka_solvent_scf_energy(solv_file)
        e_gas_values.append(e_gas)
        g_corr_values.append(g_corr)
        e_solv_values.append(e_solv)
        g_soln_values.append(e_solv + g_corr)
    n_conformers = len(g_soln_values)
    g_soln = (
        ensemble_effective_free_energy(g_soln_values, temperature)
        if n_conformers > 1
        else g_soln_values[0]
    )
    return {
        "E_gas": e_gas_values[0],
        "G_corr": g_corr_values[0],
        "E_solv": e_solv_values[0],
        "G_soln": g_soln,
        "num_conformers": n_conformers,
    }


def compute_pka(
    ha_gas_file,
    a_gas_file,
    href_gas_file=None,
    ref_gas_file=None,
    ha_solv_file=None,
    a_solv_file=None,
    href_solv_file=None,
    ref_solv_file=None,
    pka_reference=None,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    entropy_method="grimme",
    scheme="proton exchange",
    delta_G_proton=None,
):
    """Compute pKa from output files using program-independent thermochemistry.

    Each species file argument may be a path or a sequence of paths. For a
    species with more than one conformer, ``G_soln,i = E_solv,i + G_corr,i``
    is computed per conformer and replaced by the ensemble effective free
    energy ``G_eff = -RT ln Σ exp(-G_i / RT)``.

    For ``scheme='direct'``, ``delta_G_proton`` is G_soln(H+) in kcal/mol.
    If omitted, a T-dependent aqueous default is computed from Kelly,
    Cramer, and Truhlar ΔG*_solv(H+) = -265.9 kcal/mol.
    """
    proton_user_supplied = delta_G_proton is not None
    if scheme == "direct":
        if ha_solv_file is None or a_solv_file is None:
            raise ValueError(
                "ha_solv_file and a_solv_file are required for scheme='direct'."
            )
        if not proton_user_supplied:
            delta_G_proton = aqueous_proton_solution_free_energy_kcal_mol(
                temperature
            )
    elif pka_reference is None:
        raise ValueError(
            "pka_reference is required when scheme='proton exchange'."
        )
    else:
        missing = [
            name
            for name, value in (
                ("href_gas_file", href_gas_file),
                ("ref_gas_file", ref_gas_file),
                ("ha_solv_file", ha_solv_file),
                ("a_solv_file", a_solv_file),
                ("href_solv_file", href_solv_file),
                ("ref_solv_file", ref_solv_file),
            )
            if value is None
        ]
        if missing:
            raise ValueError(
                "Missing required files for proton exchange scheme: "
                + ", ".join(missing)
            )

    thermo_kwargs = dict(
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
        entropy_method=entropy_method,
    )

    ha_data = _species_solution_free_energy(
        ha_gas_file, ha_solv_file, thermo_kwargs, temperature, "HA"
    )
    a_data = _species_solution_free_energy(
        a_gas_file, a_solv_file, thermo_kwargs, temperature, "A-"
    )
    G_soln_HA_au = ha_data["G_soln"]
    G_soln_A_au = a_data["G_soln"]
    E_gas_HA_au = ha_data["E_gas"]
    E_gas_A_au = a_data["E_gas"]
    G_corr_HA_au = ha_data["G_corr"]
    G_corr_A_au = a_data["G_corr"]
    E_solv_HA_au = ha_data["E_solv"]
    E_solv_A_au = a_data["E_solv"]

    R_kcal = 0.001987204
    ln10 = 2.302585093

    if scheme == "direct":
        G_soln_HA_kcal = G_soln_HA_au * HARTREE_TO_KCAL_MOL
        G_soln_A_kcal = G_soln_A_au * HARTREE_TO_KCAL_MOL
        delta_G_diss_kcal_mol = G_soln_A_kcal + delta_G_proton - G_soln_HA_kcal
        delta_G_diss_au = delta_G_diss_kcal_mol / HARTREE_TO_KCAL_MOL
        pka = delta_G_diss_kcal_mol / (R_kcal * temperature * ln10)
        return {
            "pKa": pka,
            "scheme": "direct",
            "delta_G_proton_kcal_mol": delta_G_proton,
            "delta_G_proton_user_supplied": proton_user_supplied,
            "delta_G_diss_kcal_mol": delta_G_diss_kcal_mol,
            "delta_G_diss_au": delta_G_diss_au,
            "delta_G_soln_kcal_mol": delta_G_diss_kcal_mol,
            "delta_G_soln_au": delta_G_diss_au,
            "temperature": temperature,
            "G_soln_HA_au": G_soln_HA_au,
            "G_soln_A_au": G_soln_A_au,
            "E_solv_HA_au": E_solv_HA_au,
            "E_solv_A_au": E_solv_A_au,
            "G_corr_HA_au": G_corr_HA_au,
            "G_corr_A_au": G_corr_A_au,
            "E_gas_HA_au": E_gas_HA_au,
            "E_gas_A_au": E_gas_A_au,
            "num_conformers_HA": ha_data["num_conformers"],
            "num_conformers_A": a_data["num_conformers"],
        }

    href_data = _species_solution_free_energy(
        href_gas_file, href_solv_file, thermo_kwargs, temperature, "HRef"
    )
    ref_data = _species_solution_free_energy(
        ref_gas_file, ref_solv_file, thermo_kwargs, temperature, "Ref-"
    )
    G_soln_HRef_au = href_data["G_soln"]
    G_soln_Ref_au = ref_data["G_soln"]
    E_gas_HRef_au = href_data["E_gas"]
    E_gas_Ref_au = ref_data["E_gas"]
    G_corr_HRef_au = href_data["G_corr"]
    G_corr_Ref_au = ref_data["G_corr"]
    E_solv_HRef_au = href_data["E_solv"]
    E_solv_Ref_au = ref_data["E_solv"]

    delta_G_soln_au = (G_soln_A_au + G_soln_HRef_au) - (
        G_soln_HA_au + G_soln_Ref_au
    )
    delta_G_soln_kcal_mol = delta_G_soln_au * HARTREE_TO_KCAL_MOL
    pka = pka_reference + delta_G_soln_kcal_mol / (R_kcal * temperature * ln10)

    return {
        "pKa": pka,
        "scheme": "proton exchange",
        "pKa_reference": pka_reference,
        "delta_G_soln_kcal_mol": delta_G_soln_kcal_mol,
        "delta_G_soln_au": delta_G_soln_au,
        "temperature": temperature,
        "G_soln_HA_au": G_soln_HA_au,
        "G_soln_A_au": G_soln_A_au,
        "G_soln_HRef_au": G_soln_HRef_au,
        "G_soln_Ref_au": G_soln_Ref_au,
        "E_solv_HA_au": E_solv_HA_au,
        "E_solv_A_au": E_solv_A_au,
        "E_solv_HRef_au": E_solv_HRef_au,
        "E_solv_Ref_au": E_solv_Ref_au,
        "G_corr_HA_au": G_corr_HA_au,
        "G_corr_A_au": G_corr_A_au,
        "G_corr_HRef_au": G_corr_HRef_au,
        "G_corr_Ref_au": G_corr_Ref_au,
        "E_gas_HA_au": E_gas_HA_au,
        "E_gas_A_au": E_gas_A_au,
        "E_gas_HRef_au": E_gas_HRef_au,
        "E_gas_Ref_au": E_gas_Ref_au,
        "num_conformers_HA": ha_data["num_conformers"],
        "num_conformers_A": a_data["num_conformers"],
        "num_conformers_HRef": href_data["num_conformers"],
        "num_conformers_Ref": ref_data["num_conformers"],
    }


def _thermochemistry_value_in_units(
    thermo, filepath, attr, energy_units, label
):
    value_j_mol = _extract_thermochemistry_property(
        thermo, filepath, attr, label
    )
    return energy_conversion("j/mol", energy_units, value_j_mol)


def compute_pka_thermochemistry(
    ha_file=None,
    a_file=None,
    href_file=None,
    ref_file=None,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    energy_units="hartree",
    entropy_method="grimme",
):
    """Extract gas-phase thermochemistry for pKa species from output files.

    Each species argument may be a path or a sequence of paths. A single
    path, or a one-element sequence, returns one record. A longer sequence
    returns one record per conformer, in the given order.
    """
    results = {
        "settings": {
            "temperature": temperature,
            "concentration": concentration,
            "pressure": pressure,
            "cutoff_entropy_grimme": cutoff_entropy_grimme,
            "cutoff_enthalpy": cutoff_enthalpy,
            "energy_units": energy_units,
        }
    }
    thermo_kwargs = _pka_thermochemistry_kwargs(
        temperature,
        concentration,
        pressure,
        cutoff_entropy_grimme,
        cutoff_enthalpy,
        entropy_method,
    )

    def get_species_thermo(filepath, name):
        paths = _normalize_pka_file_arg(filepath)
        if paths is None:
            return None
        records = [_species_thermo_record(path, name) for path in paths]
        if len(records) == 1:
            return records[0]
        return records

    def _species_thermo_record(filepath, name):
        thermo = Thermochemistry(filename=filepath, **thermo_kwargs)
        return {
            "name": name,
            "E": _thermochemistry_value_in_units(
                thermo,
                filepath,
                "electronic_energy",
                energy_units,
                "SCF energy",
            ),
            "qh_G": _thermochemistry_value_in_units(
                thermo,
                filepath,
                "qrrho_gibbs_free_energy",
                energy_units,
                "quasi-harmonic Gibbs free energy",
            ),
            "ZPE": _thermochemistry_value_in_units(
                thermo,
                filepath,
                "zero_point_energy",
                energy_units,
                "zero-point energy",
            ),
            "H": _thermochemistry_value_in_units(
                thermo, filepath, "enthalpy", energy_units, "enthalpy"
            ),
            "qh_H": _thermochemistry_value_in_units(
                thermo,
                filepath,
                "qrrho_enthalpy",
                energy_units,
                "quasi-harmonic enthalpy",
            ),
            "G": _thermochemistry_value_in_units(
                thermo,
                filepath,
                "gibbs_free_energy",
                energy_units,
                "Gibbs free energy",
            ),
        }

    if ha_file is not None:
        results["HA"] = get_species_thermo(ha_file, "HA")
    if a_file is not None:
        results["A"] = get_species_thermo(a_file, "A-")
    if href_file is not None:
        results["HRef"] = get_species_thermo(href_file, "HRef")
    if ref_file is not None:
        results["Ref"] = get_species_thermo(ref_file, "Ref-")
    return results


def _print_computed_pka_pkb(pka, pkb=False, pks=None, solvent_id=None):
    """Print the computed pKa line and optional pKb conversion lines."""
    report_pkb, pks_value, pks_defaulted = resolve_pkb_reporting(
        pkb=pkb, pks=pks
    )
    warn_if_default_pks_non_aqueous(pks_defaulted, solvent_id)
    if report_pkb:
        source = "default aqueous" if pks_defaulted else "user-supplied"
        print(f"  pKs = {pks_value:.2f} ({source})")
        print()
        print(f"  *** Computed pKa(HA) = {pka:.2f} ***")
        print(f"  *** Computed pKb(B)  = {pks_to_pkb(pka, pks_value):.2f} ***")
        return
    print(f"  *** Computed pKa(HA) = {pka:.2f} ***")


def _pka_result_uses_ensemble(result):
    keys = (
        "num_conformers_HA",
        "num_conformers_A",
        "num_conformers_HRef",
        "num_conformers_Ref",
    )
    return any(result.get(key, 1) > 1 for key in keys)


def _format_pka_g_soln_line(label, value, num_conformers):
    if num_conformers is not None and num_conformers > 1:
        return (
            f"  {label} ({num_conformers} conformers, G_eff):  "
            f"{value:.10f}"
        )
    return f"  {label}:  {value:.10f}"


def _print_pka_g_soln_method_lines(result):
    print("  G_soln = E_solv + G_corr  (solution free energy)")
    if _pka_result_uses_ensemble(result):
        print(
            "  For multiple conformers, G_soln is replaced by "
            "G_eff = -RT ln Σ exp(-G_i/RT)"
        )


def print_pka_summary(
    ha_gas_file,
    a_gas_file,
    href_gas_file=None,
    ref_gas_file=None,
    ha_solv_file=None,
    a_solv_file=None,
    href_solv_file=None,
    ref_solv_file=None,
    pka_reference=None,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    entropy_method="grimme",
    scheme="proton exchange",
    delta_G_proton=None,
    pkb=False,
    pks=None,
    solvent_id=None,
):
    """Print a formatted summary of a dual-level pKa calculation."""
    result = compute_pka(
        ha_gas_file=ha_gas_file,
        a_gas_file=a_gas_file,
        href_gas_file=href_gas_file,
        ref_gas_file=ref_gas_file,
        ha_solv_file=ha_solv_file,
        a_solv_file=a_solv_file,
        href_solv_file=href_solv_file,
        ref_solv_file=ref_solv_file,
        pka_reference=pka_reference,
        temperature=temperature,
        concentration=concentration,
        pressure=pressure,
        cutoff_entropy_grimme=cutoff_entropy_grimme,
        cutoff_enthalpy=cutoff_enthalpy,
        entropy_method=entropy_method,
        scheme=scheme,
        delta_G_proton=delta_G_proton,
    )

    if scheme == "direct":
        print("=" * 78)
        print("pKa Calculation - Direct Dissociation Scheme")
        print("=" * 78)
        print("Reaction: HA → A⁻ + H⁺")
        print(f"Temperature: {temperature} K")
        print()
        print("Method:")
        print("  G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)")
        _print_pka_g_soln_method_lines(result)
        print("  ΔG_diss = G_soln(A⁻) + G_soln(H⁺) - G_soln(HA)")
        print("  pKa = ΔG_diss / (2.303 × R × T)")
        print("-" * 78)
        print()
        print("Gas-Phase Electronic Energies (E_gas, au):")
        print(f"  HA:  {result['E_gas_HA_au']:.10f}")
        print(f"  A⁻:  {result['E_gas_A_au']:.10f}")
        print()
        print("Thermal Corrections (G_corr = qh-G - E_gas, au):")
        print(f"  HA:  {result['G_corr_HA_au']:.10f}")
        print(f"  A⁻:  {result['G_corr_A_au']:.10f}")
        print()
        print("Solvent Single-Point Energies (E_solv, au):")
        print(f"  HA:  {result['E_solv_HA_au']:.10f}")
        print(f"  A⁻:  {result['E_solv_A_au']:.10f}")
        print()
        print("Solution Free Energies (G_soln = E_solv + G_corr, au):")
        print(
            _format_pka_g_soln_line(
                "HA",
                result["G_soln_HA_au"],
                result.get("num_conformers_HA", 1),
            )
        )
        print(
            _format_pka_g_soln_line(
                "A⁻",
                result["G_soln_A_au"],
                result.get("num_conformers_A", 1),
            )
        )
        print("-" * 78)
        print()
        print("pKa Calculation:")
        g_soln_h = result["delta_G_proton_kcal_mol"]
        print(f"  G_soln(H⁺) = {g_soln_h:.4f} kcal/mol")
        if result["delta_G_proton_user_supplied"]:
            print("             (user-supplied)")
        else:
            print(
                "             (computed aqueous default for water at "
                f"{temperature} K)"
            )
        print(f"  ΔG_diss = {result['delta_G_diss_au']:.10f} au")
        print(f"         = {result['delta_G_diss_kcal_mol']:.4f} kcal/mol")
        print()
        _print_computed_pka_pkb(
            result["pKa"], pkb=pkb, pks=pks, solvent_id=solvent_id
        )
        print("=" * 78)
        return result

    print("=" * 78)
    print("pKa Calculation - Dual-level Proton Exchange Scheme")
    print("=" * 78)
    print("Reaction: HA + Ref⁻ → A⁻ + HRef")
    print(f"Temperature: {temperature} K")
    print()
    print("Method:")
    print("  G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)")
    _print_pka_g_soln_method_lines(result)
    print(
        "  ΔG_soln = [G(A⁻)_soln + G(HRef)_soln] - [G(HA)_soln + G(Ref⁻)_soln]"
    )
    print("  pKa = pKa_ref + ΔG_soln / (RT × ln10)")
    print("-" * 78)
    print()
    print("Gas-Phase Electronic Energies (E_gas, au):")
    print(f"  HA:  {result['E_gas_HA_au']:.10f}")
    print(f"  A⁻:  {result['E_gas_A_au']:.10f}")
    print(f"  HRef:  {result['E_gas_HRef_au']:.10f}")
    print(f"  Ref⁻:  {result['E_gas_Ref_au']:.10f}")
    print()
    print("Thermal Corrections (G_corr = qh-G - E_gas, au):")
    print(f"  HA:  {result['G_corr_HA_au']:.10f}")
    print(f"  A⁻:  {result['G_corr_A_au']:.10f}")
    print(f"  HRef:  {result['G_corr_HRef_au']:.10f}")
    print(f"  Ref⁻:  {result['G_corr_Ref_au']:.10f}")
    print()
    print("Solvent Single-Point Energies (E_solv, au):")
    print(f"  HA:  {result['E_solv_HA_au']:.10f}")
    print(f"  A⁻:  {result['E_solv_A_au']:.10f}")
    print(f"  HRef:  {result['E_solv_HRef_au']:.10f}")
    print(f"  Ref⁻:  {result['E_solv_Ref_au']:.10f}")
    print()
    print("Solution Free Energies (G_soln = E_solv + G_corr, au):")
    print(
        _format_pka_g_soln_line(
            "HA",
            result["G_soln_HA_au"],
            result.get("num_conformers_HA", 1),
        )
    )
    print(
        _format_pka_g_soln_line(
            "A⁻",
            result["G_soln_A_au"],
            result.get("num_conformers_A", 1),
        )
    )
    print(
        _format_pka_g_soln_line(
            "HRef",
            result["G_soln_HRef_au"],
            result.get("num_conformers_HRef", 1),
        )
    )
    print(
        _format_pka_g_soln_line(
            "Ref⁻",
            result["G_soln_Ref_au"],
            result.get("num_conformers_Ref", 1),
        )
    )
    print("-" * 78)
    print()
    print("pKa Calculation:")
    print(f"  ΔG_soln = {result['delta_G_soln_au']:.10f} au")
    print(f"         = {result['delta_G_soln_kcal_mol']:.4f} kcal/mol")
    print(f"  pKa(HRef)_ref = {pka_reference:.2f}")
    print()
    _print_computed_pka_pkb(
        result["pKa"], pkb=pkb, pks=pks, solvent_id=solvent_id
    )
    print("=" * 78)
    return result


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


def build_pka_crest_job(molecule, label, settings, parent_job):
    """Build a CREST conformer-search job for one pKa species.

    CREST settings are loaded from ``settings.crest_project`` when that YAML
    exists. Otherwise default conformer settings are used. Charge and
    multiplicity come from ``molecule``; ``nprocs`` comes from the parent
    jobrunner.

    Args:
        molecule: HA or A- ``Molecule`` to sample.
        label (str): CREST job label (for example ``{label}_HA_crest``).
        settings: pKa job settings exposing ``crest_project``.
        parent_job: Parent Gaussian or ORCA pKa job.

    Returns:
        CRESTConformerSearchJob: Job with a typed CREST runner.
    """
    from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
    from chemsmart.jobs.runner import JobRunner

    parent_runner = parent_job.jobrunner
    crest_settings = _pka_crest_job_settings(molecule, settings, parent_runner)
    job = CRESTConformerSearchJob(
        molecule=molecule,
        settings=crest_settings,
        label=label,
        jobrunner=None,
        skip_completed=parent_job.skip_completed,
    )
    job.jobrunner = JobRunner.from_job(
        job,
        server=parent_runner.server,
        fake=parent_runner.FAKE,
        num_cores=parent_runner.num_cores,
    )
    return job


def select_crest_conformers(crest_job, num_conformers, fallback_molecule):
    """Select CREST geometries for DFT, or fall back to the input molecule.

    If CREST has not written an output file yet, return ``None`` so the
    caller can halt and wait (HPC resubmit). After CREST has terminated
    (normally or abnormally), always return a non-empty list: selected
    conformers, a shorter available set, or ``[fallback_molecule]``.

    ``N == 1`` uses ``crest_best.xyz`` (else the first frame of
    ``crest_conformers.xyz``). ``N > 1`` uses the N lowest frames of
    energy-sorted ``crest_conformers.xyz``. Charge and multiplicity are
    copied from ``fallback_molecule``.

    Args:
        crest_job: CREST job whose folder is parsed.
        num_conformers (int): Number of conformers requested (``N >= 1``).
        fallback_molecule: Input HA or A- geometry used if extraction fails.

    Returns:
        list or None: Selected molecules, or ``None`` if CREST has not
        finished.
    """
    if num_conformers is None or num_conformers < 1:
        num_conformers = 1

    try:
        if not os.path.exists(crest_job.outputfile):
            logger.info(
                f"CREST job {crest_job.label} has not finished; waiting."
            )
            return None
        return _select_pka_crest_conformers(
            crest_job, num_conformers, fallback_molecule
        )
    except Exception as exc:
        logger.warning(
            f"Failed to extract CREST conformers from {crest_job.label}: "
            f"{exc}. Using the input geometry."
        )
        return [_copy_pka_molecule_charge(fallback_molecule)]


def _default_pka_crest_settings():
    from chemsmart.jobs.crest.settings import CRESTJobSettings

    settings = CRESTJobSettings.default()
    settings.jobtype = "conformers"
    return settings


def _pka_crest_job_settings(molecule, settings, parent_runner):
    from chemsmart.settings.crest import CRESTProjectSettings

    crest_project = settings.crest_project
    crest_settings = None
    if crest_project is not None:
        try:
            crest_settings = CRESTProjectSettings.from_project(
                crest_project
            ).conformer_settings()
        except FileNotFoundError:
            logger.warning(
                f"No CREST project settings found for {crest_project!r}; "
                "using default CREST conformer settings."
            )

    if crest_settings is None:
        crest_settings = _default_pka_crest_settings()
    else:
        crest_settings = crest_settings.copy()
        if crest_settings.jobtype is None:
            crest_settings.jobtype = "conformers"

    crest_settings.charge = molecule.charge
    crest_settings.multiplicity = molecule.multiplicity
    if parent_runner is not None and parent_runner.num_cores is not None:
        crest_settings.nprocs = parent_runner.num_cores
    return crest_settings


def _select_pka_crest_conformers(crest_job, num_conformers, fallback_molecule):
    from chemsmart.io.crest.output import CRESTOutput

    output = CRESTOutput(folder=crest_job.folder)
    if not output.normal_termination:
        logger.warning(
            f"CREST job {crest_job.label} did not terminate normally. "
            "Using available geometries or the input structure."
        )

    if num_conformers == 1:
        selected = output.best_conformer
        if selected is None:
            logger.warning(
                f"CREST job {crest_job.label} produced no usable geometry. "
                "Using the input structure."
            )
            return [_copy_pka_molecule_charge(fallback_molecule)]
        return [_copy_pka_molecule_charge(selected, fallback_molecule)]

    conformers = list(output.conformers)
    if not conformers and output.best_conformer is not None:
        conformers = [output.best_conformer]

    if not conformers:
        logger.warning(
            f"CREST job {crest_job.label} produced no conformers. "
            "Using the input structure."
        )
        return [_copy_pka_molecule_charge(fallback_molecule)]

    if len(conformers) < num_conformers:
        logger.warning(
            f"CREST job {crest_job.label} produced {len(conformers)} "
            f"conformer(s); requested {num_conformers}. Using the available "
            "set."
        )

    return [
        _copy_pka_molecule_charge(molecule, fallback_molecule)
        for molecule in conformers[:num_conformers]
    ]


def _copy_pka_molecule_charge(molecule, charge_source=None):
    copied = molecule.copy()
    source = charge_source if charge_source is not None else molecule
    copied.charge = source.charge
    copied.multiplicity = source.multiplicity
    return copied


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


def _resolve_site_from_molecule_file(filename, mode):
    from chemsmart.io.molecules.structure import Molecule

    molecule = Molecule.from_filepath(filename)
    if molecule is None:
        raise ValueError(
            f"Could not read a molecule from {filename}. "
            "Specify -pi/--proton-index or -cc/--color-code."
        )
    return resolve_ionizable_site(molecule, mode=mode), None


def _resolve_cdxml_smarts_sites(filename, mode):
    """Resolve a unique SMARTS site for each CDXML fragment."""
    from chemsmart.io.file import PKaCDXFile

    molecules = list(PKaCDXFile(filename).molecules)
    if not molecules:
        raise ValueError(
            f"Could not read a molecule from {filename}. "
            "Specify -pi/--proton-index or -cc/--color-code."
        )
    prepared = []
    for molecule in molecules:
        site = resolve_ionizable_site(molecule, mode=mode)
        molecule.proton_index = site
        prepared.append(molecule)
    if len(prepared) > 1:
        return None, prepared
    return prepared[0].proton_index, None


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
    if is_pka_cdxml_input(filename):
        from chemsmart.io.file import PKaCDXFile

        try:
            return PKaCDXFile(filename)._resolve_proton_from_cdxml(
                color_code, mode=mode
            )
        except ValueError as exc:
            if color_code is not None or not _cdxml_has_no_colour_markup(exc):
                raise
            return _resolve_cdxml_smarts_sites(filename, mode)

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

    return _resolve_site_from_molecule_file(filename, mode)


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


def resolve_pka_batch_row(
    filepath, proton_index=None, color_code=None, pkb=False
):
    """Resolve site index and molecule for one pKa submission-table row.

    Each table row maps to a single job. An explicit ``proton_index`` always
    takes precedence (acidic H, or with ``pkb`` the heavy atom to protonate).
    When omitted, a unique SMARTS site is used for non-CDXML files; for a
    single-molecule ``.cdxml`` / ``.cdx`` filepath the coloured site is
    auto-detected, falling through to SMARTS when there is no colour markup.
    Multi-molecule CDXML files are rejected here; pass them directly as ``-f``
    with ``pka batch`` instead.

    Returns:
        tuple[int, Molecule | PKaMolecule]: Resolved site index and structure.
    """
    from chemsmart.io.file import PKaCDXFile
    from chemsmart.io.molecules.structure import Molecule

    filepath = str(filepath)
    mode = pka_submit_site_mode(pkb)
    if proton_index is not None:
        return int(proton_index), Molecule.from_filepath(filepath)

    if not is_pka_cdxml_input(filepath):
        molecule = Molecule.from_filepath(filepath)
        if molecule is None:
            raise ValueError(
                f"Could not read a molecule from {filepath}. "
                "Provide proton_index in the table, or use a single-molecule "
                ".cdxml/.cdx file with a coloured site and leave proton_index "
                "blank."
            )
        return resolve_ionizable_site(molecule, mode=mode), molecule

    if pkb:
        site, molecules = resolve_proton_index(
            filepath, None, color_code, mode="base"
        )
        if molecules is not None:
            if len(molecules) != 1:
                raise ValueError(
                    f"CDXML file {filepath} contains {len(molecules)} "
                    "molecules. Submission-table rows support "
                    "single-molecule CDXML files only. Pass a "
                    "multi-molecule CDXML file directly as -f with pka batch."
                )
            return molecules[0].proton_index, molecules[0]
        molecule = Molecule.from_filepath(filepath)
        if molecule is None:
            raise ValueError(f"Could not read a molecule from {filepath}.")
        return site, molecule

    cdx_file = PKaCDXFile(filepath)
    try:
        pka_molecules = cdx_file.get_pka_molecules(
            color_code=color_code,
            index=":",
            return_list=True,
        )
    except ValueError as exc:
        if not _cdxml_has_no_colour_markup(exc):
            raise ValueError(
                f"Could not auto-detect proton from CDXML colour for {filepath}: "
                f"{exc}"
            ) from exc
        molecule = Molecule.from_filepath(filepath)
        if molecule is None:
            raise ValueError(
                f"Could not read a molecule from {filepath}."
            ) from exc
        return resolve_ionizable_site(molecule, mode=mode), molecule

    if len(pka_molecules) != 1:
        raise ValueError(
            f"CDXML file {filepath} contains {len(pka_molecules)} molecules. "
            "Submission-table rows support single-molecule CDXML files only. "
            "Pass a multi-molecule CDXML file directly as -f with pka batch."
        )

    pka_mol = pka_molecules[0]
    return pka_mol.proton_index, pka_mol


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
