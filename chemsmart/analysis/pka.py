"""Program-independent pKa thermochemistry and ensemble free energies."""

import logging
import math
import os
import re
from pathlib import Path

from chemsmart.analysis.thermochemistry import Thermochemistry
from chemsmart.utils.constants import (
    HARTREE_TO_KCAL_MOL,
    R,
    atm_to_pa,
    energy_conversion,
)
from chemsmart.utils.repattern import conformer_index_suffix_pattern

logger = logging.getLogger(__name__)

_THERMOCHEMISTRY_TYPE = Thermochemistry
_PKA_CONFORMER_SUFFIX_RE = re.compile(conformer_index_suffix_pattern)


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


def _record_parsed_temperature(thermo, filepath, records):
    """Store a parsed thermochemistry temperature when the parser has one."""
    if records is None or type(thermo) is not _THERMOCHEMISTRY_TYPE:
        return
    temperature = thermo.file_object.temperature_in_K
    if temperature is None:
        return
    records.append((str(filepath), float(temperature)))


def _warn_inconsistent_conformer_temperatures(label, records):
    """Warn when one species' parsed conformer temperatures differ."""
    if len(records) < 2:
        return
    unique = {round(temperature, 4) for _, temperature in records}
    if len(unique) < 2:
        return
    details = ", ".join(
        f"{Path(path).name}={temperature:g} K" for path, temperature in records
    )
    logger.warning(
        "%s conformer outputs use inconsistent temperatures (%s).",
        label,
        details,
    )


def _pka_conformer_tag(filepath, index):
    stem = Path(str(filepath)).stem
    match = _PKA_CONFORMER_SUFFIX_RE.search(stem)
    if match:
        return f"c{int(match.group(1))}"
    return f"conformer {index}"


def _pka_output_error(label, filepath, index, exc):
    tag = _pka_conformer_tag(filepath, index)
    return f"{label} {tag}: {exc}"


def pka_gas_phase_data(
    filepath,
    temperature=298.15,
    concentration=1.0,
    pressure=1.0,
    cutoff_entropy_grimme=100.0,
    cutoff_enthalpy=100.0,
    entropy_method="grimme",
    temperature_records=None,
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
    _record_parsed_temperature(thermo, filepath, temperature_records)
    return electronic_energy_au, qh_gibbs_au - electronic_energy_au


def pka_solvent_scf_energy(filepath, temperature_records=None):
    """Return solvent-phase SCF energy in Hartree."""
    thermo = Thermochemistry(filename=filepath)
    electronic_energy_j_mol = _extract_thermochemistry_property(
        thermo,
        filepath,
        "electronic_energy",
        "SCF energy",
    )
    _record_parsed_temperature(thermo, filepath, temperature_records)
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
    temperature_records = []
    for index, (gas_file, solv_file) in enumerate(
        zip(gas_files, solv_files), start=1
    ):
        try:
            e_gas, g_corr = pka_gas_phase_data(
                gas_file,
                temperature_records=temperature_records,
                **thermo_kwargs,
            )
        except ValueError as exc:
            raise ValueError(
                _pka_output_error(label, gas_file, index, exc)
            ) from exc
        try:
            e_solv = pka_solvent_scf_energy(
                solv_file, temperature_records=temperature_records
            )
        except ValueError as exc:
            raise ValueError(
                _pka_output_error(label, solv_file, index, exc)
            ) from exc
        e_gas_values.append(e_gas)
        g_corr_values.append(g_corr)
        e_solv_values.append(e_solv)
        g_soln_values.append(e_solv + g_corr)
    _warn_inconsistent_conformer_temperatures(label, temperature_records)
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
