"""Shared pKa molecule and reference-acid settings.

Gaussian and ORCA pKa settings use this mixin for proton removal,
reference-acid loading, and charge/multiplicity defaults. Gas-phase
and solvent job settings stay on the program-specific classes.
"""

_PKA_SHARED_RENAME = {
    "reference": "reference_file",
    "delta_g_proton": "delta_G_proton",
}
_PKA_SHARED_IGNORE = {
    "entropy_method",
    "preview",
    "reference_color_code",
    "skip_completed",
}
_PKA_SHARED_FIELDS = {
    "concentration",
    "conjugate_base_charge",
    "conjugate_base_multiplicity",
    "crest_project",
    "cutoff_enthalpy",
    "cutoff_entropy_grimme",
    "delta_g_proton",
    "energy_units",
    "num_conformers",
    "pkb",
    "pks",
    "pressure",
    "reference",
    "reference_charge",
    "reference_conjugate_base_charge",
    "reference_conjugate_base_multiplicity",
    "reference_multiplicity",
    "reference_proton_index",
    "sampling",
    "scheme",
    "solvent_id",
    "solvent_model",
    "temperature",
}


def pka_kwargs_from_shared(
    shared,
    *,
    opt_settings,
    default_solvent_model,
    default_solvent_id,
    sp_settings=None,
):
    """Map a CLI shared dict onto pKa settings constructor kwargs.

    Known CLI-only keys are dropped. Unknown keys raise ``ValueError``.
    Solvent model and ID fall back to *opt_settings*, then *sp_settings*,
    then the program defaults.
    """
    unknown = set(shared) - _PKA_SHARED_FIELDS - _PKA_SHARED_IGNORE
    if unknown:
        raise ValueError(
            "Unknown pKa shared options: " + ", ".join(sorted(unknown))
        )

    kwargs = {}
    for key, value in shared.items():
        if key in _PKA_SHARED_IGNORE:
            continue
        kwargs[_PKA_SHARED_RENAME.get(key, key)] = value

    solvent_model = kwargs.get("solvent_model")
    if solvent_model is None:
        solvent_model = opt_settings.solvent_model
    if solvent_model is None and sp_settings is not None:
        solvent_model = sp_settings.solvent_model
    if solvent_model is None:
        solvent_model = default_solvent_model

    solvent_id = kwargs.get("solvent_id")
    if solvent_id is None:
        solvent_id = opt_settings.solvent_id
    if solvent_id is None and sp_settings is not None:
        solvent_id = sp_settings.solvent_id
    if solvent_id is None:
        solvent_id = default_solvent_id

    kwargs["solvent_model"] = solvent_model
    kwargs["solvent_id"] = solvent_id
    return kwargs


class PKaMoleculeSettingsMixin:
    """Proton removal, reference loading, and charge defaults for pKa."""

    def _store_pka_common_settings(
        self,
        *,
        proton_index,
        scheme,
        reference_file,
        reference_proton_index,
        reference_charge,
        reference_multiplicity,
        reference_conjugate_base_charge,
        reference_conjugate_base_multiplicity,
        reference_pka,
        delta_G_proton,
        solvent_model,
        solvent_id,
        conjugate_base_charge,
        conjugate_base_multiplicity,
        temperature,
        concentration,
        pressure,
        cutoff_entropy_grimme,
        cutoff_enthalpy,
        energy_units,
        pkb,
        pks,
        sampling,
        num_conformers,
        crest_project,
        default_title,
    ):
        """Store the pKa fields shared by Gaussian and ORCA settings."""
        self.proton_index = proton_index
        self.scheme = scheme
        self.solvent_model = solvent_model
        self.solvent_id = solvent_id
        self.conjugate_base_charge = conjugate_base_charge
        self.conjugate_base_multiplicity = conjugate_base_multiplicity
        self.temperature = temperature
        self.concentration = concentration
        self.pressure = pressure
        self.cutoff_entropy_grimme = cutoff_entropy_grimme
        self.cutoff_enthalpy = cutoff_enthalpy
        self.energy_units = energy_units
        self.reference_pka = reference_pka
        self.delta_G_proton = delta_G_proton
        self.pkb = bool(pkb)
        self.pks = pks
        self.sampling = bool(sampling)
        if num_conformers is None:
            num_conformers = 1
        if num_conformers < 1:
            raise ValueError("num_conformers must be >= 1.")
        self.num_conformers = int(num_conformers)
        self.crest_project = crest_project

        if not self.title:
            self.title = default_title

        if scheme == "proton exchange":
            self.reference_file = reference_file
            self.reference_proton_index = reference_proton_index
            self.reference_charge = reference_charge
            self.reference_multiplicity = reference_multiplicity
            self.reference_conjugate_base_charge = (
                reference_conjugate_base_charge
            )
            self.reference_conjugate_base_multiplicity = (
                reference_conjugate_base_multiplicity
            )
        else:
            self.reference_file = None
            self.reference_proton_index = None
            self.reference_charge = None
            self.reference_multiplicity = None
            self.reference_conjugate_base_charge = None
            self.reference_conjugate_base_multiplicity = None

        from chemsmart.analysis.pka import (
            resolve_pkb_reporting,
            warn_if_default_pks_non_aqueous,
            warn_if_non_aqueous_direct_proton_default,
        )

        warn_if_non_aqueous_direct_proton_default(
            self.scheme, self.delta_G_proton, self.solvent_id
        )
        _, _, pks_defaulted = resolve_pkb_reporting(pkb=self.pkb, pks=self.pks)
        warn_if_default_pks_non_aqueous(pks_defaulted, self.solvent_id)

    @property
    def has_reference_file(self):
        """Return True when a proton-exchange reference geometry is set."""
        return (
            self.scheme == "proton exchange"
            and self.reference_file is not None
        )

    def validate_reference_settings(self):
        """Require proton index, charge, and multiplicity for a reference file.

        Returns without error when no reference file is configured.
        """
        if not self.has_reference_file:
            return

        missing = []
        if self.reference_proton_index is None:
            missing.append("reference_proton_index")
        if self.reference_charge is None:
            missing.append("reference_charge")
        if self.reference_multiplicity is None:
            missing.append("reference_multiplicity")
        if missing:
            raise ValueError(
                "When reference_file is provided, the following must also "
                f"be specified: {', '.join(missing)}"
            )

    def get_reference_molecule(self):
        """Load the reference acid and set its charge and multiplicity."""
        if not self.has_reference_file:
            raise ValueError(
                "Reference file not provided. Cannot load reference molecule."
            )
        self.validate_reference_settings()
        from chemsmart.io.molecules.structure import Molecule

        ref_mol = Molecule.from_filepath(self.reference_file)
        ref_mol.charge = self.reference_charge
        ref_mol.multiplicity = self.reference_multiplicity
        return ref_mol

    def get_reference_conjugate_base_molecule(self):
        """Return the reference acid with its acidic proton removed."""
        return self._create_reference_conjugate_base_molecule(
            self.get_reference_molecule()
        )

    def reference_pair_molecules(self):
        """Return the reference acid and its conjugate base."""
        ref_mol = self.get_reference_molecule()
        return (
            ref_mol,
            self._create_reference_conjugate_base_molecule(ref_mol),
        )

    def _create_reference_conjugate_base_molecule(self, reference_molecule):
        ref_cb = self._remove_acidic_hydrogen(
            reference_molecule,
            self.reference_proton_index,
            "reference_proton_index",
        )
        charge, multiplicity = (
            self._reference_conjugate_base_charge_multiplicity()
        )
        ref_cb.charge = charge
        ref_cb.multiplicity = multiplicity
        return ref_cb

    @property
    def protonated_charge(self):
        """Charge of the protonated form."""
        return self.charge

    @protonated_charge.setter
    def protonated_charge(self, value):
        self.charge = value

    @property
    def protonated_multiplicity(self):
        """Multiplicity of the protonated form."""
        return self.multiplicity

    @protonated_multiplicity.setter
    def protonated_multiplicity(self, value):
        self.multiplicity = value

    def protonated_molecule(self, molecule):
        """Return a copy of HA with the protonated charge and multiplicity."""
        protonated = molecule.copy()
        charge, multiplicity = self._protonated_charge_multiplicity(molecule)
        protonated.charge = charge
        protonated.multiplicity = multiplicity
        return protonated

    def conjugate_base_molecule(self, molecule):
        """Return A- by removing the acidic proton from HA."""
        return self._create_conjugate_base_molecule(molecule)

    def conjugate_pair_molecules(self, molecule):
        """Return HA and A-."""
        return molecule, self._create_conjugate_base_molecule(molecule)

    def _create_molecules(self, molecule):
        """Return a charge-updated HA copy and A-."""
        return (
            self.protonated_molecule(molecule),
            self._create_conjugate_base_molecule(molecule),
        )

    def _create_conjugate_base_molecule(self, molecule):
        conjugate = self._remove_acidic_hydrogen(
            molecule, self.proton_index, "proton_index"
        )
        original_charge = 0 if molecule.charge is None else molecule.charge
        original_mult = (
            1 if molecule.multiplicity is None else molecule.multiplicity
        )
        if self.conjugate_base_charge is not None:
            conjugate.charge = self.conjugate_base_charge
        else:
            conjugate.charge = original_charge - 1
        if self.conjugate_base_multiplicity is not None:
            conjugate.multiplicity = self.conjugate_base_multiplicity
        else:
            conjugate.multiplicity = original_mult
        return conjugate

    def _protonated_charge_multiplicity(self, molecule):
        """Return HA charge and multiplicity from settings, then the molecule."""
        if self.charge is not None:
            charge = self.charge
        elif molecule.charge is not None:
            charge = molecule.charge
        else:
            charge = 0
        if self.multiplicity is not None:
            multiplicity = self.multiplicity
        elif molecule.multiplicity is not None:
            multiplicity = molecule.multiplicity
        else:
            multiplicity = 1
        return charge, multiplicity

    def _conjugate_base_charge_multiplicity(self, prot_charge, prot_mult):
        """Return A- charge and multiplicity for sub-job settings."""
        if self.conjugate_base_charge is not None:
            charge = self.conjugate_base_charge
        else:
            charge = prot_charge - 1
        if self.conjugate_base_multiplicity is not None:
            multiplicity = self.conjugate_base_multiplicity
        else:
            multiplicity = prot_mult
        return charge, multiplicity

    def _reference_conjugate_base_charge_multiplicity(self):
        """Return Ref- charge and multiplicity."""
        if self.reference_conjugate_base_charge is not None:
            charge = self.reference_conjugate_base_charge
        else:
            charge = self.reference_charge - 1
        if self.reference_conjugate_base_multiplicity is not None:
            multiplicity = self.reference_conjugate_base_multiplicity
        else:
            multiplicity = self.reference_multiplicity
        return charge, multiplicity

    def _remove_acidic_hydrogen(self, molecule, index, index_name):
        """Remove one hydrogen with ``Molecule.delete_atoms_by_indices``."""
        if index is None:
            raise ValueError(
                f"{index_name} must be specified to create the conjugate "
                "base molecule. Use 1-based indexing."
            )
        if index < 1 or index > len(molecule):
            raise ValueError(
                f"{index_name} {index} is out of range. "
                f"Molecule has {len(molecule)} atoms "
                f"(1-indexed: 1 to {len(molecule)})."
            )
        symbol = molecule.symbols[index - 1]
        if symbol not in ("H", "h"):
            raise ValueError(
                f"Atom at {index_name} {index} is '{symbol}', not hydrogen."
            )
        return molecule.delete_atoms_by_indices(index, one_based=True)
