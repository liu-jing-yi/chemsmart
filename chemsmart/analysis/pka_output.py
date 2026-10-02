"""Shared thermochemistry accessors for Gaussian and ORCA pKa outputs."""

import os

from chemsmart.utils.constants import energy_conversion


class PKaOutputThermochemistryMixin:
    """Thermochemistry properties and pKa pass-throughs for one output file.

    Subclasses call :meth:`_init_pka_output_thermochemistry` after their
    file-parser ``__init__``. ``Thermochemistry`` is imported lazily so
    Gaussian and ORCA output modules can inherit this mixin without a
    circular import.
    """

    def _init_pka_output_thermochemistry(
        self,
        temperature=298.15,
        concentration=1.0,
        pressure=1.0,
        cutoff_entropy_grimme=100.0,
        cutoff_enthalpy=100.0,
        entropy_method="grimme",
        energy_units="hartree",
    ):
        self.temperature = temperature
        self.concentration = concentration
        self.pressure = pressure
        self.cutoff_entropy_grimme = cutoff_entropy_grimme
        self.cutoff_enthalpy = cutoff_enthalpy
        self.entropy_method = entropy_method
        self.energy_units = energy_units.lower()
        self._thermochemistry = None

    @property
    def thermochemistry(self):
        """Configured ``Thermochemistry`` object for this output file."""
        if self._thermochemistry is None:
            from chemsmart.analysis.thermochemistry import Thermochemistry

            self._thermochemistry = Thermochemistry(
                filename=self.filename,
                temperature=self.temperature,
                concentration=self.concentration,
                pressure=self.pressure,
                use_weighted_mass=False,
                alpha=4,
                s_freq_cutoff=self.cutoff_entropy_grimme,
                entropy_method=self.entropy_method,
                h_freq_cutoff=self.cutoff_enthalpy,
                energy_units=self.energy_units,
                check_imaginary_frequencies=True,
            )
        return self._thermochemistry

    def _energy_in_units(self, value, label):
        if value is None:
            raise ValueError(
                f"Cannot compute {label} for {self.filename}. "
                "The file may not contain frequency calculation data."
            )
        return energy_conversion("j/mol", self.energy_units, value)

    @property
    def electronic_energy_in_units(self):
        """Electronic energy (E) in the configured units."""
        return energy_conversion(
            "j/mol",
            self.energy_units,
            self.thermochemistry.electronic_energy,
        )

    @property
    def qh_gibbs_free_energy(self):
        """Quasi-harmonic Gibbs free energy qh-G(T) in the configured units."""
        return self._energy_in_units(
            self.thermochemistry.qrrho_gibbs_free_energy,
            "qh-Gibbs free energy",
        )

    @property
    def zero_point_energy_in_units(self):
        """Zero-point energy in the configured units."""
        return self._energy_in_units(
            self.thermochemistry.zero_point_energy, "zero-point energy"
        )

    @property
    def enthalpy_in_units(self):
        """Enthalpy (H) in the configured units."""
        return self._energy_in_units(self.thermochemistry.enthalpy, "enthalpy")

    @property
    def qh_enthalpy_in_units(self):
        """Quasi-harmonic enthalpy qh-H(T) in the configured units."""
        return self._energy_in_units(
            self.thermochemistry.qrrho_enthalpy, "qh-enthalpy"
        )

    @property
    def gibbs_free_energy_in_units(self):
        """Uncorrected Gibbs free energy G(T) in the configured units."""
        return self._energy_in_units(
            self.thermochemistry.gibbs_free_energy, "Gibbs free energy"
        )

    @property
    def thermochemical_properties(self):
        """Return the standard set of thermochemical energies."""
        return {
            "electronic_energy": self.electronic_energy_in_units,
            "zero_point_energy": self.zero_point_energy_in_units,
            "enthalpy": self.enthalpy_in_units,
            "qh_enthalpy": self.qh_enthalpy_in_units,
            "gibbs_free_energy": self.gibbs_free_energy_in_units,
            "qh_gibbs_free_energy": self.qh_gibbs_free_energy,
        }

    def compute_thermochemistry(self):
        """Return thermochemistry values for this output file."""
        thermo = self.thermochemistry
        return {
            "structure": os.path.splitext(os.path.basename(self.filename))[0],
            "electronic_energy": self.electronic_energy_in_units,
            "zero_point_energy": self.zero_point_energy_in_units,
            "enthalpy": self.enthalpy_in_units,
            "qh_enthalpy": self.qh_enthalpy_in_units,
            "entropy_times_temperature": (
                energy_conversion(
                    "j/mol",
                    self.energy_units,
                    thermo.entropy_times_temperature,
                )
                if thermo.entropy_times_temperature
                else None
            ),
            "qh_entropy_times_temperature": (
                energy_conversion(
                    "j/mol",
                    self.energy_units,
                    thermo.qrrho_entropy_times_temperature,
                )
                if thermo.qrrho_entropy_times_temperature
                else None
            ),
            "gibbs_free_energy": self.gibbs_free_energy_in_units,
            "qh_gibbs_free_energy": self.qh_gibbs_free_energy,
        }

    @staticmethod
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
    ):
        """Compute thermochemistry for pKa species (HA, A-, HRef, Ref-)."""
        from chemsmart.analysis.pka import compute_pka_thermochemistry

        return compute_pka_thermochemistry(
            ha_file=ha_file,
            a_file=a_file,
            href_file=href_file,
            ref_file=ref_file,
            temperature=temperature,
            concentration=concentration,
            pressure=pressure,
            cutoff_entropy_grimme=cutoff_entropy_grimme,
            cutoff_enthalpy=cutoff_enthalpy,
            energy_units=energy_units,
        )

    @staticmethod
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
        """Compute pKa using a dual-level thermodynamic cycle."""
        from chemsmart.analysis.pka import compute_pka

        return compute_pka(
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

    @staticmethod
    def print_pka_summary(
        ha_gas_file,
        a_gas_file,
        href_gas_file,
        ref_gas_file,
        ha_solv_file,
        a_solv_file,
        href_solv_file,
        ref_solv_file,
        pka_reference,
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
        from chemsmart.analysis.pka import print_pka_summary

        return print_pka_summary(
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
            pkb=pkb,
            pks=pks,
            solvent_id=solvent_id,
        )
