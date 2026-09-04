"""Shared pKa thermodynamic-cycle chain for Gaussian and ORCA jobs.

``PkaChainMixin`` builds Opt (HA, A−) → Ref Opt → SP → Ref SP.
Program jobs set ``_opt_job_class`` and ``_sp_job_class``.
"""

import os

from chemsmart.jobs.chain import ChainMixin, JobPhase


class PkaChainMixin(ChainMixin):
    """HA/A− and optional reference Opt → SP workflow on top of ``ChainMixin``.

    Subclasses set ``_opt_job_class`` and ``_sp_job_class``.
    """

    _opt_job_class = None
    _sp_job_class = None
    _shared_reference_molecule_cache = {}

    @classmethod
    def _reference_cache_key(cls, settings):
        if settings is None or not settings.has_reference_file:
            return None
        return (
            settings.scheme,
            os.path.abspath(settings.reference_file),
            settings.reference_proton_index,
            settings.reference_charge,
            settings.reference_multiplicity,
            settings.reference_conjugate_base_charge,
            settings.reference_conjugate_base_multiplicity,
        )

    @classmethod
    def _get_cached_reference_pair(cls, settings):
        cache_key = cls._reference_cache_key(settings)
        if cache_key is None:
            return None
        if cache_key not in cls._shared_reference_molecule_cache:
            cls._shared_reference_molecule_cache[cache_key] = (
                settings.reference_pair_molecules()
            )
        return cls._shared_reference_molecule_cache[cache_key]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.opt_jobs = []
        self.ref_opt_jobs = []
        self.sp_jobs = None
        self.ref_sp_jobs = None

        self.protonated_job = None
        self.conjugate_base_job = None
        self.protonated_sp_job = None
        self.conjugate_base_sp_job = None

        self.ref_acid_job = None
        self.ref_conjugate_base_job = None
        self.ref_acid_sp_job = None
        self.ref_conjugate_base_sp_job = None

        self._prepare_pka_jobs()
        self.phases = self._build_phases()

    @property
    def has_reference_jobs(self):
        return bool(self.settings and self.settings.has_reference_file)

    def _ha_opt_label(self):
        return f"{self.label}_HA_opt"

    def _a_opt_label(self):
        return f"{self.label}_A_opt"

    def _href_opt_label(self):
        return f"{self.label}_HRef_opt"

    def _ref_opt_label(self):
        return f"{self.label}_Ref_opt"

    def _ha_sp_label(self):
        return f"{self.label}_HA_sp"

    def _a_sp_label(self):
        return f"{self.label}_A_sp"

    def _href_sp_label(self):
        return f"{self.label}_HRef_sp"

    def _ref_sp_label(self):
        return f"{self.label}_Ref_sp"

    def _finalize_child_job(self, job, legacy_label=None):
        return job

    def _child_job(
        self, job_class, molecule, settings, label, legacy_label=None
    ):
        job = job_class(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self._finalize_child_job(job, legacy_label=legacy_label)
        return job

    def _prepare_pka_jobs(self):
        if self.settings is None:
            return

        prot_mol, conj_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        prot_opt_settings, conj_opt_settings = (
            self.settings.conjugate_pair_job_settings(self.molecule)
        )

        self.protonated_job = self._child_job(
            self._opt_job_class,
            prot_mol,
            prot_opt_settings,
            self._ha_opt_label(),
            legacy_label=self.label,
        )
        self.conjugate_base_job = self._child_job(
            self._opt_job_class,
            conj_mol,
            conj_opt_settings,
            self._a_opt_label(),
            legacy_label=f"{self.label}_cb",
        )
        self.opt_jobs = [self.protonated_job, self.conjugate_base_job]

        if self.has_reference_jobs:
            self.ref_opt_jobs = self._prepare_ref_opt_jobs()
            self.ref_acid_job, self.ref_conjugate_base_job = self.ref_opt_jobs

    def _prepare_ref_opt_jobs(self):
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            ref_acid_mol, ref_cb_mol = self.settings.reference_pair_molecules()
        else:
            ref_acid_mol, ref_cb_mol = reference_pair
        ref_acid_settings, ref_cb_settings = (
            self.settings.reference_pair_job_settings()
        )

        href_label = self._href_opt_label()
        ref_label = self._ref_opt_label()
        ref_acid_job = self._child_job(
            self._opt_job_class,
            ref_acid_mol,
            ref_acid_settings,
            href_label,
            legacy_label=href_label,
        )
        ref_cb_job = self._child_job(
            self._opt_job_class,
            ref_cb_mol,
            ref_cb_settings,
            ref_label,
            legacy_label=ref_label,
        )
        return [ref_acid_job, ref_cb_job]

    def _build_phases(self):
        return [
            JobPhase(
                name="Opt",
                jobs_factory=lambda: self.opt_jobs,
                stop_on_incomplete=True,
                require_complete=True,
                stop_message="Opt jobs incomplete, halting serial execution.",
            ),
            JobPhase(
                name="Ref Opt",
                jobs_factory=lambda: self.ref_opt_jobs,
                skip_if=lambda: not self.has_reference_jobs,
                stop_on_incomplete=True,
                require_complete=True,
                stop_message=(
                    "Ref Opt jobs incomplete, halting serial execution."
                ),
            ),
            JobPhase(
                name="SP",
                jobs_factory=self._sp_jobs_for_phase,
                stop_on_incomplete=True,
                require_complete=True,
                stop_message="SP jobs incomplete, halting serial execution.",
            ),
            JobPhase(
                name="Ref SP",
                jobs_factory=self._ref_sp_jobs_for_phase,
                skip_if=lambda: not self.has_reference_jobs,
                stop_on_incomplete=True,
                require_complete=True,
                stop_message=(
                    "Ref SP jobs incomplete, halting serial execution."
                ),
            ),
        ]

    def _sp_jobs_for_phase(self):
        if self.sp_jobs is None:
            self._create_sp_jobs()
        return self.sp_jobs

    def _ref_sp_jobs_for_phase(self):
        if not self.has_reference_jobs:
            return []
        if self.ref_sp_jobs is None:
            self._create_ref_sp_jobs()
        return self.ref_sp_jobs

    def _create_sp_jobs(self):
        prot_opt_mol = self._output_molecule(
            self.protonated_job, self.protonated_job.molecule
        )
        conj_opt_mol = self._output_molecule(
            self.conjugate_base_job, self.conjugate_base_job.molecule
        )
        prot_sp_settings, conj_sp_settings = (
            self.settings.conjugate_pair_sp_job_settings(self.molecule)
        )

        self.protonated_sp_job = self._child_job(
            self._sp_job_class,
            prot_opt_mol,
            prot_sp_settings,
            self._ha_sp_label(),
            legacy_label=f"{self.label}_sp",
        )
        self.conjugate_base_sp_job = self._child_job(
            self._sp_job_class,
            conj_opt_mol,
            conj_sp_settings,
            self._a_sp_label(),
            legacy_label=f"{self.label}_cb_sp",
        )
        self.sp_jobs = [self.protonated_sp_job, self.conjugate_base_sp_job]

    def _create_ref_sp_jobs(self):
        ref_acid_opt_mol = self._output_molecule(
            self.ref_acid_job, self.ref_acid_job.molecule
        )
        ref_cb_opt_mol = self._output_molecule(
            self.ref_conjugate_base_job, self.ref_conjugate_base_job.molecule
        )
        ref_acid_sp_settings, ref_cb_sp_settings = (
            self.settings.reference_pair_sp_job_settings()
        )

        href_sp_label = self._href_sp_label()
        ref_sp_label = self._ref_sp_label()
        self.ref_acid_sp_job = self._child_job(
            self._sp_job_class,
            ref_acid_opt_mol,
            ref_acid_sp_settings,
            href_sp_label,
            legacy_label=href_sp_label,
        )
        self.ref_conjugate_base_sp_job = self._child_job(
            self._sp_job_class,
            ref_cb_opt_mol,
            ref_cb_sp_settings,
            ref_sp_label,
            legacy_label=ref_sp_label,
        )
        self.ref_sp_jobs = [
            self.ref_acid_sp_job,
            self.ref_conjugate_base_sp_job,
        ]

    def _opt_jobs_are_complete(self):
        if not self.opt_jobs:
            return False
        return all(job.is_complete() for job in self.opt_jobs)

    def _ref_opt_jobs_are_complete(self):
        if not self.has_reference_jobs:
            return True
        if not self.ref_opt_jobs:
            return False
        return all(job.is_complete() for job in self.ref_opt_jobs)

    @property
    def protonated_output(self):
        if self.protonated_job is None:
            return None
        return self.protonated_job._output()

    @property
    def conjugate_base_output(self):
        if self.conjugate_base_job is None:
            return None
        return self.conjugate_base_job._output()

    @property
    def protonated_sp_output(self):
        if self.protonated_sp_job is None:
            return None
        return self.protonated_sp_job._output()

    @property
    def conjugate_base_sp_output(self):
        if self.conjugate_base_sp_job is None:
            return None
        return self.conjugate_base_sp_job._output()

    @property
    def ref_acid_output(self):
        if not self.has_reference_jobs or self.ref_acid_job is None:
            return None
        return self.ref_acid_job._output()

    @property
    def ref_conjugate_base_output(self):
        if not self.has_reference_jobs or self.ref_conjugate_base_job is None:
            return None
        return self.ref_conjugate_base_job._output()

    @property
    def ref_acid_sp_output(self):
        if not self.has_reference_jobs or self.ref_acid_sp_job is None:
            return None
        return self.ref_acid_sp_job._output()

    @property
    def ref_conjugate_base_sp_output(self):
        if (
            not self.has_reference_jobs
            or self.ref_conjugate_base_sp_job is None
        ):
            return None
        return self.ref_conjugate_base_sp_job._output()

    def _pka_output_files(self):
        ha_file = (
            self.protonated_job.outputfile if self.protonated_job else None
        )
        a_file = (
            self.conjugate_base_job.outputfile
            if self.conjugate_base_job
            else None
        )
        href_file = None
        ref_file = None
        if self.has_reference_jobs and self._ref_opt_jobs_are_complete():
            href_file = (
                self.ref_acid_job.outputfile if self.ref_acid_job else None
            )
            ref_file = (
                self.ref_conjugate_base_job.outputfile
                if self.ref_conjugate_base_job
                else None
            )
        return ha_file, a_file, href_file, ref_file

    def _reference_pka(self):
        return None

    def _thermo_kwargs(self):
        return dict(
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
        )

    def compute_thermochemistry(self):
        """Compute and return thermochemistry results for all species."""
        from chemsmart.analysis.pka import compute_pka_thermochemistry

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot compute thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        ha_file, a_file, href_file, ref_file = self._pka_output_files()
        return compute_pka_thermochemistry(
            ha_file=ha_file,
            a_file=a_file,
            href_file=href_file,
            ref_file=ref_file,
            energy_units=self.settings.energy_units,
            **self._thermo_kwargs(),
        )

    def print_thermochemistry(self):
        """Print formatted thermochemistry summary to stdout."""
        from chemsmart.analysis.pka import print_pka_summary

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot print thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        ha_gas, a_gas, href_gas, ref_gas = self._pka_output_files()
        ha_solv = (
            self.protonated_sp_job.outputfile
            if self.protonated_sp_job
            else None
        )
        a_solv = (
            self.conjugate_base_sp_job.outputfile
            if self.conjugate_base_sp_job
            else None
        )
        href_solv = ref_solv = None
        if self.has_reference_jobs:
            href_solv = (
                self.ref_acid_sp_job.outputfile
                if self.ref_acid_sp_job
                else None
            )
            ref_solv = (
                self.ref_conjugate_base_sp_job.outputfile
                if self.ref_conjugate_base_sp_job
                else None
            )

        print_pka_summary(
            ha_gas_file=ha_gas,
            a_gas_file=a_gas,
            href_gas_file=href_gas,
            ref_gas_file=ref_gas,
            ha_solv_file=ha_solv,
            a_solv_file=a_solv,
            href_solv_file=href_solv,
            ref_solv_file=ref_solv,
            pka_reference=self._reference_pka(),
            scheme=self.settings.scheme,
            delta_G_proton=self.settings.delta_G_proton,
            **self._thermo_kwargs(),
        )
