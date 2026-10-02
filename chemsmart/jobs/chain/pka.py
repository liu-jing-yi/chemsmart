"""Shared pKa phase machine for program-specific job subclasses.

Subclasses supply ``opt_job_class`` and ``sp_job_class``. Optional hooks
are ``_bind_subjob``,
``_sync_subjob_folders``, ``_subjob_label``, ``_subjob_legacy_label``,
and ``_crest_job_label``.
"""

import logging
import os

from chemsmart.jobs.chain.sampling import (
    build_pka_crest_job,
    select_crest_conformers,
)
from chemsmart.jobs.runner import decide_phase_transition, run_phase_jobs
from chemsmart.utils.datasets import pka_job_species_outputs, pka_subjob_label

logger = logging.getLogger(__name__)


class PKaJob:
    """Gas-phase opt, solvent single-point, and optional CREST phases.

    Program jobs mix this in ahead of their ``Job`` base so these methods
    replace that base's ``_run``. Each subclass keeps its own reference
    molecule cache.
    """

    opt_job_class = None
    sp_job_class = None

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        cls._shared_reference_molecule_cache = {}

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
        cache = cls._shared_reference_molecule_cache
        if cache_key not in cache:
            cache[cache_key] = settings.reference_pair_molecules()
        return cache[cache_key]

    def _reset_pka_child_jobs(self):
        """Clear child-job state and build the initial opt/SP phase."""
        self.crest_jobs = []
        self.protonated_crest_job = None
        self.conjugate_base_crest_job = None
        self.ref_acid_crest_job = None
        self.ref_conjugate_base_crest_job = None

        self.ref_opt_jobs = []
        self.ref_acid_opt_jobs = []
        self.ref_conjugate_base_opt_jobs = []
        self.ref_sp_jobs = None
        self.ref_acid_sp_jobs = None
        self.ref_conjugate_base_sp_jobs = None
        self.ref_acid_job = None
        self.ref_conjugate_base_job = None
        self.ref_acid_sp_job = None
        self.ref_conjugate_base_sp_job = None

        self.has_reference_jobs = bool(self.settings.has_reference_file)
        self._prepare_pka_jobs()

    def _bind_subjob(self, job, legacy_label=None):
        """Attach program-specific output lookup to a child job."""

    def _sync_subjob_folders(self, jobs):
        """Point child jobs at the parent folder before a phase runs."""

    def _subjob_legacy_label(self, species, stage):
        """Return a pre-rename output stem for one conformer, or None.

        ``species`` is ``HA``, ``A``, ``HRef``, or ``Ref``. ``stage`` is
        ``opt`` or ``sp``. The stem is used only when the ensemble has one
        conformer.
        """
        return None

    def _subjob_label(self, species, stage, index, num_conformers):
        """Return the DFT label for one species, stage, and conformer."""
        return pka_subjob_label(
            self.label, species, stage, index, num_conformers
        )

    def _crest_job_label(self, species):
        """Return the CREST label for ``HA``, ``A``, ``HRef``, or ``Ref``."""
        return f"{self.label}_{species}_crest"

    def _optimized_molecule_from_job(
        self, job, fallback_molecule, legacy_label=None
    ):
        """Return a finished optimized geometry, or the fallback molecule."""
        out = job._output()
        if out is not None and out.normal_termination:
            return out.molecule
        return fallback_molecule

    def _reference_pair_molecules(self):
        """Return the HRef and Ref- molecules for the configured reference acid."""
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            return self.settings.reference_pair_molecules()
        return reference_pair

    def _prepare_pka_jobs(self):
        """Prepare CREST, gas-phase opt, and solvent SP jobs."""
        prot_mol, conj_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        if self.settings.sampling:
            self._attach_crest_job("HA", prot_mol)
            self._attach_crest_job("A", conj_mol)
            n = self.settings.num_conformers
            ha_molecules = [prot_mol] * n
            a_molecules = [conj_mol] * n
        else:
            ha_molecules = [prot_mol]
            a_molecules = [conj_mol]

        self._prepare_target_opt_jobs(ha_molecules, a_molecules)
        self._create_sp_jobs()

        if not self.has_reference_jobs:
            return

        href_mol, ref_mol = self._reference_pair_molecules()
        if self.settings.sampling:
            self._attach_crest_job("HRef", href_mol)
            self._attach_crest_job("Ref", ref_mol)
            n = self.settings.num_conformers
            href_molecules = [href_mol] * n
            ref_molecules = [ref_mol] * n
        else:
            href_molecules = [href_mol]
            ref_molecules = [ref_mol]
        self._prepare_ref_opt_jobs(href_molecules, ref_molecules)
        self._create_ref_sp_jobs()

    def _attach_crest_job(self, species, molecule):
        job = build_pka_crest_job(
            molecule,
            self._crest_job_label(species),
            self.settings,
            self,
        )
        if species == "HA":
            self.protonated_crest_job = job
        elif species == "A":
            self.conjugate_base_crest_job = job
        elif species == "HRef":
            self.ref_acid_crest_job = job
        elif species == "Ref":
            self.ref_conjugate_base_crest_job = job
        else:
            raise ValueError(f"Unknown pKa CREST species {species!r}.")
        self.crest_jobs.append(job)
        return job

    def _prepare_target_opt_jobs(self, ha_molecules, a_molecules):
        """Build gas-phase opt+freq jobs for HA and A- conformers."""
        prot_opt_settings, conj_opt_settings = (
            self.settings.conjugate_pair_job_settings(self.molecule)
        )
        self.protonated_opt_jobs = self._make_species_opt_jobs(
            ha_molecules, prot_opt_settings, "HA"
        )
        self.conjugate_base_opt_jobs = self._make_species_opt_jobs(
            a_molecules, conj_opt_settings, "A"
        )
        self.protonated_job = self.protonated_opt_jobs[0]
        self.conjugate_base_job = self.conjugate_base_opt_jobs[0]
        self.opt_jobs = self.protonated_opt_jobs + self.conjugate_base_opt_jobs
        self._clear_target_sp_jobs()

    def _prepare_ref_opt_jobs(self, href_molecules, ref_molecules):
        """Build gas-phase opt+freq jobs for HRef and Ref- conformers."""
        ref_acid_settings, ref_conjugate_base_settings = (
            self.settings.reference_pair_job_settings()
        )
        self.ref_acid_opt_jobs = self._make_species_opt_jobs(
            href_molecules, ref_acid_settings, "HRef"
        )
        self.ref_conjugate_base_opt_jobs = self._make_species_opt_jobs(
            ref_molecules, ref_conjugate_base_settings, "Ref"
        )
        self.ref_acid_job = self.ref_acid_opt_jobs[0]
        self.ref_conjugate_base_job = self.ref_conjugate_base_opt_jobs[0]
        self.ref_opt_jobs = (
            self.ref_acid_opt_jobs + self.ref_conjugate_base_opt_jobs
        )
        self._clear_reference_sp_jobs()

    def _make_species_opt_jobs(self, molecules, settings, species):
        num_conformers = len(molecules)
        legacy_label = (
            self._subjob_legacy_label(species, "opt")
            if num_conformers == 1
            else None
        )
        jobs = []
        for index, molecule in enumerate(molecules, start=1):
            jobs.append(
                self._make_child_job(
                    self.opt_job_class,
                    molecule,
                    settings,
                    self._subjob_label(species, "opt", index, num_conformers),
                    legacy_label=legacy_label,
                )
            )
        return jobs

    def _create_sp_jobs(self):
        """Create one solution-phase SP job per optimized conformer."""
        prot_sp_settings, conj_sp_settings = (
            self.settings.conjugate_pair_sp_job_settings(self.molecule)
        )
        self.protonated_sp_jobs = self._make_species_sp_jobs(
            self.protonated_opt_jobs, prot_sp_settings, "HA"
        )
        self.conjugate_base_sp_jobs = self._make_species_sp_jobs(
            self.conjugate_base_opt_jobs, conj_sp_settings, "A"
        )
        self.protonated_sp_job = self.protonated_sp_jobs[0]
        self.conjugate_base_sp_job = self.conjugate_base_sp_jobs[0]
        self.sp_jobs = self.protonated_sp_jobs + self.conjugate_base_sp_jobs

    def _create_ref_sp_jobs(self):
        """Create one reference solvent SP job per optimized conformer."""
        ref_acid_sp_settings, ref_conjugate_base_sp_settings = (
            self.settings.reference_pair_sp_job_settings()
        )
        self.ref_acid_sp_jobs = self._make_species_sp_jobs(
            self.ref_acid_opt_jobs, ref_acid_sp_settings, "HRef"
        )
        self.ref_conjugate_base_sp_jobs = self._make_species_sp_jobs(
            self.ref_conjugate_base_opt_jobs,
            ref_conjugate_base_sp_settings,
            "Ref",
        )
        self.ref_acid_sp_job = self.ref_acid_sp_jobs[0]
        self.ref_conjugate_base_sp_job = self.ref_conjugate_base_sp_jobs[0]
        self.ref_sp_jobs = (
            self.ref_acid_sp_jobs + self.ref_conjugate_base_sp_jobs
        )

    def _make_species_sp_jobs(self, opt_jobs, settings, species):
        num_conformers = len(opt_jobs)
        legacy_label = (
            self._subjob_legacy_label(species, "sp")
            if num_conformers == 1
            else None
        )
        opt_legacy_label = (
            self._subjob_legacy_label(species, "opt")
            if num_conformers == 1
            else None
        )
        jobs = []
        for index, opt_job in enumerate(opt_jobs, start=1):
            molecule = self._optimized_molecule_from_job(
                opt_job, opt_job.molecule, legacy_label=opt_legacy_label
            )
            jobs.append(
                self._make_child_job(
                    self.sp_job_class,
                    molecule,
                    settings,
                    self._subjob_label(species, "sp", index, num_conformers),
                    legacy_label=legacy_label,
                )
            )
        return jobs

    def _make_child_job(
        self, job_class, molecule, settings, label, legacy_label=None
    ):
        if job_class is None:
            raise TypeError(
                f"{type(self).__name__} must set opt_job_class and "
                "sp_job_class."
            )
        job = job_class(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self._bind_subjob(job, legacy_label=legacy_label)
        return job

    def _clear_target_sp_jobs(self):
        self.sp_jobs = None
        self.protonated_sp_jobs = None
        self.conjugate_base_sp_jobs = None
        self.protonated_sp_job = None
        self.conjugate_base_sp_job = None

    def _clear_reference_sp_jobs(self):
        self.ref_sp_jobs = None
        self.ref_acid_sp_jobs = None
        self.ref_conjugate_base_sp_jobs = None
        self.ref_acid_sp_job = None
        self.ref_conjugate_base_sp_job = None

    def _run_opt_jobs(self):
        """Run gas phase optimization jobs."""
        self._sync_subjob_folders(self.opt_jobs)
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.opt_jobs,
            stop_on_incomplete=True,
            logger_obj=logger,
            phase_label="gas phase optimization",
        )

    def _run_ref_opt_jobs(self):
        """Run reference gas phase optimization jobs."""
        if not self.has_reference_jobs:
            return
        self._sync_subjob_folders(self.ref_opt_jobs)
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.ref_opt_jobs,
            stop_on_incomplete=True,
            logger_obj=logger,
            phase_label="reference gas phase optimization",
        )

    def _run_sp_jobs(self):
        """Rebuild and run solution-phase SP jobs from optimized geometries."""
        if not self._opt_jobs_are_complete():
            logger.warning(
                "Optimization jobs not complete. Cannot run SP jobs."
            )
            return
        self._create_sp_jobs()
        self._sync_subjob_folders(self.sp_jobs)
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.sp_jobs,
            stop_on_incomplete=True,
            logger_obj=logger,
            phase_label="solution phase SP",
        )

    def _run_ref_sp_jobs(self):
        """Rebuild and run reference solution-phase SP jobs."""
        if not self.has_reference_jobs:
            return
        if not self._ref_opt_jobs_are_complete():
            logger.warning(
                "Reference optimization jobs not complete. "
                "Cannot run reference SP jobs."
            )
            return
        self._create_ref_sp_jobs()
        self._sync_subjob_folders(self.ref_sp_jobs)
        run_phase_jobs(
            parent_runner=self.jobrunner,
            jobs=self.ref_sp_jobs,
            stop_on_incomplete=True,
            logger_obj=logger,
            phase_label="reference solution phase SP",
        )

    def _run_crest_jobs(self):
        """Run CREST sampling jobs for HA, A-, and any reference acid."""
        run_phase_jobs(
            parent_runner=None,
            jobs=self.crest_jobs,
            stop_on_incomplete=False,
            logger_obj=logger,
            phase_label="CREST sampling",
        )

    def _select_species_conformers(self, crest_job, fallback_molecule):
        return select_crest_conformers(
            crest_job, self.settings.num_conformers, fallback_molecule
        )

    def _selected_crest_conformers(self):
        """Return selected conformers, or None if any CREST job has not finished.

        The tuple is ``(HA, A-, HRef, Ref-)``. Reference entries are ``None``
        when no reference acid is configured.
        """
        prot_mol, conj_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        ha_confs = self._select_species_conformers(
            self.protonated_crest_job, prot_mol
        )
        a_confs = self._select_species_conformers(
            self.conjugate_base_crest_job, conj_mol
        )
        if ha_confs is None or a_confs is None:
            return None
        if not self.has_reference_jobs:
            return ha_confs, a_confs, None, None
        href_mol, ref_mol = self._reference_pair_molecules()
        href_confs = self._select_species_conformers(
            self.ref_acid_crest_job, href_mol
        )
        ref_confs = self._select_species_conformers(
            self.ref_conjugate_base_crest_job, ref_mol
        )
        if href_confs is None or ref_confs is None:
            return None
        return ha_confs, a_confs, href_confs, ref_confs

    def _run(self, **kwargs):
        """Run CREST when requested, then gas-phase opt and solvent SP."""
        if self.settings.sampling:
            self._run_crest_jobs()
            selected = self._selected_crest_conformers()
            crest_transition = decide_phase_transition(
                phase_name="CREST",
                require_complete=True,
                is_complete=selected is not None,
                stop_message="CREST jobs incomplete, halting serial execution.",
            )
            if not crest_transition.proceed:
                logger.info(crest_transition.message)
                return
            self._prepare_target_opt_jobs(selected[0], selected[1])
            if self.has_reference_jobs:
                self._prepare_ref_opt_jobs(selected[2], selected[3])

        self._run_opt_jobs()

        opt_transition = decide_phase_transition(
            phase_name="Opt",
            require_complete=True,
            is_complete=self._opt_jobs_are_complete(),
            stop_message="Opt jobs incomplete, halting serial execution.",
        )
        if not opt_transition.proceed:
            logger.info(opt_transition.message)
            return

        if self.has_reference_jobs:
            self._run_ref_opt_jobs()
            ref_opt_transition = decide_phase_transition(
                phase_name="Ref Opt",
                require_complete=True,
                is_complete=self._ref_opt_jobs_are_complete(),
                stop_message="Ref Opt jobs incomplete, halting serial execution.",
            )
            if not ref_opt_transition.proceed:
                logger.info(ref_opt_transition.message)
                return

        self._run_sp_jobs()

        sp_transition = decide_phase_transition(
            phase_name="SP",
            require_complete=True,
            is_complete=self._sp_jobs_are_complete(),
            stop_message="SP jobs incomplete, halting serial execution.",
        )
        if not sp_transition.proceed:
            logger.info(sp_transition.message)
            return

        if self.has_reference_jobs:
            self._run_ref_sp_jobs()

    def is_complete(self):
        """Return True when target and reference opt and SP jobs have finished.

        CREST sampling is not required for parent completion.
        """
        if not self._opt_jobs_are_complete():
            return False
        if not self._sp_jobs_are_complete():
            return False
        if self.has_reference_jobs:
            if not self._ref_opt_jobs_are_complete():
                return False
            if not self._ref_sp_jobs_are_complete():
                return False
        return True

    def _opt_jobs_are_complete(self):
        """Return True when every target gas-phase optimization job finished."""
        if not self.opt_jobs:
            return False
        return all(job.is_complete() for job in self.opt_jobs)

    def _ref_opt_jobs_are_complete(self):
        """Return True when reference opt jobs finished or are not configured."""
        if not self.has_reference_jobs:
            return True
        if not self.ref_opt_jobs:
            return False
        return all(job.is_complete() for job in self.ref_opt_jobs)

    def _sp_jobs_are_complete(self):
        """Return True when every target solution-phase SP job finished."""
        if not self.sp_jobs:
            return False
        return all(job.is_complete() for job in self.sp_jobs)

    def _ref_sp_jobs_are_complete(self):
        """Return True when reference SP jobs finished or are not configured."""
        if not self.has_reference_jobs:
            return True
        if not self.ref_sp_jobs:
            return False
        return all(job.is_complete() for job in self.ref_sp_jobs)

    @property
    def protonated_output(self):
        """Parsed output for the protonated gas-phase optimization job."""
        return self.protonated_job._output()

    @property
    def conjugate_base_output(self):
        """Parsed output for the conjugate-base gas-phase optimization job."""
        return self.conjugate_base_job._output()

    @property
    def protonated_sp_output(self):
        """Parsed output for the protonated solution-phase SP job."""
        return self.protonated_sp_job._output()

    @property
    def conjugate_base_sp_output(self):
        """Parsed output for the conjugate-base solution-phase SP job."""
        return self.conjugate_base_sp_job._output()

    @property
    def ref_acid_output(self):
        """Parsed output for the reference-acid gas-phase optimization job."""
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_job._output()

    @property
    def ref_conjugate_base_output(self):
        """Parsed output for the reference conjugate-base optimization job."""
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_job._output()

    @property
    def ref_acid_sp_output(self):
        """Parsed output for the reference-acid solution-phase SP job."""
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_sp_job._output()

    @property
    def ref_conjugate_base_sp_output(self):
        """Parsed output for the reference conjugate-base solution-phase SP job."""
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_sp_job._output()

    def _pka_output_files(self):
        """Return gas and solvent output paths for every pKa species.

        Each species maps to ``{"gas": [...], "solv": [...]}`` in conformer
        order. Reference species are included when a reference acid is set.
        """
        return pka_job_species_outputs(self)

    def _require_completed_opt_jobs(self, action):
        if self._opt_jobs_are_complete():
            return
        raise ValueError(
            f"Cannot {action}: optimization jobs are not complete. "
            "Run the pKa jobs first using job.run()."
        )

    def compute_thermochemistry(self):
        """Return thermochemistry for every pKa species in the job."""
        from chemsmart.analysis.pka import compute_pka_thermochemistry

        self._require_completed_opt_jobs("compute thermochemistry")
        files = self._pka_output_files()
        return compute_pka_thermochemistry(
            ha_file=files["HA"]["gas"],
            a_file=files["A-"]["gas"],
            href_file=files["HRef"]["gas"] if "HRef" in files else None,
            ref_file=files["Ref-"]["gas"] if "Ref-" in files else None,
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
            energy_units=self.settings.energy_units,
        )

    def print_thermochemistry(self):
        """Print the pKa summary and return the result dictionary."""
        from chemsmart.analysis.pka import print_pka_summary

        self._require_completed_opt_jobs("print thermochemistry")
        files = self._pka_output_files()
        return print_pka_summary(
            ha_gas_file=files["HA"]["gas"],
            a_gas_file=files["A-"]["gas"],
            href_gas_file=files["HRef"]["gas"] if "HRef" in files else None,
            ref_gas_file=files["Ref-"]["gas"] if "Ref-" in files else None,
            ha_solv_file=files["HA"]["solv"],
            a_solv_file=files["A-"]["solv"],
            href_solv_file=(
                files["HRef"]["solv"] if "HRef" in files else None
            ),
            ref_solv_file=files["Ref-"]["solv"] if "Ref-" in files else None,
            pka_reference=self.settings.reference_pka,
            temperature=self.settings.temperature,
            concentration=self.settings.concentration,
            pressure=self.settings.pressure,
            cutoff_entropy_grimme=self.settings.cutoff_entropy_grimme,
            cutoff_enthalpy=self.settings.cutoff_enthalpy,
            scheme=self.settings.scheme,
            delta_G_proton=self.settings.delta_G_proton,
            pkb=self.settings.pkb,
            pks=self.settings.pks,
            solvent_id=self.settings.solvent_id,
        )
