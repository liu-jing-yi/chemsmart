"""
ORCA pKa calculation job implementation.

This module provides the ORCApKaJob class for performing pKa
calculations using ORCA with a proper thermodynamic cycle:
1. Optional CREST conformational sampling of HA, A-, and any reference acid
2. Gas phase optimization + frequency for each selected conformer of HA and A-
3. Solution phase single point for each optimized conformer at the same level of theory

Using the same level of theory ensures proper error cancellation for
solvation free energy calculations.
"""

import logging
import os

from chemsmart.cli.pka import build_pka_crest_job, select_crest_conformers
from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import ORCApKaJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.runner import decide_phase_transition, run_phase_jobs
from chemsmart.utils.datasets import pka_job_species_outputs, pka_subjob_label

logger = logging.getLogger(__name__)


class ORCApKaJob(ORCAJob):
    """
    ORCA job class for pKa calculations using the dual-level proton exchange cycle.

    Performs pKa calculations using the following workflow:
    1. Optionally run CREST conformational sampling on HA, A-, and,
       when a reference acid is set, HRef and Ref-
    2. Optimize HA in gas phase (opt + freq)
    3. Optimize A- in gas phase (opt + freq)
    4. Run SP on optimized HA in solution
    5. Run SP on optimized A- in solution
    6. (Optional) Same for reference acid Href and Ref-

    When sampling is enabled, each sampled species yields N gas-phase
    opt+freq jobs and N solvent single-point jobs, one pair per conformer.
    Parent completion requires DFT opt+SP only.

    Attributes:
        TYPE (str): Job type identifier ('orcapka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (ORCApKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcapka"
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

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        skip_completed=True,
        **kwargs,
    ):
        if not isinstance(settings, ORCApKaJobSettings):
            raise ValueError(
                f"Settings must be instance of ORCApKaJobSettings, "
                f"but got {type(settings).__name__} instead!"
            )

        if settings.proton_index is None:
            raise ValueError(
                "proton_index must be specified in ORCApKaJobSettings "
                "to identify which proton to remove for the conjugate base."
            )

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            **kwargs,
        )

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

    # ------------------------------------------------------------------
    # Basename helpers for label derivation
    # ------------------------------------------------------------------

    @property
    def _ref_basename(self):
        """Basename for the reference acid (Href), derived from the
        reference geometry filename so it stays unique when multiple HA
        share one Href."""
        if not self.has_reference_jobs:
            return None
        return os.path.splitext(
            os.path.basename(self.settings.reference_file)
        )[0]

    @property
    def _ref_conjugate_base_label(self):
        """Label for the reference conjugate base (Ref⁻)."""
        ref = self._ref_basename
        if ref is None:
            return None
        return f"{ref}_cb"

    @classmethod
    def settings_class(cls):
        return ORCApKaJobSettings

    # ------------------------------------------------------------------
    # Molecule properties
    # ------------------------------------------------------------------

    @property
    def protonated_molecule(self):
        """Get the protonated molecule (HA)."""
        protonated_mol, _ = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        return protonated_mol

    @property
    def conjugate_base_molecule(self):
        """Get the conjugate base molecule (A-)."""
        _, conjugate_base_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        return conjugate_base_mol

    @property
    def reference_molecule(self):
        """Get the reference acid molecule (Href)."""
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            return self.settings.get_reference_molecule()
        return reference_pair[0]

    @property
    def reference_conjugate_base_molecule(self):
        """Get the reference conjugate base molecule (Ref-)."""
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            return self.settings.get_reference_conjugate_base_molecule()
        return reference_pair[1]

    # ------------------------------------------------------------------
    # Job preparation
    # ------------------------------------------------------------------

    def _subjob_output_paths(self, job, legacy_label=None):
        """Candidate ORCA output files for a pKa sub-job."""
        job.folder = self.folder
        paths = []
        runner = job.jobrunner
        if runner is not None:
            runner_out = getattr(runner, "job_outputfile", None)
            if runner_out:
                paths.append(runner_out)
        paths.append(job.outputfile)
        if legacy_label is not None:
            paths.append(os.path.join(self.folder, f"{legacy_label}.out"))
        seen = set()
        ordered = []
        for path in paths:
            if path and path not in seen:
                seen.add(path)
                ordered.append(path)
        return ordered

    def _subjob_is_complete(self, job, legacy_label=None):
        from chemsmart.io.orca.output import ORCAOutput

        for path in self._subjob_output_paths(job, legacy_label):
            if not path or not os.path.exists(path):
                continue
            try:
                if ORCAOutput(path).normal_termination:
                    return True
            except Exception:
                continue
        return False

    def _subjob_output(self, job, legacy_label=None):
        from chemsmart.io.orca.output import ORCAOutput

        for path in self._subjob_output_paths(job, legacy_label):
            if not os.path.exists(path):
                continue
            try:
                output = ORCAOutput(path)
            except Exception:
                continue
            if output.normal_termination:
                return output
        return None

    def _bind_subjob(self, job, legacy_label=None):
        """Keep sub-jobs in the parent folder and resolve scratch/legacy outputs."""
        job.folder = self.folder
        parent = self

        def is_complete():
            job.folder = parent.folder
            return parent._subjob_is_complete(job, legacy_label)

        job.is_complete = is_complete

    def _sync_subjob_folders(self, jobs):
        for job in jobs or []:
            job.folder = self.folder

    def _prepare_pka_jobs(self):
        """Prepare optimization jobs for target and reference acids."""
        prot_mol, conj_mol = self.settings.conjugate_pair_molecules(
            self.molecule
        )
        if self.settings.sampling:
            self.protonated_crest_job = build_pka_crest_job(
                prot_mol,
                f"{self.label}_HA_crest",
                self.settings,
                self,
            )
            self.conjugate_base_crest_job = build_pka_crest_job(
                conj_mol,
                f"{self.label}_A_crest",
                self.settings,
                self,
            )
            self.crest_jobs = [
                self.protonated_crest_job,
                self.conjugate_base_crest_job,
            ]
            n = self.settings.num_conformers
            ha_molecules = [prot_mol] * n
            a_molecules = [conj_mol] * n
        else:
            ha_molecules = [prot_mol]
            a_molecules = [conj_mol]

        self._prepare_target_opt_jobs(ha_molecules, a_molecules)
        self._create_sp_jobs()

        if self.has_reference_jobs:
            href_mol, ref_mol = self._reference_pair_molecules()
            if self.settings.sampling:
                self.ref_acid_crest_job = build_pka_crest_job(
                    href_mol,
                    f"{self._ref_basename}_crest",
                    self.settings,
                    self,
                )
                self.ref_conjugate_base_crest_job = build_pka_crest_job(
                    ref_mol,
                    f"{self._ref_conjugate_base_label}_crest",
                    self.settings,
                    self,
                )
                self.crest_jobs.extend(
                    [
                        self.ref_acid_crest_job,
                        self.ref_conjugate_base_crest_job,
                    ]
                )
                n = self.settings.num_conformers
                href_molecules = [href_mol] * n
                ref_molecules = [ref_mol] * n
            else:
                href_molecules = [href_mol]
                ref_molecules = [ref_mol]
            self._prepare_ref_opt_jobs(href_molecules, ref_molecules)
            self._create_ref_sp_jobs()

    def _reference_pair_molecules(self):
        """Return the HRef and Ref- molecules for the configured reference acid."""
        reference_pair = self._get_cached_reference_pair(self.settings)
        if reference_pair is None:
            return self.settings.reference_pair_molecules()
        return reference_pair

    def _ref_conformer_label(self, base_label, index, num_conformers):
        if num_conformers > 1:
            return f"{base_label}_c{index}"
        return base_label

    def _prepare_target_opt_jobs(self, ha_molecules, a_molecules):
        """Build gas-phase opt+freq jobs for HA and A- conformers."""
        protonated_settings, conjugate_base_settings = (
            self.settings.conjugate_pair_job_settings(self.molecule)
        )
        self.protonated_opt_jobs = self._make_species_opt_jobs(
            ha_molecules,
            protonated_settings,
            "HA",
            legacy_label=self.label,
        )
        self.conjugate_base_opt_jobs = self._make_species_opt_jobs(
            a_molecules,
            conjugate_base_settings,
            "A",
            legacy_label=f"{self.label}_cb",
        )
        self.protonated_job = self.protonated_opt_jobs[0]
        self.conjugate_base_job = self.conjugate_base_opt_jobs[0]
        self.opt_jobs = self.protonated_opt_jobs + self.conjugate_base_opt_jobs
        self.sp_jobs = None
        self.protonated_sp_jobs = None
        self.conjugate_base_sp_jobs = None
        self.protonated_sp_job = None
        self.conjugate_base_sp_job = None

    def _make_species_opt_jobs(
        self, molecules, settings, species, legacy_label=None
    ):
        num_conformers = len(molecules)
        jobs = []
        bind_legacy = legacy_label if num_conformers == 1 else None
        for index, molecule in enumerate(molecules, start=1):
            job = ORCAOptJob(
                molecule=molecule,
                settings=settings,
                label=pka_subjob_label(
                    self.label, species, "opt", index, num_conformers
                ),
                jobrunner=self.jobrunner,
                skip_completed=self.skip_completed,
            )
            self._bind_subjob(job, legacy_label=bind_legacy)
            jobs.append(job)
        return jobs

    def _create_sp_jobs(self):
        """Create one solution-phase SP job per optimized conformer."""
        protonated_sp_settings, conjugate_base_sp_settings = (
            self.settings._create_solution_phase_sp_settings(self.molecule)
        )
        self.protonated_sp_jobs = self._make_species_sp_jobs(
            self.protonated_opt_jobs,
            protonated_sp_settings,
            "HA",
            opt_legacy_label=self.label,
            sp_legacy_label=f"{self.label}_sp",
        )
        self.conjugate_base_sp_jobs = self._make_species_sp_jobs(
            self.conjugate_base_opt_jobs,
            conjugate_base_sp_settings,
            "A",
            opt_legacy_label=f"{self.label}_cb",
            sp_legacy_label=f"{self.label}_cb_sp",
        )
        self.protonated_sp_job = self.protonated_sp_jobs[0]
        self.conjugate_base_sp_job = self.conjugate_base_sp_jobs[0]
        self.sp_jobs = self.protonated_sp_jobs + self.conjugate_base_sp_jobs

    def _make_species_sp_jobs(
        self,
        opt_jobs,
        settings,
        species,
        opt_legacy_label=None,
        sp_legacy_label=None,
    ):
        num_conformers = len(opt_jobs)
        jobs = []
        opt_legacy = opt_legacy_label if num_conformers == 1 else None
        sp_legacy = sp_legacy_label if num_conformers == 1 else None
        for index, opt_job in enumerate(opt_jobs, start=1):
            jobs.append(
                self._make_sp_job(
                    opt_job,
                    opt_job.molecule,
                    settings,
                    pka_subjob_label(
                        self.label, species, "sp", index, num_conformers
                    ),
                    opt_legacy_label=opt_legacy,
                    sp_legacy_label=sp_legacy,
                )
            )
        return jobs

    def _prepare_ref_opt_jobs(self, href_molecules, ref_molecules):
        """Build gas-phase opt+freq jobs for HRef and Ref- conformers."""
        ref_acid_settings, ref_cb_settings = (
            self.settings.reference_pair_job_settings()
        )
        self.ref_acid_opt_jobs = self._make_ref_opt_jobs(
            href_molecules, ref_acid_settings, self._ref_basename
        )
        self.ref_conjugate_base_opt_jobs = self._make_ref_opt_jobs(
            ref_molecules, ref_cb_settings, self._ref_conjugate_base_label
        )
        self.ref_acid_job = self.ref_acid_opt_jobs[0]
        self.ref_conjugate_base_job = self.ref_conjugate_base_opt_jobs[0]
        self.ref_opt_jobs = (
            self.ref_acid_opt_jobs + self.ref_conjugate_base_opt_jobs
        )
        self.ref_sp_jobs = None
        self.ref_acid_sp_jobs = None
        self.ref_conjugate_base_sp_jobs = None
        self.ref_acid_sp_job = None
        self.ref_conjugate_base_sp_job = None

    def _make_ref_opt_jobs(self, molecules, settings, base_label):
        num_conformers = len(molecules)
        legacy_label = base_label if num_conformers == 1 else None
        jobs = []
        for index, molecule in enumerate(molecules, start=1):
            job = ORCAOptJob(
                molecule=molecule,
                settings=settings,
                label=self._ref_conformer_label(
                    base_label, index, num_conformers
                ),
                jobrunner=self.jobrunner,
                skip_completed=self.skip_completed,
            )
            self._bind_subjob(job, legacy_label=legacy_label)
            jobs.append(job)
        return jobs

    def _create_ref_sp_jobs(self):
        """Create one reference solvent SP job per optimized conformer."""
        ref_acid_sp_settings, ref_cb_sp_settings = (
            self.settings.reference_pair_sp_job_settings()
        )
        self.ref_acid_sp_jobs = self._make_ref_sp_jobs(
            self.ref_acid_opt_jobs,
            ref_acid_sp_settings,
            f"{self._ref_basename}_sp",
            opt_legacy_label=self._ref_basename,
        )
        self.ref_conjugate_base_sp_jobs = self._make_ref_sp_jobs(
            self.ref_conjugate_base_opt_jobs,
            ref_cb_sp_settings,
            f"{self._ref_conjugate_base_label}_sp",
            opt_legacy_label=self._ref_conjugate_base_label,
        )
        self.ref_acid_sp_job = self.ref_acid_sp_jobs[0]
        self.ref_conjugate_base_sp_job = self.ref_conjugate_base_sp_jobs[0]
        self.ref_sp_jobs = (
            self.ref_acid_sp_jobs + self.ref_conjugate_base_sp_jobs
        )

    def _make_ref_sp_jobs(
        self, opt_jobs, settings, base_label, opt_legacy_label=None
    ):
        num_conformers = len(opt_jobs)
        opt_legacy = opt_legacy_label if num_conformers == 1 else None
        sp_legacy = base_label if num_conformers == 1 else None
        jobs = []
        for index, opt_job in enumerate(opt_jobs, start=1):
            jobs.append(
                self._make_sp_job(
                    opt_job,
                    opt_job.molecule,
                    settings,
                    self._ref_conformer_label(
                        base_label, index, num_conformers
                    ),
                    opt_legacy_label=opt_legacy,
                    sp_legacy_label=sp_legacy,
                )
            )
        return jobs

    # ------------------------------------------------------------------
    # Execution
    # ------------------------------------------------------------------

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
        """Run solution phase single point jobs using optimized geometries."""
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
        """Run reference solution phase single point jobs."""
        if not self.has_reference_jobs:
            return
        if not self._ref_opt_jobs_are_complete():
            logger.warning(
                "Reference optimization jobs not complete. Cannot run reference SP jobs."
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

    def _make_sp_job(
        self,
        opt_job,
        fallback_molecule,
        sp_settings,
        sp_label,
        opt_legacy_label=None,
        sp_legacy_label=None,
    ):
        """Create SP job using optimized geometry if available."""
        out = self._subjob_output(opt_job, legacy_label=opt_legacy_label)
        if out is not None:
            mol = out.molecule
        else:
            mol = fallback_molecule

        sp_job = ORCASinglePointJob(
            molecule=mol,
            settings=sp_settings,
            label=sp_label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )
        self._bind_subjob(sp_job, legacy_label=sp_legacy_label)
        return sp_job

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
        """
        Execute the pKa calculation.

        Optionally runs CREST sampling, then gas phase optimization jobs
        followed by solution phase SP jobs.
        """
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
        """
        Check if all pKa jobs are complete.

        Returns:
            bool: True if all optimization jobs and SP jobs
                have completed successfully (including reference jobs if provided).
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

    # ------------------------------------------------------------------
    # Output accessors
    # ------------------------------------------------------------------

    @property
    def protonated_output(self):
        return self.protonated_job._output()

    @property
    def conjugate_base_output(self):
        return self.conjugate_base_job._output()

    @property
    def protonated_sp_output(self):
        return self.protonated_sp_job._output()

    @property
    def conjugate_base_sp_output(self):
        return self.conjugate_base_sp_job._output()

    @property
    def ref_acid_output(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_job._output()

    @property
    def ref_conjugate_base_output(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_job._output()

    @property
    def ref_acid_sp_output(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_acid_sp_job._output()

    @property
    def ref_conjugate_base_sp_output(self):
        if not self.has_reference_jobs:
            return None
        return self.ref_conjugate_base_sp_job._output()

    # ------------------------------------------------------------------
    # Thermochemistry
    # ------------------------------------------------------------------

    def _opt_jobs_are_complete(self):
        """Return True when both target gas-phase optimization jobs finished."""
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
        """Return True when both target solution-phase SP jobs finished."""
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

    def _pka_output_files(self):
        """Return gas and solvent output paths for every pKa species.

        Each species maps to ``{"gas": [...], "solv": [...]}`` in conformer
        order. Reference species are included when a reference acid is set.
        """
        return pka_job_species_outputs(self)

    def compute_thermochemistry(self):
        """Compute and return thermochemistry results for all species."""
        from chemsmart.io.orca.output import ORCApKaOutput

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot compute thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        files = self._pka_output_files()

        return ORCApKaOutput.compute_pka_thermochemistry(
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
        """Print formatted thermochemistry summary to stdout."""
        from chemsmart.io.orca.output import ORCApKaOutput

        if not self._opt_jobs_are_complete():
            raise ValueError(
                "Cannot print thermochemistry: optimization jobs are not complete. "
                "Run the pKa jobs first using job.run()."
            )

        files = self._pka_output_files()

        return ORCApKaOutput.print_pka_summary(
            ha_gas_file=files["HA"]["gas"],
            a_gas_file=files["A-"]["gas"],
            href_gas_file=files["HRef"]["gas"] if "HRef" in files else None,
            ref_gas_file=files["Ref-"]["gas"] if "Ref-" in files else None,
            ha_solv_file=files["HA"]["solv"],
            a_solv_file=files["A-"]["solv"],
            href_solv_file=files["HRef"]["solv"] if "HRef" in files else None,
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
