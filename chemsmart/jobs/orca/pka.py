"""
ORCA pKa calculation job implementation.

ORCApKaJob runs the shared pKa phase machine: optional CREST sampling,
gas-phase optimization plus frequency, then a solution-phase single point
for each selected conformer. Reference-acid labels stay derived from the
reference geometry filename so one HRef can be shared by several HA jobs.
Child jobs are ORCAOptJob and ORCASinglePointJob.
"""

import logging
import os

from chemsmart.jobs.chain.pka import PKaJob
from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import ORCApKaJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.utils.datasets import pka_subjob_label

logger = logging.getLogger(__name__)


class ORCApKaJob(PKaJob, ORCAJob):
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
    Parent completion requires DFT opt+SP only. Solvent jobs are rebuilt
    from the optimized geometries when the SP phase runs.

    Attributes:
        TYPE (str): Job type identifier ('orcapka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (ORCApKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcapka"
    opt_job_class = ORCAOptJob
    sp_job_class = ORCASinglePointJob

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
                "Settings must be instance of ORCApKaJobSettings, "
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
        self._reset_pka_child_jobs()

    @property
    def _ref_basename(self):
        """Basename of the reference geometry, shared when several HA use one HRef."""
        if not self.has_reference_jobs:
            return None
        return os.path.splitext(
            os.path.basename(self.settings.reference_file)
        )[0]

    @property
    def _ref_conjugate_base_label(self):
        """Label stem for the reference conjugate base (Ref-)."""
        ref = self._ref_basename
        if ref is None:
            return None
        return f"{ref}_cb"

    @classmethod
    def settings_class(cls):
        return ORCApKaJobSettings

    def _pka_output_class(self):
        from chemsmart.io.orca.output import ORCApKaOutput

        return ORCApKaOutput

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

    def _crest_job_label(self, species):
        if species == "HRef":
            return f"{self._ref_basename}_crest"
        if species == "Ref":
            return f"{self._ref_conjugate_base_label}_crest"
        return super()._crest_job_label(species)

    def _subjob_label(self, species, stage, index, num_conformers):
        if species == "HRef":
            base = self._ref_basename
        elif species == "Ref":
            base = self._ref_conjugate_base_label
        else:
            return pka_subjob_label(
                self.label, species, stage, index, num_conformers
            )
        if stage == "sp":
            base = f"{base}_sp"
        if num_conformers > 1:
            return f"{base}_c{index}"
        return base

    def _subjob_legacy_label(self, species, stage):
        if species == "HA":
            if stage == "opt":
                return self.label
            return f"{self.label}_sp"
        if species == "A":
            if stage == "opt":
                return f"{self.label}_cb"
            return f"{self.label}_cb_sp"
        if species == "HRef":
            base = self._ref_basename
        elif species == "Ref":
            base = self._ref_conjugate_base_label
        else:
            return None
        if stage == "sp":
            return f"{base}_sp"
        return base

    def _optimized_molecule_from_job(
        self, job, fallback_molecule, legacy_label=None
    ):
        out = self._subjob_output(job, legacy_label=legacy_label)
        if out is not None:
            return out.molecule
        return fallback_molecule

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
