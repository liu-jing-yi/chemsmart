"""
ORCA pKa calculation job implementation.

This module provides the ORCApKaJob class for performing pKa
calculations using ORCA with a proper thermodynamic cycle:
1. Gas phase optimization + frequency for both HA and A-
2. Solution phase single point for both HA and A- at the same level of theory

Using the same level of theory ensures proper error cancellation for
solvation free energy calculations.
"""

import os

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import ORCApKaJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.pka import PkaChainMixin


class ORCApKaJob(PkaChainMixin, ORCAJob):
    """
    ORCA job class for pKa calculations using the dual-level proton exchange cycle.

    Performs pKa calculations using the following workflow:
    1. Optimize HA in gas phase (opt + freq)
    2. Optimize A- in gas phase (opt + freq)
    3. Run SP on optimized HA in solution
    4. Run SP on optimized A- in solution
    5. (Optional) Same for reference acid Href and Ref-

    Attributes:
        TYPE (str): Job type identifier ('orcapka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (ORCApKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "orcapka"
    _opt_job_class = ORCAOptJob
    _sp_job_class = ORCASinglePointJob
    _shared_reference_molecule_cache = {}

    def __init__(self, molecule, settings=None, **kwargs):
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

        super().__init__(molecule=molecule, settings=settings, **kwargs)

    @classmethod
    def settings_class(cls):
        return ORCApKaJobSettings

    @property
    def _ref_basename(self):
        """Basename for the reference acid, from the reference geometry file."""
        if not self.settings.has_reference_file:
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

    def _href_opt_label(self):
        return self._ref_basename

    def _ref_opt_label(self):
        return self._ref_conjugate_base_label

    def _href_sp_label(self):
        ref = self._ref_basename
        if ref is None:
            return None
        return f"{ref}_sp"

    def _ref_sp_label(self):
        ref_cb = self._ref_conjugate_base_label
        if ref_cb is None:
            return None
        return f"{ref_cb}_sp"

    def _reference_pka(self):
        return self.settings.reference_pka

    def _subjob_output_paths(self, job, legacy_label=None):
        """Candidate ORCA output files for a pKa sub-job."""
        paths = []
        runner = job.jobrunner
        if runner is not None:
            runner_out = getattr(runner, "job_outputfile", None)
            if runner_out:
                paths.append(runner_out)
        paths.append(job.outputfile)
        paths.append(os.path.join(self.folder, f"{job.label}.out"))
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

    def _finalize_child_job(self, job, legacy_label=None):
        """Keep sub-jobs in the parent folder and resolve scratch/legacy outputs."""
        job.folder = self.folder

        def is_complete():
            return self._subjob_is_complete(job, legacy_label)

        def _output():
            return self._subjob_output(job, legacy_label)

        job.is_complete = is_complete
        job._output = _output
        return job
