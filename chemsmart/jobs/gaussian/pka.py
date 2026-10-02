"""
Gaussian pKa calculation job implementation.

GaussianpKaJob runs the shared pKa phase machine: optional CREST sampling,
gas-phase optimization plus frequency, then a solution-phase single point
for each selected conformer. Child jobs are GaussianOptJob and
GaussianSinglePointJob.
"""

from chemsmart.jobs.chain.pka import PKaJob
from chemsmart.jobs.gaussian.job import GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob


class GaussianpKaJob(PKaJob, GaussianJob):
    """
    Gaussian job class for pKa calculations using direct thermodynamic cycle.

    Performs pKa calculations using the following workflow:
    1. Optionally run CREST conformational sampling on HA, A-, and,
       when a reference acid is set, HRef and Ref-
    2. Optimize HA in gas phase (opt + freq) - get G(HA)_gas
    3. Optimize A- in gas phase (opt + freq) - get G(A-)_gas
    4. Run SP on optimized HA in solution - get E(HA)_aq
    5. Run SP on optimized A- in solution - get E(A-)_aq
    6. Calculate solvation free energies and pKa

    When sampling is enabled, each sampled species yields N gas-phase
    opt+freq jobs and N solvent single-point jobs, one pair per conformer.
    Parent completion requires DFT opt+SP only. Solvent jobs are rebuilt
    from the optimized geometries when the SP phase runs.

    Attributes:
        TYPE (str): Job type identifier ('g16pka').
        molecule (Molecule): Protonated molecular structure (HA).
        settings (GaussianpKaJobSettings): pKa calculation configuration.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    TYPE = "g16pka"
    opt_job_class = GaussianOptJob
    sp_job_class = GaussianSinglePointJob

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        skip_completed=True,
        **kwargs,
    ):
        if not isinstance(settings, GaussianpKaJobSettings):
            raise ValueError(
                "Settings must be instance of GaussianpKaJobSettings for "
                f"{self.__class__.__name__}, but is {settings} instead!"
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
    def original_mol(self):
        """Original molecule used to initialize the job (usually HA)."""
        return self.molecule

    @property
    def conjugate_base_mol(self):
        """Conjugate base molecule (A-)."""
        if self.conjugate_base_job:
            return self.conjugate_base_job.molecule
        _, conj_mol = self.settings.conjugate_pair_molecules(self.molecule)
        return conj_mol
