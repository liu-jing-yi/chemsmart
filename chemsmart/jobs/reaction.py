"""Shared Guess → Opt → SP chain for Gaussian and ORCA reaction jobs."""

import logging

from chemsmart.jobs.chain import ChainMixin, JobPhase

logger = logging.getLogger(__name__)

PATH_SEARCH_SINGLE_STRUCTURE = (
    "Path search requires a single reactant geometry and a single product "
    "geometry. Combine fragments into one structure, or provide a TS guess "
    "with extra fragments to optimize without path search."
)


class ReactionChainMixin(ChainMixin):
    """Reaction child-job construction on top of ``ChainMixin``.

    Subclasses set ``_opt_job_class``, ``_ts_job_class``, and
    ``_sp_job_class``, and implement ``_make_guess_job``.

    Path search (QST/NEB) runs only for exactly one reactant geometry and
    one product geometry. Extra fragments with a TS parent skip Guess.
    """

    _opt_job_class = None
    _ts_job_class = None
    _sp_job_class = None
    uses_neb = False

    def __init__(
        self,
        *args,
        reactants=None,
        products=None,
        ts_guess=None,
        no_path_search=False,
        no_sp=False,
        opt_settings=None,
        ts_settings=None,
        sp_settings=None,
        neb_settings=None,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)
        self.reactants = tuple(reactants or ())
        self.products = tuple(products or ())
        self.ts_guess = ts_guess
        self.no_path_search = bool(no_path_search)
        self.no_sp = bool(no_sp)
        self.opt_settings = opt_settings
        self.ts_settings = ts_settings
        self.sp_settings = sp_settings
        self.neb_settings = neb_settings
        self._configure_path_search()

        self.guess_job = (
            self._make_guess_job() if self.uses_path_search else None
        )
        self.reactant_sp_jobs = None
        self.product_sp_jobs = None
        self.ts_sp_job = None
        self.sp_jobs = None
        self._create_opt_jobs()
        self.phases = self._build_reaction_phases()

    def _configure_path_search(self):
        """Skip or reject Guess when fragment counts are not a single R/P pair."""
        if self.no_path_search or not self.products:
            return
        n_reactants = len(self.reactants) if self.reactants else 1
        n_products = len(self.products)
        if n_reactants == 1 and n_products == 1:
            return
        if self.reactants:
            logger.info(
                "Skipping path search for %s: extra fragments with a TS "
                "geometry are optimized, not used as QST/NEB endpoints.",
                self.label,
            )
            self.no_path_search = True
            return
        raise ValueError(PATH_SEARCH_SINGLE_STRUCTURE)

    @property
    def uses_path_search(self):
        return bool(self.products) and not self.no_path_search

    @property
    def path_reactant(self):
        return self.reactants[0] if self.reactants else self.molecule

    @property
    def path_product(self):
        return self.products[0] if self.products else None

    @property
    def path_ts_guess(self):
        if self.ts_guess is not None:
            return self.ts_guess
        if self.reactants:
            return self.molecule
        return None

    @property
    def reactant_molecules(self):
        if self.reactants:
            return self.reactants
        if self.uses_path_search:
            return (self.molecule,)
        return ()

    def _make_guess_job(self):
        raise NotImplementedError

    def _prepare_guess_inputs(self):
        return None

    def _apply_charge_mult(self, settings, molecule):
        if molecule.charge is not None:
            settings.charge = molecule.charge
        if molecule.multiplicity is not None:
            settings.multiplicity = molecule.multiplicity
        return settings

    def _child_settings(self, template, jobtype, freq, molecule):
        base = template if template is not None else self.settings
        settings = base.copy()
        settings.jobtype = jobtype
        settings.freq = freq
        return self._apply_charge_mult(settings, molecule)

    def _ts_settings_for(self, molecule):
        return self._child_settings(self.ts_settings, "ts", True, molecule)

    def _output_molecule(self, job, fallback):
        output = job._output()
        if output is not None and output.normal_termination:
            return output.molecule
        return fallback

    def _ts_molecule(self):
        if self.uses_path_search and self.guess_job is not None:
            return self._output_molecule(self.guess_job, self.molecule)
        return self.molecule

    def _child_job(self, job_class, molecule, settings, label):
        return job_class(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _opt_jobs_for_role(self, molecules, letter):
        n = len(molecules)
        suffixes = (
            (letter,)
            if n == 1
            else tuple(f"{letter}{i}" for i in range(1, n + 1))
        )
        return [
            self._child_job(
                self._opt_job_class,
                molecule,
                self._child_settings(self.opt_settings, "opt", True, molecule),
                f"{self.label}_{suffix}_opt",
            )
            for molecule, suffix in zip(molecules, suffixes)
        ]

    def _set_opt_jobs(self):
        self.opt_jobs = (
            list(self.reactant_opt_jobs)
            + [self.ts_opt_job]
            + list(self.product_opt_jobs)
        )

    def _create_opt_jobs(self):
        self.reactant_opt_jobs = self._opt_jobs_for_role(
            self.reactant_molecules, "R"
        )
        ts_molecule = self._ts_molecule()
        self.ts_opt_job = self._child_job(
            self._ts_job_class,
            ts_molecule,
            self._ts_settings_for(ts_molecule),
            f"{self.label}_TS_opt",
        )
        self.product_opt_jobs = self._opt_jobs_for_role(self.products, "P")
        self._set_opt_jobs()

    def _opt_jobs_for_phase(self):
        if self.uses_path_search:
            ts_molecule = self._ts_molecule()
            self.ts_opt_job = self._child_job(
                self._ts_job_class,
                ts_molecule,
                self._ts_settings_for(ts_molecule),
                f"{self.label}_TS_opt",
            )
            self._set_opt_jobs()
        return self.opt_jobs

    def _make_sp_job(self, opt_job):
        molecule = self._output_molecule(opt_job, opt_job.molecule)
        return self._child_job(
            self._sp_job_class,
            molecule,
            self._child_settings(self.sp_settings, "sp", False, molecule),
            f"{opt_job.label.removesuffix('_opt')}_sp",
        )

    def _sp_jobs_for_phase(self):
        if self.sp_jobs is None:
            self.reactant_sp_jobs = [
                self._make_sp_job(job) for job in self.reactant_opt_jobs
            ]
            self.ts_sp_job = self._make_sp_job(self.ts_opt_job)
            self.product_sp_jobs = [
                self._make_sp_job(job) for job in self.product_opt_jobs
            ]
            self.sp_jobs = (
                list(self.reactant_sp_jobs)
                + [self.ts_sp_job]
                + list(self.product_sp_jobs)
            )
        return self.sp_jobs

    def _build_reaction_phases(self):
        return [
            JobPhase(
                name="Guess",
                jobs=([self.guess_job] if self.guess_job is not None else []),
                skip_if=lambda: not self.uses_path_search,
                before_run=self._prepare_guess_inputs,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message=(
                    "Guess jobs incomplete, halting serial execution."
                ),
            ),
            JobPhase(
                name="Opt",
                jobs_factory=self._opt_jobs_for_phase,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message="Opt jobs incomplete, halting serial execution.",
            ),
            JobPhase(
                name="SP",
                jobs_factory=self._sp_jobs_for_phase,
                skip_if=lambda: self.no_sp,
                require_complete=True,
                stop_on_incomplete=True,
                stop_message="SP jobs incomplete, halting serial execution.",
            ),
        ]
