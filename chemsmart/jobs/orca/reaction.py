"""
ORCA reaction workflow: optional NEB-TS guess, then TS/R/P opt and SP.

Guess reuses ``ORCANEBJob`` with project ``neb_settings()`` (default
``NEB-TS``). Child Opt/TS/SP jobs are stock classes sequenced with
``ChainMixin``.
"""

import os

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.neb import ORCANEBJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.settings import (
    ORCAJobSettings,
    ORCANEBJobSettings,
    ORCATSJobSettings,
)
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.orca.ts import ORCATSJob
from chemsmart.jobs.reaction import ReactionChainMixin

_TS_ONLY_KEYS = (
    "inhess",
    "inhess_filename",
    "hybrid_hess",
    "hybrid_hess_atoms",
    "numhess",
    "recalc_hess",
    "trust_radius",
    "tssearch_type",
    "scants_modred",
    "full_scan",
)
_NEB_ONLY_KEYS = (
    "joboption",
    "nimages",
    "ending_xyzfile",
    "intermediate_xyzfile",
    "restarting_xyzfile",
    "preopt_ends",
)


def _settings_as(settings, cls, drop_keys):
    """Copy ``settings`` as ``cls``, dropping attributes from other jobs."""
    if isinstance(settings, cls):
        return settings.copy()
    data = dict(settings.copy().__dict__)
    for key in drop_keys:
        data.pop(key, None)
    return cls(**data)


class ORCAReactionJob(ReactionChainMixin, ORCAJob):
    """ORCA R/TS/P chain: optional NEB-TS guess, then opt+freq and SP.

    Case 1 (no product, or ``no_path_search``): skip Guess. Parent
    ``molecule`` is the TS. Optional reactants and products are extra
    minima.

    Case 2 (products given): run ``ORCANEBJob`` (reactant start, product
    ``ending_xyzfile``, optional TS guess ``intermediate_xyzfile``), then
    the same case-1 characterization. Parent ``molecule`` is the reactant
    unless ``reactants`` is also set, in which case parent ``molecule`` is
    the NEB intermediate.

    Attributes:
        TYPE (str): Job type identifier ('orcareaction').
        molecule (Molecule): TS (case 1) or reactant (case 2).
        settings (ORCAJobSettings): Parent settings used as the default
            child template when opt/ts/sp/neb settings are omitted.
        label (str): Base job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the jobs.
        skip_completed (bool): If True, completed jobs are not rerun.
        reactants (tuple): Extra reactant minima to optimize.
        products (tuple): Product minima; presence selects case 2.
        ts_guess (Molecule): Optional NEB intermediate when ``-f`` is the
            reactant.
        no_path_search (bool): Force case 1 even when products are set.
        no_sp (bool): Omit solution-phase single-point children.
    """

    TYPE = "orcareaction"
    _opt_job_class = ORCAOptJob
    _ts_job_class = ORCATSJob
    _sp_job_class = ORCASinglePointJob

    def __init__(
        self,
        molecule,
        settings=None,
        label=None,
        jobrunner=None,
        skip_completed=True,
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
        if not isinstance(settings, ORCAJobSettings):
            raise ValueError(
                f"Settings must be instance of ORCAJobSettings for "
                f"{self.__class__.__name__}, but is {settings} instead!"
            )

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            skip_completed=skip_completed,
            reactants=reactants,
            products=products,
            ts_guess=ts_guess,
            no_path_search=no_path_search,
            no_sp=no_sp,
            opt_settings=opt_settings,
            ts_settings=ts_settings,
            sp_settings=sp_settings,
            neb_settings=neb_settings,
            **kwargs,
        )

    def _ts_settings_for(self, molecule):
        settings = super()._ts_settings_for(molecule)
        if isinstance(settings, ORCATSJobSettings):
            return settings
        converted = _settings_as(settings, ORCATSJobSettings, _NEB_ONLY_KEYS)
        converted.jobtype = "ts"
        converted.freq = True
        return self._apply_charge_mult(converted, molecule)

    def _neb_settings_for(self, molecule):
        if self.neb_settings is not None:
            base = self.neb_settings
        elif self.opt_settings is not None:
            base = self.opt_settings
        else:
            base = self.settings
        settings = _settings_as(base, ORCANEBJobSettings, _TS_ONLY_KEYS)
        settings.jobtype = "neb"
        settings.freq = False
        if settings.joboption is None:
            settings.joboption = "NEB-TS"
        return self._apply_charge_mult(settings, molecule)

    def _make_guess_job(self):
        settings = self._neb_settings_for(self.path_reactant)
        settings.ending_xyzfile = os.path.join(
            self.folder, f"{self.label}_neb_P.xyz"
        )
        if self.path_ts_guess is not None:
            settings.intermediate_xyzfile = os.path.join(
                self.folder, f"{self.label}_neb_TS.xyz"
            )
        return ORCANEBJob(
            molecule=self.path_reactant,
            settings=settings,
            label=f"{self.label}_neb",
            jobrunner=self.jobrunner,
            skip_completed=self.skip_completed,
        )

    def _prepare_guess_inputs(self):
        if self.guess_job is None:
            return
        settings = self.guess_job.settings
        product_path = settings.ending_xyzfile
        if product_path:
            self.path_product.write_xyz(product_path, mode="w")
        ts_path = settings.intermediate_xyzfile
        if ts_path and self.path_ts_guess is not None:
            self.path_ts_guess.write_xyz(ts_path, mode="w")
