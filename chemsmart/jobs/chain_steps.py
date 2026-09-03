"""YAML chain-step registry and ChainJob builder."""

from dataclasses import dataclass

from chemsmart.jobs.chain import ChainJob, JobPhase
from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
from chemsmart.jobs.crest.job import CRESTJob
from chemsmart.jobs.gaussian.fukui import GaussianFukuiJob
from chemsmart.jobs.gaussian.irc import GaussianIRCJob
from chemsmart.jobs.gaussian.modred import GaussianModredJob
from chemsmart.jobs.gaussian.nci import GaussianNCIJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.qrc import GaussianQRCJob
from chemsmart.jobs.gaussian.resp import GaussianRESPJob
from chemsmart.jobs.gaussian.scan import GaussianScanJob
from chemsmart.jobs.gaussian.settings import GaussianTDDFTJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.tddft import GaussianTDDFTJob
from chemsmart.jobs.gaussian.ts import GaussianTSJob
from chemsmart.jobs.gaussian.wbi import GaussianWBIJob
from chemsmart.jobs.orca.fukui import ORCAFukuiJob
from chemsmart.jobs.orca.irc import ORCAIRCJob
from chemsmart.jobs.orca.modred import ORCAModredJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.qrc import ORCAQRCJob
from chemsmart.jobs.orca.scan import ORCAScanJob
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.orca.ts import ORCATSJob
from chemsmart.jobs.xtb.hess import XTBHessJob
from chemsmart.jobs.xtb.opt import XTBOptJob
from chemsmart.jobs.xtb.singlepoint import XTBSinglePointJob
from chemsmart.settings.chain import ChainProjectSettings
from chemsmart.settings.crest import CRESTProjectSettings
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.orca import ORCAProjectSettings
from chemsmart.settings.xtb import XTBProjectSettings

__all__ = [
    "ChainStepSpec",
    "CHAIN_PROGRAM_SETTINGS",
    "CHAIN_STEP_SPECS",
    "CHAIN_NESTED_ONLY_JOBS",
    "get_chain_step_spec",
    "molecule_from_completed_job",
    "build_chain_job",
]

_RESP_ROUTE = (
    "HF/6-31+G(d) SCF=Tight Pop=MK IOp(6/33=2,6/41=10,6/42=17,6/50=1)"
)

# CLI jobs that need extra structures or option surfaces. Run them as
# nested program commands, e.g. ``chain -p NAME gaussian pka ...``.
CHAIN_NESTED_ONLY_JOBS = {
    "gaussian": (
        "pka",
        "reaction",
        "qmmm",
        "com",
        "userjob",
        "link",
        "traj",
        "dias",
        "crest",
    ),
    "orca": ("pka", "reaction", "qmmm", "inp", "neb"),
}


@dataclass(frozen=True)
class ChainStepSpec:
    """Resolved YAML chain step: job class and project-settings class."""

    job_class: type
    project_settings_class: type


_GAUSSIAN_ORCA_JOBS = (
    ("opt", GaussianOptJob, ORCAOptJob),
    ("ts", GaussianTSJob, ORCATSJob),
    ("sp", GaussianSinglePointJob, ORCASinglePointJob),
    ("modred", GaussianModredJob, ORCAModredJob),
    ("irc", GaussianIRCJob, ORCAIRCJob),
    ("scan", GaussianScanJob, ORCAScanJob),
    ("fukui", GaussianFukuiJob, ORCAFukuiJob),
    ("qrc", GaussianQRCJob, ORCAQRCJob),
)

_CHAIN_PROGRAMS = {
    "crest": (
        CRESTProjectSettings,
        {"conformers": CRESTConformerSearchJob},
    ),
    "xtb": (
        XTBProjectSettings,
        {
            "opt": XTBOptJob,
            "sp": XTBSinglePointJob,
            "hess": XTBHessJob,
        },
    ),
    "gaussian": (
        GaussianProjectSettings,
        {
            **{
                name: gaussian for name, gaussian, _orca in _GAUSSIAN_ORCA_JOBS
            },
            "nci": GaussianNCIJob,
            "wbi": GaussianWBIJob,
            "td": GaussianTDDFTJob,
            "resp": GaussianRESPJob,
        },
    ),
    "orca": (
        ORCAProjectSettings,
        {name: orca for name, _gaussian, orca in _GAUSSIAN_ORCA_JOBS},
    ),
}

CHAIN_PROGRAM_SETTINGS = {
    program: settings_class
    for program, (settings_class, _jobs) in _CHAIN_PROGRAMS.items()
}
CHAIN_STEP_SPECS = {
    program: jobs
    for program, (_settings_class, jobs) in _CHAIN_PROGRAMS.items()
}


def _gaussian_td_settings(project_settings):
    settings = project_settings.td_settings()
    if settings is None:
        settings = project_settings.sp_settings()
    if isinstance(settings, GaussianTDDFTJobSettings):
        return settings
    return GaussianTDDFTJobSettings(**settings.__dict__)


def _gaussian_resp_settings(project_settings):
    settings = project_settings.sp_settings().copy()
    settings.route_to_be_written = _RESP_ROUTE
    return settings


_JOB_SETTINGS_ALIASES = {"fukui": "sp", "qrc": "opt"}
_JOB_SETTINGS = {
    "conformers": lambda settings: settings.conformer_settings(),
    "opt": lambda settings: settings.opt_settings(),
    "ts": lambda settings: settings.ts_settings(),
    "sp": lambda settings: settings.sp_settings(),
    "hess": lambda settings: settings.hess_settings(),
    "modred": lambda settings: settings.modred_settings(),
    "irc": lambda settings: settings.irc_settings(),
    "scan": lambda settings: settings.scan_settings(),
    "nci": lambda settings: settings.nci_settings(),
    "wbi": lambda settings: settings.wbi_settings(),
}
_GAUSSIAN_JOB_SETTINGS = {
    "td": _gaussian_td_settings,
    "resp": _gaussian_resp_settings,
}


def _settings_for_chain_step(program, job, project_settings):
    if program == "gaussian":
        gaussian_settings = _GAUSSIAN_JOB_SETTINGS.get(job)
        if gaussian_settings is not None:
            return gaussian_settings(project_settings)
    name = _JOB_SETTINGS_ALIASES.get(job, job)
    settings_for = _JOB_SETTINGS.get(name)
    if settings_for is None:
        raise ValueError(
            f"No project settings method for chain step {program} {job!r}."
        )
    return settings_for(project_settings)


def get_chain_step_spec(program, job):
    """Return the registry entry for a YAML ``(program, job)`` pair.

    Args:
        program (str): Chain program name (``crest``, ``xtb``,
            ``gaussian``, or ``orca``).
        job (str): YAML job name (e.g. ``opt``, ``irc``, ``fukui``).

    Returns:
        ChainStepSpec: Job class and project-settings class for the pair.

    Raises:
        ValueError: If the pair is not a YAML chain step. Jobs that need
            extra structures or option surfaces (pKa, reaction, QM/MM)
            must be run as nested program commands.
    """
    jobs = CHAIN_STEP_SPECS.get(program)
    job_class = None if jobs is None else jobs.get(job)
    if job_class is not None:
        return ChainStepSpec(
            job_class=job_class,
            project_settings_class=CHAIN_PROGRAM_SETTINGS[program],
        )
    nested_jobs = CHAIN_NESTED_ONLY_JOBS.get(program, ())
    if job in nested_jobs:
        raise ValueError(
            f"Chain YAML steps do not support {program} {job!r}; "
            f"run it as a nested command "
            f"(chain -p NAME {program} {job} ...)."
        )
    allowed = ", ".join(jobs or ())
    if allowed:
        raise ValueError(
            f"Unsupported chain step job {job!r} for program "
            f"{program!r}. Supported jobs: {allowed}."
        )
    allowed_programs = ", ".join(ChainProjectSettings.PROGRAMS)
    raise ValueError(
        f"Unsupported chain step program {program!r}. "
        f"Allowed programs: {allowed_programs}."
    )


def molecule_from_completed_job(job):
    """Geometry to pass from a completed child to the next chain step.

    CREST uses the best conformer. Other jobs use
    ``optimized_structure()``, falling back to the output molecule when
    the child did not optimize (e.g. single-point).
    """
    if not job.is_complete():
        return None
    if isinstance(job, CRESTJob):
        output = job._output()
        return None if output is None else output.best_conformer
    mol = job.optimized_structure()
    if mol is not None:
        return mol
    output = job._output()
    return None if output is None else output.molecule


def _first_set(*values):
    for value in values:
        if value is not None:
            return value
    return None


def _charge_multiplicity_overrides(molecule, charge, multiplicity):
    mol_charge = None if molecule is None else molecule.charge
    mol_mult = None if molecule is None else molecule.multiplicity
    overrides = {}
    for key, value in (
        ("charge", _first_set(charge, mol_charge)),
        ("multiplicity", _first_set(multiplicity, mol_mult)),
    ):
        if value is not None:
            overrides[key] = value
    return overrides


def _molecule_for_step(index, phases, initial_molecule):
    if index == 0:
        return initial_molecule
    previous_jobs = phases[index - 1].resolve_jobs()
    if not previous_jobs or not all(
        job.is_complete() for job in previous_jobs
    ):
        return None
    return molecule_from_completed_job(previous_jobs[0])


def build_chain_job(
    chain_settings,
    molecule,
    label,
    charge=None,
    multiplicity=None,
    jobrunner=None,
    skip_completed=True,
    **kwargs,
):
    """Build a ``ChainJob`` from YAML pipeline steps.

    Each step becomes a ``JobPhase`` whose factory reads geometry from the
    previous completed child. Charge and multiplicity from ``charge`` /
    ``multiplicity`` (and otherwise from ``molecule``) are merged into
    each child with keywords ``("charge", "multiplicity")``.

    Args:
        chain_settings: Loaded chain project settings with ``steps``.
        molecule: Input structure for the first step.
        label (str): Chain label; children are
            ``{label}_{program}_{job}``.
        charge (int, optional): Chain-level charge override.
        multiplicity (int, optional): Chain-level multiplicity override.
        jobrunner: Runner assigned to the parent ``ChainJob``.
        skip_completed (bool): Forwarded to child jobs.

    Returns:
        ChainJob: Phased pipeline ready for ``run`` / ``sub``.

    Raises:
        ValueError: If there are no steps, ``molecule`` is missing, or a
            step's ``(program, job)`` pair is not a YAML chain step.
    """
    if not chain_settings.steps:
        raise ValueError(
            f"Chain project {chain_settings.PROJECT_NAME!r} has no steps."
        )
    if molecule is None:
        raise ValueError("Chain pipeline requires a molecule.")

    prepared = []
    project_settings_by_program = {}
    for step in chain_settings.steps:
        spec = get_chain_step_spec(step.program, step.job)
        if step.program not in project_settings_by_program:
            project_name = chain_settings.project_for(step.program)
            project_settings_by_program[step.program] = (
                spec.project_settings_class.from_project(project_name)
            )
        prepared.append((step, spec))

    overrides = _charge_multiplicity_overrides(molecule, charge, multiplicity)
    phases = []

    def make_factory(index):
        def factory():
            mol = _molecule_for_step(index, phases, molecule)
            if mol is None:
                return []
            step, spec = prepared[index]
            settings = _settings_for_chain_step(
                step.program,
                step.job,
                project_settings_by_program[step.program],
            )
            if overrides:
                settings = settings.merge(overrides, keywords=tuple(overrides))
            child_label = f"{label}_{step.program}_{step.job}"
            return [
                spec.job_class(
                    molecule=mol,
                    settings=settings,
                    label=child_label,
                    jobrunner=None,
                    skip_completed=skip_completed,
                )
            ]

        return factory

    for index, (step, _spec) in enumerate(prepared):
        phase_name = f"{step.program}_{step.job}"
        phases.append(
            JobPhase(
                name=phase_name,
                jobs_factory=make_factory(index),
                require_complete=True,
                stop_on_incomplete=True,
                stop_message=(
                    f"{phase_name} incomplete, halting serial execution."
                ),
            )
        )

    return ChainJob(
        molecule=molecule,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        phases=phases,
        **kwargs,
    )
