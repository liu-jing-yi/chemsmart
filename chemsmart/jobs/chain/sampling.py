"""CREST conformer selection used by pKa jobs."""

import logging
import os

logger = logging.getLogger(__name__)


def build_pka_crest_job(molecule, label, settings, parent_job):
    """Build a CREST conformer-search job for one pKa species.

    CREST settings are loaded from ``settings.crest_project`` when that YAML
    exists. Otherwise default conformer settings are used. Charge and
    multiplicity come from ``molecule``; ``nprocs`` comes from the parent
    jobrunner.

    Args:
        molecule: HA or A- ``Molecule`` to sample.
        label (str): CREST job label (for example ``{label}_HA_crest``).
        settings: pKa job settings exposing ``crest_project``.
        parent_job: Parent Gaussian or ORCA pKa job.

    Returns:
        CRESTConformerSearchJob: Job with a typed CREST runner.
    """
    from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
    from chemsmart.jobs.runner import JobRunner

    parent_runner = parent_job.jobrunner
    crest_settings = _pka_crest_job_settings(molecule, settings, parent_runner)
    job = CRESTConformerSearchJob(
        molecule=molecule,
        settings=crest_settings,
        label=label,
        jobrunner=None,
        skip_completed=parent_job.skip_completed,
    )
    job.jobrunner = JobRunner.from_job(
        job,
        server=parent_runner.server,
        fake=parent_runner.FAKE,
        num_cores=parent_runner.num_cores,
    )
    return job


def select_crest_conformers(crest_job, num_conformers, fallback_molecule):
    """Select CREST geometries for DFT, or fall back to the input molecule.

    If CREST has not written an output file yet, return ``None`` so the
    caller can halt and wait (HPC resubmit). After CREST has terminated
    (normally or abnormally), always return a non-empty list: selected
    conformers, a shorter available set, or ``[fallback_molecule]``.

    ``N == 1`` uses ``crest_best.xyz`` (else the first frame of
    ``crest_conformers.xyz``). ``N > 1`` uses the N lowest frames of
    energy-sorted ``crest_conformers.xyz``. Charge and multiplicity are
    copied from ``fallback_molecule``.

    Args:
        crest_job: CREST job whose folder is parsed.
        num_conformers (int): Number of conformers requested (``N >= 1``).
        fallback_molecule: Input HA or A- geometry used if extraction fails.

    Returns:
        list or None: Selected molecules, or ``None`` if CREST has not
        finished.
    """
    if num_conformers is None or num_conformers < 1:
        num_conformers = 1

    try:
        if not os.path.exists(crest_job.outputfile):
            logger.info(
                f"CREST job {crest_job.label} has not finished; waiting."
            )
            return None
        return _select_pka_crest_conformers(
            crest_job, num_conformers, fallback_molecule
        )
    except Exception as exc:
        logger.warning(
            f"Failed to extract CREST conformers from {crest_job.label}: "
            f"{exc}. Using the input geometry."
        )
        return [_copy_pka_molecule_charge(fallback_molecule)]


def _default_pka_crest_settings():
    from chemsmart.jobs.crest.settings import CRESTJobSettings

    settings = CRESTJobSettings.default()
    settings.jobtype = "conformers"
    return settings


def _pka_crest_job_settings(molecule, settings, parent_runner):
    from chemsmart.settings.crest import CRESTProjectSettings

    crest_project = settings.crest_project
    crest_settings = None
    if crest_project is not None:
        try:
            crest_settings = CRESTProjectSettings.from_project(
                crest_project
            ).conformer_settings()
        except FileNotFoundError:
            logger.warning(
                f"No CREST project settings found for {crest_project!r}; "
                "using default CREST conformer settings."
            )

    if crest_settings is None:
        crest_settings = _default_pka_crest_settings()
    else:
        crest_settings = crest_settings.copy()
        if crest_settings.jobtype is None:
            crest_settings.jobtype = "conformers"

    crest_settings.charge = molecule.charge
    crest_settings.multiplicity = molecule.multiplicity
    if parent_runner is not None and parent_runner.num_cores is not None:
        crest_settings.nprocs = parent_runner.num_cores
    return crest_settings


def _select_pka_crest_conformers(crest_job, num_conformers, fallback_molecule):
    from chemsmart.io.crest.output import CRESTOutput

    output = CRESTOutput(folder=crest_job.folder)
    if not output.normal_termination:
        logger.warning(
            f"CREST job {crest_job.label} did not terminate normally. "
            "Using available geometries or the input structure."
        )

    if num_conformers == 1:
        selected = output.best_conformer
        if selected is None:
            logger.warning(
                f"CREST job {crest_job.label} produced no usable geometry. "
                "Using the input structure."
            )
            return [_copy_pka_molecule_charge(fallback_molecule)]
        return [_copy_pka_molecule_charge(selected, fallback_molecule)]

    conformers = list(output.conformers)
    if not conformers and output.best_conformer is not None:
        conformers = [output.best_conformer]

    if not conformers:
        logger.warning(
            f"CREST job {crest_job.label} produced no conformers. "
            "Using the input structure."
        )
        return [_copy_pka_molecule_charge(fallback_molecule)]

    if len(conformers) < num_conformers:
        logger.warning(
            f"CREST job {crest_job.label} produced {len(conformers)} "
            f"conformer(s); requested {num_conformers}. Using the available "
            "set."
        )

    return [
        _copy_pka_molecule_charge(molecule, fallback_molecule)
        for molecule in conformers[:num_conformers]
    ]


def _copy_pka_molecule_charge(molecule, charge_source=None):
    copied = molecule.copy()
    source = charge_source if charge_source is not None else molecule
    copied.charge = source.charge
    copied.multiplicity = source.multiplicity
    return copied
