import logging

import click

logger = logging.getLogger(__name__)

CHAIN_PROJECT_SETTINGS_KEY = "chain_project_settings"


def resolve_program_project(ctx, program, project):
    """Resolve a program's project name for CLI group callbacks.

    An explicit ``-p`` on the program group is used as-is. Otherwise, if
    chain project settings are present on ``ctx.obj``, that file's alias
    for ``program`` is used. A missing alias raises ``click.UsageError``.
    Without chain settings, ``project`` is returned unchanged so
    ``from_project`` keeps its current behavior.
    """
    if project is not None:
        return project
    if ctx.obj is None or CHAIN_PROJECT_SETTINGS_KEY not in ctx.obj:
        return project
    chain_settings = ctx.obj[CHAIN_PROJECT_SETTINGS_KEY]
    try:
        return chain_settings.project_for(program)
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


def build_jobs(ctx, job_cls, settings, skip_completed, kwargs):
    """Build one or more jobs from molecules stored on the Click context."""
    jobrunner = ctx.obj["jobrunner"]
    molecules = ctx.obj["molecules"]
    molecule_indices = ctx.obj["molecule_indices"]
    label = ctx.obj["label"]

    if len(molecules) > 1 and molecule_indices is not None:
        jobs = []
        for molecule, idx in zip(molecules, molecule_indices):
            # Preserve one output directory per selected structure.
            molecule_label = f"{label}_idx{idx}"
            logger.info(f"Creating job {molecule_label}")
            jobs.append(
                job_cls(
                    molecule=molecule,
                    settings=settings,
                    label=molecule_label,
                    jobrunner=jobrunner,
                    skip_completed=skip_completed,
                    **kwargs,
                )
            )
        return jobs

    molecule = molecules[-1]
    logger.info(f"Creating job {label}")
    return job_cls(
        molecule=molecule,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
        skip_completed=skip_completed,
        **kwargs,
    )
