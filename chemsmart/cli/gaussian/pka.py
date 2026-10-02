"""
CLI for Gaussian pKa input generation (job submission).

Subcommands
-----------
submit         Submit single-molecule (or multi-fragment CDXML) pKa jobs.
batch          Table-driven batch job submission.

When ``pka`` is invoked without an explicit subcommand the ``submit``
path is executed automatically for backward compatibility.

Output analysis (analyze, batch-analyze, thermo) lives in the
backend-independent ``chemsmart run pka`` command.
"""

import click

from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.job import click_job_options
from chemsmart.cli.pka import (
    batch_pka_jobs,
    click_pka_proton_options,
    click_pka_shared_options,
    configure_pka_submission,
    submit_pka_jobs,
)
from chemsmart.jobs.gaussian.pka import GaussianpKaJob
from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings
from chemsmart.utils.cli import MyCommand, MyGroup


def _gaussian_pka_settings(
    proton_index, shared, opt_settings, project_settings
):
    return GaussianpKaJobSettings.build_gaussian_pka_settings(
        proton_index,
        shared,
        opt_settings,
        project_settings.sp_settings(),
    )


def _gaussian_pka_label(raw):
    return raw


@gaussian.group("pka", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_pka_shared_options
@click_pka_proton_options
@click.pass_context
def pka(ctx, **kwargs):
    """Gaussian pKa job submission.

    For output analysis (backend-independent), use:
      chemsmart run pka analyze ...
      chemsmart run pka batch-analyze ...
    """
    return configure_pka_submission(
        ctx,
        submit_command=submit,
        batch_command=batch,
        **kwargs,
    )


@pka.command("submit", cls=MyCommand)
@click_job_options
@click_pka_proton_options
@click.pass_context
def submit(ctx, skip_completed, proton_index, color_code, **kwargs):
    """Submit a single-molecule Gaussian pKa calculation.

    Builds Gaussian pKa settings from CLI/project context, validates required
    inputs (charge/multiplicity), and returns the created job(s). If the input
    is a multi-fragment CDXML, this will expand to one job per molecule.

    Example:
        chemsmart sub gaussian <gaussian_options> pka <pka_options> submit

    Args:
        ctx: Click context containing shared options and objects from the
            parent `pka` command (settings, molecules, jobrunner).
        skip_completed (bool): Skip execution for jobs already completed.
        **kwargs: Additional CLI options forwarded to job creation.

    Returns:
        GaussianpKaJob: Single job when one molecule is processed.
        list[GaussianpKaJob]: List of jobs when multiple molecules detected.
    """
    return submit_pka_jobs(
        ctx,
        skip_completed,
        proton_index,
        color_code,
        job_class=GaussianpKaJob,
        settings_builder=_gaussian_pka_settings,
        label_for=_gaussian_pka_label,
        batch_command=batch,
        program="Gaussian",
        **kwargs,
    )


@pka.command("batch", cls=MyCommand)
@click_job_options
@click_pka_proton_options
@click.pass_context
def batch(ctx, skip_completed, proton_index, color_code, **kwargs):
    """Batch pKa job submission from a CSV table or multi-molecule CDXML.

    CSV tables provide filepath, proton_index, charge, and multiplicity per row.
    CDXML files create one job per fragment using coloured-proton detection.
    """
    return batch_pka_jobs(
        ctx,
        skip_completed,
        proton_index,
        color_code,
        job_class=GaussianpKaJob,
        settings_builder=_gaussian_pka_settings,
        label_for=_gaussian_pka_label,
        submit_command=submit,
        program="Gaussian",
        **kwargs,
    )
