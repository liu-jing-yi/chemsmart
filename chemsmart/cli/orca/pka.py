"""
CLI for ORCA pKa input generation (job submission).

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

from chemsmart.cli.job import click_job_options
from chemsmart.cli.orca.orca import orca
from chemsmart.cli.pka import (
    batch_pka_jobs,
    click_pka_proton_options,
    click_pka_shared_options,
    configure_pka_submission,
    ensure_pka_label_suffix,
    submit_pka_jobs,
)
from chemsmart.jobs.orca.settings import ORCApKaJobSettings
from chemsmart.utils.cli import MyCommand, MyGroup


def _orca_pka_settings(proton_index, shared, opt_settings, _project_settings):
    return ORCApKaJobSettings.build_orca_pka_settings(
        proton_index, shared, opt_settings
    )


def _orca_job_class():
    from chemsmart.jobs.orca.pka import ORCApKaJob

    return ORCApKaJob


@orca.group("pka", cls=MyGroup, invoke_without_command=True)
@click_job_options
@click_pka_shared_options
@click_pka_proton_options
@click.pass_context
def pka(ctx, **kwargs):
    """ORCA pKa job submission.

    \b
    Subcommands:
      submit         Single-molecule job submission (default).
      batch          Table-driven batch submission.

    When invoked without a subcommand the ``submit`` path runs
    automatically.

    \b
    For output analysis (backend-independent):
      chemsmart run pka analyze ...
      chemsmart run pka batch-analyze ...
      chemsmart run pka thermo ...

    \b
    Thermodynamic cycles:
      proton exchange (default): HA + Ref- -> A- + HRef
      direct: uses G_soln(H+) in water
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
    """Submit a single-molecule ORCA pKa calculation.

    \b
    Examples:
      chemsmart run orca -f acid.xyz -c 0 -m 1 pka -pi 10 \\
          -r ref.xyz -rpi 1 -rc 0 -rm 1 submit

      chemsmart run orca -f acid.xyz -c 0 -m 1 pka -pi 10 \\
          -s direct submit
    """
    return submit_pka_jobs(
        ctx,
        skip_completed,
        proton_index,
        color_code,
        job_class=_orca_job_class(),
        settings_builder=_orca_pka_settings,
        label_for=ensure_pka_label_suffix,
        batch_command=batch,
        program="ORCA",
        **kwargs,
    )


@pka.command("batch", cls=MyCommand)
@click_job_options
@click_pka_proton_options
@click.pass_context
def batch(ctx, skip_completed, proton_index, color_code, **kwargs):
    """Batch ORCA pKa job submission from a CSV table or multi-molecule CDXML.

    \b
    CSV table format (4 columns, whitespace or comma-delimited):
        filepath    proton_index    charge    multiplicity

    CDXML files create one job per fragment using coloured-proton detection.

    \b
    Examples:
      chemsmart run orca -p myproject -f molecules.txt pka \\
          -s "proton exchange" -r ref.xyz -rpi 5 -rc 0 -rm 1 batch
    """
    return batch_pka_jobs(
        ctx,
        skip_completed,
        proton_index,
        color_code,
        job_class=_orca_job_class(),
        settings_builder=_orca_pka_settings,
        label_for=ensure_pka_label_suffix,
        submit_command=submit,
        program="ORCA",
        **kwargs,
    )
