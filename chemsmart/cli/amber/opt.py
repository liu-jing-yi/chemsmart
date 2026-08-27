"""
Amber geometry optimization CLI command.

This module provides the command-line interface for running Amber
geometry optimization calculations.
"""

import logging

import click

from chemsmart.cli.amber.amber import (
    amber,
    click_amber_jobtype_options,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@amber.command("opt", cls=MyCommand)
@click_job_options
@click_amber_jobtype_options
@click.option(
    "--topology",
    "-top",
    type=str,
    default=None,
    help="Topology file (.prmtop).",
)
@click.option(
    "--coordinates",
    "-crd",
    type=str,
    default=None,
    help="Coordinates file (.inpcrd/.rst7).",
)
@click.option(
    "--max-cycles",
    "-mc",
    type=int,
    default=1000,
    help="Maximum optimization cycles.",
)
@click.option(
    "--convergence",
    "-conv",
    type=float,
    default=1e-4,
    help="Convergence criterion for optimization.",
)
@click.pass_context
def opt(
    ctx,
    topology,
    coordinates,
    max_cycles,
    convergence,
    job_name,
    max_runtime,
    num_nodes,
    ppn,
    queue,
    **kwargs,
):
    """
    Perform geometry optimization using Amber.

    Optimizes the molecular geometry using Amber force fields to find
    the minimum energy structure.

    Examples:
        Basic optimization:
        $ chemsmart run amber -f molecule.pdb opt

        With custom parameters:
        $ chemsmart run amber -f structure.pdb opt -mc 2000 -conv 1e-6

        With topology files:
        $ chemsmart run amber opt -top system.prmtop -crd system.inpcrd -mc 500

    Args:
        topology: Topology file (.prmtop)
        coordinates: Coordinates file (.inpcrd/.rst7)
        max_cycles: Maximum optimization cycles
        convergence: Convergence criterion
        job_name: Name for the job
        max_runtime: Maximum runtime in HH:MM:SS
        num_nodes: Number of compute nodes
        ppn: Processors per node
        queue: Queue/partition name
    """
    # Get settings from project
    project_settings = ctx.obj["project_settings"]
    opt_project_settings = project_settings.opt_settings()

    # Get job settings and keywords from context
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # Merge project optimization settings with job settings
    opt_settings = opt_project_settings.merge(job_settings, keywords=keywords)

    # Update settings with CLI parameters
    opt_settings.simulation_type = "opt"
    opt_settings.steps = max_cycles
    opt_settings.job_name = job_name
    opt_settings.max_runtime = max_runtime
    opt_settings.num_nodes = num_nodes
    opt_settings.ppn = ppn
    opt_settings.queue = queue

    # Set topology and coordinate files if provided
    if topology:
        opt_settings.topology_file = topology
    if coordinates:
        opt_settings.coordinate_file = coordinates

    # Store convergence criterion (could extend settings class for this)
    opt_settings.convergence = convergence

    # Get label and molecules
    label = ctx.obj["label"]
    molecules = ctx.obj["molecules"]

    # Use first molecule if available
    molecule = molecules[0] if molecules else None

    logger.info(f"Optimization job settings: {opt_settings.__dict__}")

    # Create and return Amber optimization job
    from chemsmart.jobs.amber.job import AmberOptJob

    return AmberOptJob(
        molecule=molecule, settings=opt_settings, label=label, **kwargs
    )
