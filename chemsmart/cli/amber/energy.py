"""
Amber energy calculation CLI command.

This module provides the command-line interface for running Amber
single-point energy calculations.
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


@amber.command("energy", cls=MyCommand)
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
@click.pass_context
def energy(
    ctx,
    topology,
    coordinates,
    job_name,
    max_runtime,
    num_nodes,
    ppn,
    queue,
    **kwargs,
):
    """
    Calculate single-point energy using Amber.

    Performs a single-point energy calculation using Amber force fields
    without any dynamics or optimization.

    Examples:
        Basic energy calculation:
        $ chemsmart run amber -f molecule.pdb energy

        With topology files:
        $ chemsmart run amber energy -top system.prmtop -crd system.inpcrd

    Args:
        topology: Topology file (.prmtop)
        coordinates: Coordinates file (.inpcrd/.rst7)
        job_name: Name for the job
        max_runtime: Maximum runtime in HH:MM:SS
        num_nodes: Number of compute nodes
        ppn: Processors per node
        queue: Queue/partition name
    """
    # Get settings from project
    project_settings = ctx.obj["project_settings"]
    energy_project_settings = project_settings.energy_settings()

    # Get job settings and keywords from context
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # Merge project energy settings with job settings
    energy_settings = energy_project_settings.merge(
        job_settings, keywords=keywords
    )

    # Update settings with CLI parameters
    energy_settings.simulation_type = "energy"
    energy_settings.steps = 0  # No dynamics for energy calculation
    energy_settings.job_name = job_name
    energy_settings.max_runtime = max_runtime
    energy_settings.num_nodes = num_nodes
    energy_settings.ppn = ppn
    energy_settings.queue = queue

    # Set topology and coordinate files if provided
    if topology:
        energy_settings.topology_file = topology
    if coordinates:
        energy_settings.coordinate_file = coordinates

    # Get label and molecules
    label = ctx.obj["label"]
    molecules = ctx.obj["molecules"]

    # Use first molecule if available
    molecule = molecules[0] if molecules else None

    logger.info(f"Energy job settings: {energy_settings.__dict__}")

    # Create and return Amber energy job
    from chemsmart.jobs.amber.job import AmberEnergyJob

    return AmberEnergyJob(
        molecule=molecule, settings=energy_settings, label=label, **kwargs
    )
