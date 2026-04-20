"""
Amber molecular dynamics CLI command.

This module provides the command-line interface for running Amber
molecular dynamics simulations.
"""

import logging

import click

from chemsmart.cli.amber.amber import (
    amber,
    click_amber_jobtype_options,
    click_amber_settings_options,
)
from chemsmart.cli.job import click_job_options
from chemsmart.utils.cli import MyCommand

logger = logging.getLogger(__name__)


@amber.command("md", cls=MyCommand)
@click_job_options
@click_amber_settings_options
@click_amber_jobtype_options
@click.pass_context
def md(
    ctx,
    ensemble,
    temperature,
    pressure,
    steps,
    timestep,
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
    Run Amber molecular dynamics simulation.

    Performs molecular dynamics simulation using Amber with the specified
    ensemble, temperature, pressure, and time parameters.

    Examples:
        NVT simulation:
        $ chemsmart run amber -f system.pdb md -s 100000 -T 300 -e NVT

        NPT simulation:
        $ chemsmart run amber -f system.pdb md -s 50000 -T 310 -P 1.0 -e NPT

        With custom topology files:
        $ chemsmart run amber md -top system.prmtop -crd system.inpcrd -s 200000

    Args:
        ensemble: Statistical ensemble (NVE, NVT, NPT)
        temperature: Target temperature in Kelvin
        pressure: Target pressure in bar (for NPT)
        steps: Number of simulation steps
        timestep: Integration timestep in ps
        topology: Topology file (.prmtop)
        coordinates: Initial coordinates (.inpcrd/.rst7)
        job_name: Name for the simulation job
        max_runtime: Maximum runtime in HH:MM:SS
        num_nodes: Number of compute nodes
        ppn: Processors per node
        queue: Queue/partition name
    """
    # Get settings from project
    project_settings = ctx.obj["project_settings"]
    md_project_settings = project_settings.md_settings()

    # Get job settings and keywords from context
    job_settings = ctx.obj["job_settings"]
    keywords = ctx.obj["keywords"]

    # Merge project MD settings with job settings
    md_settings = md_project_settings.merge(job_settings, keywords=keywords)

    # Update settings with CLI parameters
    md_settings.simulation_type = "md"
    md_settings.ensemble = ensemble
    md_settings.temperature = temperature
    if pressure is not None:
        md_settings.pressure = pressure
    md_settings.steps = steps
    md_settings.job_name = job_name
    md_settings.max_runtime = max_runtime
    md_settings.num_nodes = num_nodes
    md_settings.ppn = ppn
    md_settings.queue = queue

    # Set topology and coordinate files if provided
    if topology:
        md_settings.topology_file = topology
    if coordinates:
        md_settings.coordinate_file = coordinates

    # Get label and molecules
    label = ctx.obj["label"]
    molecules = ctx.obj["molecules"]

    # Use first molecule if available
    molecule = molecules[0] if molecules else None

    logger.info(f"MD job settings: {md_settings.__dict__}")

    # Create and return Amber MD job
    from chemsmart.jobs.amber.job import AmberMDJob

    return AmberMDJob(
        molecule=molecule, settings=md_settings, label=label, **kwargs
    )
