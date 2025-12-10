"""
Amber Command Line Interface

This module provides the main CLI interface for Amber molecular dynamics
calculations. It defines common options, settings configurations, and the
main Amber command group that serves as the entry point for all Amber-related
operations.
"""

import functools
import logging
import os

import click

from chemsmart.cli.amber.energy import energy

# Import and register subcommands
from chemsmart.cli.amber.md import md
from chemsmart.cli.amber.opt import opt
from chemsmart.cli.job import (
    click_file_label_and_index_options,
    click_filename_options,
    click_pubchem_options,
)
from chemsmart.io.molecules.structure import Molecule
from chemsmart.utils.cli import MyGroup
from chemsmart.utils.io import clean_label
from chemsmart.utils.utils import return_objects_from_string_index

logger = logging.getLogger(__name__)


def click_amber_options(f):
    """
    Common click options decorator for Amber jobs.

    This decorator adds common command-line options that are shared across
    different Amber job types, specifically project settings.
    """

    @click.option(
        "--project", "-p", type=str, default=None, help="Project settings."
    )
    @functools.wraps(f)
    def wrapper_common_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_common_options


def click_amber_settings_options(f):
    """
    Common click options decorator for Amber computational settings.

    This decorator adds comprehensive command-line options for configuring
    Amber molecular dynamics calculations including ensemble, temperature,
    pressure, and simulation parameters.
    """

    @click.option(
        "--ensemble",
        "-e",
        type=click.Choice(["NVE", "NVT", "NPT"], case_sensitive=False),
        default="NVT",
        help="Statistical ensemble for simulation.",
    )
    @click.option(
        "--temperature",
        "-T",
        type=float,
        default=300.0,
        help="Target temperature in Kelvin.",
    )
    @click.option(
        "--pressure",
        "-P",
        type=float,
        default=None,
        help="Target pressure in bar (for NPT ensemble).",
    )
    @click.option(
        "--steps",
        "-s",
        type=int,
        default=50000,
        help="Number of simulation steps.",
    )
    @click.option(
        "--timestep",
        "-dt",
        type=float,
        default=0.002,
        help="Integration timestep in ps.",
    )
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
        help="Initial coordinates file (.inpcrd/.rst7).",
    )
    @functools.wraps(f)
    def wrapper_amber_settings(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_amber_settings


def click_amber_jobtype_options(f):
    """
    Common click options decorator for Amber job type settings.

    This decorator adds options for basic computational parameters like
    job name, runtime, and resource allocation.
    """

    @click.option(
        "--job-name",
        "-jn",
        type=str,
        default="Amber Simulation",
        help="Name for the simulation job.",
    )
    @click.option(
        "--max-runtime",
        "-rt",
        type=str,
        default="48:00:00",
        help="Maximum runtime in HH:MM:SS format.",
    )
    @click.option(
        "--num-nodes",
        "-nn",
        type=int,
        default=1,
        help="Number of compute nodes to request.",
    )
    @click.option("--ppn", type=int, default=16, help="Processors per node.")
    @click.option(
        "--queue",
        "-q",
        type=str,
        default="standard",
        help="Queue/partition name for job submission.",
    )
    @functools.wraps(f)
    def wrapper_jobtype_options(*args, **kwargs):
        return f(*args, **kwargs)

    return wrapper_jobtype_options


@click.group(name="amber", cls=MyGroup)
@click.pass_context
@click_amber_options
@click_filename_options
@click_file_label_and_index_options
@click_pubchem_options
def amber(
    ctx,
    project,
    filename,
    pubchem,
    index,
    label,
):
    """
    Amber molecular dynamics calculations.

    This command group provides access to Amber MD simulations, energy
    calculations, and geometry optimizations. Amber is a suite of programs
    for molecular dynamics simulations of proteins and nucleic acids.

    Examples:
        Run MD simulation:
        $ chemsmart run amber -f system.pdb md -s 100000 -T 300 -e NVT

        Energy calculation:
        $ chemsmart run amber -f molecule.pdb energy

        Geometry optimization:
        $ chemsmart run amber -f structure.pdb opt
    """
    logger.info("Entering Amber command group")

    # Import here to avoid circular imports
    from chemsmart.settings.amber import AmberProjectSettings

    # Initialize context object if needed
    ctx.ensure_object(dict)

    # Get project settings
    if project:
        project_settings = AmberProjectSettings.from_project(project)
    else:
        project_settings = AmberProjectSettings()

    ctx.obj["project_settings"] = project_settings

    # Handle molecule input
    molecules = []

    if pubchem:
        # Load from PubChem
        molecule = Molecule.from_pubchem(pubchem)
        molecules.append(molecule)
        if not label:
            label = clean_label(pubchem)
    elif filename:
        # Load from file
        if index:
            # Load specific indices
            all_molecules = Molecule.from_filename(filename)
            indices = return_objects_from_string_index(all_molecules, index)
            molecules = [all_molecules[i] for i in indices]
        else:
            # Load all molecules
            molecules = Molecule.from_filename(filename)
            if not isinstance(molecules, list):
                molecules = [molecules]

        if not label:
            label = clean_label(
                os.path.splitext(os.path.basename(filename))[0]
            )

    # Store molecules and label in context
    ctx.obj["molecules"] = molecules
    ctx.obj["label"] = label

    # Get job settings from project
    ctx.obj["job_settings"] = project_settings.main_settings()
    ctx.obj["keywords"] = {}

    logger.debug(f"Loaded {len(molecules)} molecule(s)")
    logger.debug(f"Job label: {label}")


# Import and register subcommands will be done after the main command is defined


# Import and register subcommands

amber.add_command(md)
amber.add_command(energy)
amber.add_command(opt)
