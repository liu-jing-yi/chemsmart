"""
Amber input file writer.

This module provides functionality for writing Amber input files
including MD input, minimization, and other simulation types.
"""

import logging
import os

logger = logging.getLogger(__name__)


class AmberInputWriter:
    """
    Writer class for Amber input files.

    Handles generation of Amber input files for various simulation types
    including molecular dynamics, energy calculations, and optimization.
    """

    def __init__(self, job):
        """
        Initialize Amber input writer.

        Args:
            job: AmberJob instance containing settings and molecule
        """
        self.job = job
        self.settings = job.settings
        self.molecule = job.molecule

    def write(self, target_directory: str = ".") -> str:
        """
        Write complete Amber input file.

        Args:
            target_directory: Directory to write input files

        Returns:
            Path to main input file
        """
        input_path = os.path.join(target_directory, self.job.input_filename)

        logger.debug(f"Writing Amber input file: {input_path}")

        with open(input_path, "w") as f:
            self._write_header(f)
            self._write_control_section(f)

        logger.info(f"Amber input file written: {input_path}")
        return input_path

    def _write_header(self, f):
        """Write header comment."""
        f.write(f"Amber {self.settings.simulation_type} calculation\n")

    def _write_control_section(self, f):
        """Write the control (&cntrl) section."""
        f.write("&cntrl\n")

        # Basic simulation parameters
        if self.settings.simulation_type == "md":
            f.write("  imin=0,\n")  # Run molecular dynamics
            f.write(f"  nstlim={self.settings.steps},\n")
            f.write("  dt=0.002,\n")  # 2 fs timestep
        elif self.settings.simulation_type == "energy":
            f.write("  imin=0,\n")  # No minimization
            f.write("  nstlim=0,\n")  # No dynamics
            f.write("  maxcyc=1,\n")  # Single energy evaluation
        elif self.settings.simulation_type == "opt":
            f.write("  imin=1,\n")  # Minimization
            f.write(f"  maxcyc={self.settings.steps},\n")
            f.write("  ncyc=100,\n")  # Steepest descent steps

        # Temperature control
        if (
            hasattr(self.settings, "temperature")
            and self.settings.temperature > 0
        ):
            if (
                self.settings.ensemble in ["NVT", "NPT"]
                and self.settings.simulation_type == "md"
            ):
                f.write(f"  temp0={self.settings.temperature},\n")
                f.write("  ntt=3,\n")  # Langevin thermostat
                f.write("  gamma_ln=2.0,\n")

        # Pressure control for NPT
        if (
            hasattr(self.settings, "ensemble")
            and self.settings.ensemble == "NPT"
        ):
            if self.settings.simulation_type == "md":
                f.write("  ntp=1,\n")  # Isotropic pressure scaling
                f.write(f"  pres0={self.settings.pressure or 1.0},\n")
                f.write("  taup=2.0,\n")

        # Output control
        if self.settings.simulation_type == "md":
            f.write("  ntpr=1000,\n")  # Print frequency
            f.write("  ntwx=1000,\n")  # Trajectory write frequency
            f.write("  ntwr=5000,\n")  # Restart file frequency
        elif self.settings.simulation_type == "energy":
            f.write("  ntpr=1,\n")  # Print energy
        elif self.settings.simulation_type == "opt":
            f.write("  ntpr=10,\n")  # Print optimization progress

        # Constraints (if any)
        if hasattr(self.molecule, "constraints") and self.molecule:
            f.write("  ntc=2,\n")  # SHAKE constraints on bonds with hydrogen
            f.write(
                "  ntf=2,\n"
            )  # Do not calculate forces for constrained bonds

        f.write("&end\n")

    def write_topology_file(self, target_directory: str = ".") -> str:
        """
        Write Amber topology file if molecule is available.

        Args:
            target_directory: Directory to write files

        Returns:
            Path to topology file
        """
        if not self.molecule:
            logger.warning("No molecule available for topology generation")
            return None

        topology_path = os.path.join(
            target_directory, self.job.topology_filename
        )

        # This is a placeholder - in practice you would use AmberTools
        # to generate proper topology files
        logger.info(f"Topology file would be generated at: {topology_path}")
        return topology_path

    def write_coordinate_file(self, target_directory: str = ".") -> str:
        """
        Write Amber coordinate file if molecule is available.

        Args:
            target_directory: Directory to write files

        Returns:
            Path to coordinate file
        """
        if not self.molecule:
            logger.warning("No molecule available for coordinate generation")
            return None

        coordinate_path = os.path.join(
            target_directory, self.job.coordinate_filename
        )

        # This is a placeholder - in practice you would convert
        # from the molecule format to Amber coordinate format
        logger.info(
            f"Coordinate file would be generated at: {coordinate_path}"
        )
        return coordinate_path
