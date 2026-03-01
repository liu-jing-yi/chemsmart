"""
Amber job implementation.

This module contains the main Amber job classes for running molecular dynamics
and energy calculations using the Amber program package.
"""

import logging
import os
from typing import Optional

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.amber.settings import AmberJobSettings
from chemsmart.jobs.job import Job
from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class AmberJob(Job):
    """
    Base Amber job class.

    This class provides the foundation for all Amber molecular dynamics
    calculations including setup, execution, and output handling.

    Attributes:
        PROGRAM (str): Program identifier ('Amber').
        TYPE (str): Job type identifier.
        molecule (Molecule): Molecular structure used for the calculation.
        settings (AmberJobSettings): Configuration options for the job.
        label (str): Job identifier used for file naming.
        jobrunner (JobRunner): Execution backend that runs the job.
        skip_completed (bool): If True, completed jobs are not rerun.
    """

    PROGRAM = "Amber"
    TYPE = "amber"

    def __init__(
        self,
        molecule: Optional[Molecule] = None,
        settings: Optional[AmberJobSettings] = None,
        label: Optional[str] = None,
        jobrunner: Optional[JobRunner] = None,
        **kwargs,
    ):
        """
        Initialize AmberJob.

        Args:
            molecule: Molecule object for the calculation
            settings: AmberJobSettings instance
            label: Job label for identification
            jobrunner: Job runner instance
            **kwargs: Additional arguments passed to parent class
        """
        if settings is None:
            settings = AmberJobSettings()

        # Set default label if not provided
        if label is None:
            if molecule is not None:
                label = f"amber_{molecule.formula}"
            else:
                label = "amber_job"

        super().__init__(
            molecule=molecule,
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )

        # Validate settings
        if hasattr(self.settings, "validate_settings"):
            self.settings.validate_settings()

        logger.debug(f"Created AmberJob with label: {self.label}")

    @property
    def input_filename(self) -> str:
        """Get the main Amber input filename."""
        return f"{self.label}.in"

    @property
    def output_filename(self) -> str:
        """Get the main Amber output filename."""
        return f"{self.label}.out"

    @property
    def topology_filename(self) -> str:
        """Get the topology filename."""
        return f"{self.label}.prmtop"

    @property
    def coordinate_filename(self) -> str:
        """Get the coordinate filename."""
        return f"{self.label}.inpcrd"

    @property
    def restart_filename(self) -> str:
        """Get the restart filename."""
        return f"{self.label}.rst7"

    @property
    def trajectory_filename(self) -> str:
        """Get the trajectory filename."""
        return f"{self.label}.mdcrd"

    def write_input(self, target_directory: str = ".") -> str:
        """
        Write Amber input file.

        Args:
            target_directory: Directory to write input files

        Returns:
            Path to written input file
        """
        input_path = os.path.join(target_directory, self.input_filename)

        with open(input_path, "w") as f:
            f.write(self._generate_input_content())

        logger.info(f"Written Amber input file: {input_path}")
        return input_path

    def _generate_input_content(self) -> str:
        """
        Generate Amber input file content.

        Returns:
            String content for Amber input file
        """
        settings = self.settings

        # Basic input structure for Amber
        lines = [
            f"Amber {settings.simulation_type} simulation",
            "&cntrl",
            f"  nstlim={settings.steps},",
            "  dt=0.002,",  # 2 fs timestep
        ]

        # Temperature control
        if settings.ensemble in ["NVT", "NPT"]:
            lines.extend(
                [
                    f"  temp0={settings.temperature},",
                    "  ntt=3,",  # Langevin thermostat
                    "  gamma_ln=2.0,",
                ]
            )

        # Pressure control for NPT
        if settings.ensemble == "NPT":
            lines.extend(
                [
                    "  ntp=1,",  # Isotropic pressure scaling
                    f"  pres0={settings.pressure or 1.0},",
                    "  taup=2.0,",
                ]
            )

        # Output control
        lines.extend(
            [
                "  ntpr=1000,",  # Print frequency
                "  ntwx=1000,",  # Trajectory write frequency
                "  ntwr=5000,",  # Restart file frequency
            ]
        )

        lines.append("&end")

        return "\n".join(lines)

    def run(self, **kwargs) -> None:
        """
        Execute the Amber job.

        Args:
            **kwargs: Additional arguments for job execution
        """
        logger.info(f"Running Amber job: {self.label}")

        # Write input files
        from chemsmart.jobs.amber.writer import AmberInputWriter

        writer = AmberInputWriter(self)
        writer.write(self.folder)

        # Set up execution command
        if hasattr(self.jobrunner, "run_amber"):
            # Use specialized Amber runner if available
            self.jobrunner.run_amber(self)
        else:
            # Use generic runner
            super().run(**kwargs)

        logger.info(f"Amber job {self.label} completed")

    @classmethod
    def from_files(
        cls,
        topology_file: str,
        coordinate_file: str,
        settings: Optional[AmberJobSettings] = None,
        label: Optional[str] = None,
        jobrunner: Optional[JobRunner] = None,
        **kwargs,
    ) -> "AmberJob":
        """
        Create AmberJob from topology and coordinate files.

        Args:
            topology_file: Path to topology (.prmtop) file
            coordinate_file: Path to coordinate (.inpcrd/.rst7) file
            settings: AmberJobSettings instance
            label: Job label
            jobrunner: Job runner instance
            **kwargs: Additional arguments

        Returns:
            Configured AmberJob instance
        """
        if settings is None:
            settings = AmberJobSettings()

        settings.topology_file = topology_file
        settings.coordinate_file = coordinate_file

        if label is None:
            label = os.path.splitext(os.path.basename(topology_file))[0]

        return cls(
            molecule=None,  # Will be loaded from files
            settings=settings,
            label=label,
            jobrunner=jobrunner,
            **kwargs,
        )


class AmberMDJob(AmberJob):
    """
    Specialized Amber job for molecular dynamics simulations.
    """

    TYPE = "amber_md"

    def __init__(self, **kwargs):
        """Initialize Amber MD job."""
        super().__init__(**kwargs)

        # Ensure simulation type is set to MD
        if self.settings.simulation_type != "md":
            self.settings.simulation_type = "md"


class AmberEnergyJob(AmberJob):
    """
    Specialized Amber job for single-point energy calculations.
    """

    TYPE = "amber_energy"

    def __init__(self, **kwargs):
        """Initialize Amber energy job."""
        super().__init__(**kwargs)

        # Set up for single-point energy calculation
        self.settings.simulation_type = "energy"
        self.settings.steps = 0  # No dynamics


class AmberOptJob(AmberJob):
    """
    Specialized Amber job for geometry optimization.
    """

    TYPE = "amber_opt"

    def __init__(self, **kwargs):
        """Initialize Amber optimization job."""
        super().__init__(**kwargs)

        # Set up for geometry optimization
        self.settings.simulation_type = "opt"
        self.settings.steps = 1000  # Optimization steps
