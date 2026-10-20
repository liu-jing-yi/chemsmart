"""
Amber job runner implementation.

This module provides the runner class for executing Amber molecular dynamics
calculations, handling job submission, monitoring, and output processing.
"""

import logging
import os

from chemsmart.jobs.runner import JobRunner

logger = logging.getLogger(__name__)


class AmberJobRunner(JobRunner):
    """
    Job runner specialized for Amber calculations.

    Provides Amber-specific execution logic including proper command
    construction, environment setup, and output handling.
    """

    def __init__(self, **kwargs):
        """
        Initialize Amber job runner.

        Args:
            **kwargs: Arguments passed to parent JobRunner class
        """
        super().__init__(**kwargs)
        self.program = "Amber"

    def run_amber(self, job) -> None:
        """
        Execute an Amber job.

        Args:
            job: AmberJob instance to execute
        """
        logger.info(f"Running Amber job: {job.label}")

        # Construct Amber command
        command = self._build_amber_command(job)

        # Set up environment
        env = self._setup_amber_environment()

        # Execute job
        self._execute_command(command, job.folder, env)

        logger.info(f"Amber job {job.label} execution completed")

    def _build_amber_command(self, job) -> str:
        """
        Build the Amber execution command.

        Args:
            job: AmberJob instance

        Returns:
            Command string to execute
        """
        settings = job.settings

        # Determine appropriate Amber executable
        if settings.simulation_type == "md":
            executable = "pmemd" if self._has_pmemd() else "sander"
        elif settings.simulation_type in ["energy", "opt"]:
            executable = "sander"
        else:
            executable = "sander"  # Default fallback

        # Build command arguments
        args = [
            executable,
            "-O",  # Overwrite output files
            f"-i {job.input_filename}",
            f"-p {job.topology_filename}",
            f"-c {job.coordinate_filename}",
            f"-o {job.output_filename}",
        ]

        # Add trajectory output for MD
        if settings.simulation_type == "md":
            args.extend(
                [
                    f"-x {job.trajectory_filename}",
                    f"-r {job.restart_filename}",
                ]
            )

        # Add parallelization if multiple cores
        if hasattr(self, "num_cores") and self.num_cores > 1:
            if executable == "pmemd":
                # PMEMD uses MPI
                return f"mpirun -np {self.num_cores} {' '.join(args)}"
            else:
                # SANDER can use OpenMP
                args.insert(0, f"export OMP_NUM_THREADS={self.num_cores} &&")

        return " ".join(args)

    def _setup_amber_environment(self) -> dict:
        """
        Set up environment variables for Amber execution.

        Returns:
            Dictionary of environment variables
        """
        env = os.environ.copy()

        # Add common Amber environment variables
        if "AMBERHOME" in env:
            amber_home = env["AMBERHOME"]
            env["PATH"] = f"{amber_home}/bin:{env.get('PATH', '')}"
            env["LD_LIBRARY_PATH"] = (
                f"{amber_home}/lib:{env.get('LD_LIBRARY_PATH', '')}"
            )

        # Set threading for OpenMP if not using MPI
        if not self._has_mpi():
            env["OMP_NUM_THREADS"] = str(getattr(self, "num_cores", 1))

        return env

    def _has_pmemd(self) -> bool:
        """
        Check if PMEMD is available.

        Returns:
            True if PMEMD executable is found
        """
        try:
            import shutil

            return shutil.which("pmemd") is not None
        except Exception:
            return False

    def _has_mpi(self) -> bool:
        """
        Check if MPI is available.

        Returns:
            True if MPI is available
        """
        try:
            import shutil

            return (
                shutil.which("mpirun") is not None
                or shutil.which("mpiexec") is not None
            )
        except Exception:
            return False

    def _execute_command(
        self, command: str, working_dir: str, env: dict
    ) -> None:
        """
        Execute the Amber command.

        Args:
            command: Command string to execute
            working_dir: Directory to run command in
            env: Environment variables
        """
        logger.debug(f"Executing Amber command: {command}")
        logger.debug(f"Working directory: {working_dir}")

        # Use parent class execution logic
        # This would typically interface with the job scheduler/runner
        super().run_command(command, working_dir, env)
