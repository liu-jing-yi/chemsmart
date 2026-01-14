import inspect
import logging
from abc import abstractmethod
from typing import Optional

from chemsmart.settings.executable import (
    AmberExecutable,
    GaussianExecutable,
    NCIPLOTExecutable,
    ORCAExecutable,
)
from chemsmart.settings.user import ChemsmartUserSettings
from chemsmart.utils.mixins import RegistryMixin

user_settings = ChemsmartUserSettings()


logger = logging.getLogger(__name__)


class RunScript:
    """
    Script generator for computational job execution.

    Creates Python scripts that handle job execution with proper environment
    setup and command line argument passing. Manages the execution context
    for computational chemistry jobs.

    Attributes:
        filename (str): Path to the output script file.
        batch (bool): Whether this is a batch job execution.
        cli_args: Command line arguments to pass to the job.
    """

    def __init__(self, filename, cli_args, batch=False):
        """
        Initialize the run script generator.

        Args:
            filename (str): Path where the script will be written.
            cli_args: Command line arguments for job execution.
            batch (bool): Whether this is a batch job. Defaults to False.
        """
        self.filename = filename
        self.batch = batch
        self.cli_args = cli_args

    def write(self):
        """
        Write the run script to the specified file.

        Creates a Python script file that can be executed to run the
        computational job with the specified arguments.
        """
        with open(self.filename, "w") as f:
            self._write(f)

    def _write(self, f):
        """
        Write the script contents to the file handle.

        Generates a Python script with proper environment setup and
        job execution commands.

        Args:
            f: File handle to write the script contents to.
        """
        contents = f"""\
        #!/usr/bin/env python
        import os
        os.environ['OMP_NUM_THREADS'] = '1'
        
        from chemsmart.cli.run import run

        def run_job():
            run({self.cli_args!r})

        if __name__ == '__main__':
            run_job()
        """

        # Needed to remove leading whitespace in the docstring
        contents = inspect.cleandoc(contents)

        f.write(contents)


class Submitter(RegistryMixin):
    """
    Abstract base class for job submission systems.

    Provides the foundation for scheduler-specific job submitters that handle
    the creation and submission of computational chemistry jobs to various
    cluster management systems.

    Attributes:
        NAME (str): Class-level identifier for the submitter type.
        name (str): Instance identifier for this submitter (often same as NAME).
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission parameters passed through to subclasses.
    """

    NAME: Optional[str] = None

    def __init__(self, name, job, server, **kwargs):
        """
        Initialize the job submitter.

        Args:
            name (str): Name identifier for this submitter instance.
            job: Job instance to be submitted.
            server: Server configuration for submission.
            **kwargs: Additional submission parameters.
        """
        self.name = name
        self.job = job
        self.server = server
        self.kwargs = kwargs

    def __str__(self):
        """
        String representation of the submitter.

        Returns:
            str: Human-readable submitter description.
        """
        return f"Submitter: {self.name}"

    def __eq__(self, other):
        """
        Check equality based on submitter name.

        Args:
            other (Submitter): Another submitter instance to compare.

        Returns:
            bool: True if submitter names are equal.
        """
        return self.name == other.name

    def __hash__(self):
        """
        Generate hash based on submitter name.

        Returns:
            int: Hash value for submitter name.
        """
        return hash(self.name)

    def __repr__(self):
        """
        Developer representation of the submitter.

        Returns:
            str: Detailed submitter representation for debugging.
        """
        return f"Submitter(name={self.name})"

    def __call__(self):
        """
        Create a submitter instance based on the name.

        Returns:
            Submitter: Configured submitter instance.

        Raises:
            ValueError: If no submitter is defined for the specified name.
        """
        submitter_cls = [
            s for s in Submitter.subclasses() if self.name == s.NAME
        ]
        if len(submitter_cls) == 0:
            raise ValueError(
                f"No submitter of defined name: {self.name}.\n"
                f"Available submitters: {Submitter.subclasses()}"
            )

        assert len(submitter_cls) == 1
        submitter_cls = submitter_cls[0]
        return submitter_cls(**self.kwargs)

    @classmethod
    def from_dict(cls, d):
        """
        Create a submitter instance from a dictionary.

        Args:
            d (dict): Dictionary containing submitter configuration.

        Returns:
            Submitter: Configured submitter instance.
        """
        return cls(**d)

    @property
    def submit_folder(self):
        """
        Get the submission folder for the job.

        Returns:
            str: Path to the job submission folder.
        """
        return self.job.folder

    @property
    def submit_script(self):
        """
        Get the submission script filename.

        Returns:
            str: Filename for the job submission script.
        """
        if self.job.label is not None:
            return f"chemsmart_sub_{self.job.label}.sh"
        return "chemsmart_sub.sh"

    @property
    def run_script(self):
        """
        Get the run script filename.

        Returns:
            str: Filename for the job execution script.
        """
        if self.job.label is not None:
            return f"chemsmart_run_{self.job.label}.py"
        return "chemsmart_run.py"

    @property
    def executable(self):
        """
        Get the executable configuration for the job's program.

        Returns:
            Executable: Instance of the appropriate executable handler
            (GaussianExecutable, ORCAExecutable, or NCIPLOTExecutable)
            based on `job.PROGRAM`.

        Raises:
            ValueError: If the job's program is not supported.
        """
        if self.job.PROGRAM.lower() == "gaussian":
            executable = GaussianExecutable.from_servername(self.server.name)
        elif self.job.PROGRAM.lower() == "orca":
            executable = ORCAExecutable.from_servername(self.server.name)
        elif self.job.PROGRAM.lower() == "nciplot":
            executable = NCIPLOTExecutable.from_servername(self.server.name)
        elif self.job.PROGRAM.lower() == "amber":
            executable = AmberExecutable.from_servername(self.server.name)
        else:
            # Need to add programs here to be supported for other types of programs
            raise ValueError(f"Program {self.job.PROGRAM} not supported.")
        return executable

    def write(self, cli_args):
        """
        Write the submission and run scripts for the job.

        Creates both the scheduler-specific submission script and the Python
        run script that will execute the computational job.

        Args:
            cli_args: Command line arguments for job execution.
        """
        if self.job.is_complete():
            logger.warning("Submitting an already complete job.")
        self._write_runscript(cli_args)
        self._write_submitscript()

    def _write_runscript(self, cli_args):
        """
        Write the Python run script for job execution.

        Creates a Python script that handles the actual job execution
        with proper environment setup and argument passing.

        Args:
            cli_args: Command line arguments for the job.
        """
        runscript = RunScript(self.run_script, cli_args)
        logger.debug(f"Writing run script to: {runscript.filename}")
        runscript.write()

    def _write_submitscript(self):
        """
        Write the scheduler submission script.

        Creates a shell script with scheduler directives and job execution
        commands appropriate for the target cluster management system.
        """
        with open(self.submit_script, "w") as f:
            logger.debug(f"Writing submission script to: {self.submit_script}")
            self._write_bash_header(f)
            self._write_scheduler_options(f)
            self._write_program_specifics(f)
            self._write_extra_commands(f)
            self._write_change_to_job_directory(f)
            self._write_job_command(f)

    @staticmethod
    def _write_bash_header(f):
        """
        Write the bash shebang header to the script.

        Args:
            f: File handle for writing the script.
        """
        f.write("#!/bin/bash\n\n")

    @abstractmethod
    def _write_scheduler_options(self, f):
        """
        Write scheduler-specific options to the submission script.

        This method must be implemented by subclasses to provide
        scheduler-specific directives and resource requests.

        Args:
            f: File handle for writing scheduler options.

        Raises:
            NotImplementedError: If not implemented by subclass.
        """
        raise NotImplementedError

    def _write_program_specifics(self, f):
        """
        Write program-specific environment setup to the script.

        Includes conda environment activation, module loading, script
        sourcing, and environment variable configuration specific to
        the computational program being used.

        Args:
            f: File handle for writing program-specific setup.
        """
        self._write_program_specific_conda_env(f)
        self._write_load_program_specific_modules(f)
        self._write_source_program_specific_script(f)
        self._write_program_specific_environment_variables(f)

    def _write_program_specific_conda_env(self, f):
        """
        Write conda environment activation commands.

        Different computational programs may require different conda
        environments for proper execution. This method writes the
        necessary activation commands.

        Args:
            f: File handle for writing conda environment setup.
        """
        if self.executable.conda_env is not None:
            logger.debug(
                f"Writing conda environment: {self.executable.conda_env}"
            )
            f.write("# conda environment\n")
            for line in self.executable.conda_env:
                f.write(line)
            f.write("\n")

    def _write_load_program_specific_modules(self, f):
        """
        Write module loading commands for program dependencies.

        Different computational programs may require loading different
        environment modules for proper execution.

        Args:
            f: File handle for writing module loading commands.
        """
        if self.executable.modules is not None:
            logger.debug(f"Writing modules: {self.executable.modules}")
            f.write("# modules\n")
            for line in self.executable.modules:
                f.write(line)
            f.write("\n")

    def _write_source_program_specific_script(self, f):
        """
        Write script sourcing commands for program setup.

        Different computational programs may require sourcing specific
        setup scripts for proper environment configuration.

        Args:
            f: File handle for writing script sourcing commands.
        """
        if self.executable.scripts is not None:
            logger.debug(f"Writing scripts: {self.executable.scripts}")
            f.write("# program specific scripts\n")
            for line in self.executable.scripts:
                f.write(line)
            f.write("\n")

    def _write_extra_commands(self, f):
        """
        Write additional server-specific commands.

        Extra commands that may be required for the job execution.
        These commands are needed for all jobs across all programs
        and are specific to the server configuration.

        Args:
            f: File handle for writing extra commands.
        """
        if self.server.extra_commands is not None:
            for line in self.server.extra_commands:
                f.write(line)
            f.write("\n")

    def _write_program_specific_environment_variables(self, f):
        """
        Write program-specific environment variables.

        Different computational programs may require different environment
        variables for proper execution. May need to configure different
        scratch folders for different programs (e.g., Gaussian vs ORCA).

        Args:
            f: File handle for writing environment variable exports.
        """
        if self.executable.envars is not None:
            f.write("# Writing program specific environment variables\n")
            for key, value in self.executable.env.items():
                f.write(f"export {key}={value}\n")
            f.write("\n")

    @abstractmethod
    def _write_change_to_job_directory(self, f):
        """
        Write scheduler-specific directory change command.

        Each scheduler system has different environment variables
        for the job submission directory. This method must be
        implemented by subclasses.

        Args:
            f: File handle for writing directory change command.

        Raises:
            NotImplementedError: If not implemented by subclass.
        """
        raise NotImplementedError

    def _write_job_command(self, f):
        """
        Write the final job execution commands.

        Makes the run script executable and executes it in the background,
        then waits for completion.

        Args:
            f: File handle for writing job execution commands.
        """
        f.write(f"chmod +x ./{self.run_script}\n")
        f.write(f"./{self.run_script} &\n")
        f.write("wait\n")

    @classmethod
    def from_scheduler_type(cls, scheduler_type, **kwargs):
        """
        Create a submitter instance for the specified scheduler type.

        Factory method that finds and instantiates the appropriate
        submitter subclass based on the scheduler type name.

        Args:
            scheduler_type (str): Name of the scheduler system
                (e.g., "PBS", "SLURM", "SLF", "FUGAKU").
            **kwargs: Additional arguments passed to the submitter constructor.

        Returns:
            Submitter: Instance of the appropriate submitter subclass
            (one of PBSSubmitter, SLURMSubmitter, SLFSubmitter, or
            FUGAKUSubmitter) configured with the provided kwargs
            (e.g., job and server).

        Raises:
            ValueError: If no submitter is found for the specified scheduler type.
        """
        submitters = cls.subclasses()
        for submitter in submitters:
            if submitter.NAME == scheduler_type:
                return submitter(**kwargs)
        raise ValueError(
            f"Could not find any submitters for scheduler type: {scheduler_type}."
        )


class PBSSubmitter(Submitter):
    """
    PBS (Portable Batch System) job submitter.

    Handles job submission to PBS/Torque cluster management systems.
    Creates PBS-specific submission scripts with appropriate resource
    requests and scheduler directives.

    Attributes:
        NAME (str): Identifier for PBS scheduler type ('PBS').
        name (str): Inherited; instance identifier (often 'PBS').
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission parameters passed to the base class.
    """

    NAME = "PBS"

    def __init__(self, name="PBS", job=None, server=None, **kwargs):
        """
        Initialize PBS submitter.

        Args:
            name (str): Name identifier for this submitter. Defaults to "PBS".
            job: Job instance to be submitted.
            server: Server configuration for PBS submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write PBS-specific scheduler directives.

        Writes PBS directives for output files, resource requests,
        queue selection, walltime, and user notification settings.

        Args:
            f: File handle for writing PBS directives.
        """
        f.write(f"#PBS -o {self.job.label}.pbsout\n")
        f.write(f"#PBS -e {self.job.label}.pbserr\n")
        if self.server.num_gpus > 0:
            f.write(f"#PBS -l gpus={self.server.num_gpus}\n")
        f.write(
            f"#PBS -l select=1:ncpus={self.server.num_cores}:"
            f"mpiprocs={self.server.num_cores}:mem={self.server.mem_gb}G\n"
        )
        # using only one node here
        if self.server.queue_name:
            f.write(f"#PBS -q {self.server.queue_name}\n")
        if self.server.num_hours:
            f.write(f"#PBS -l walltime={self.server.num_hours}:00:00\n")
        if user_settings is not None:
            if user_settings.data.get("PROJECT"):
                f.write(f"#PBS -P {user_settings.data['PROJECT']}\n")
            if user_settings.data.get("EMAIL"):
                f.write(f"#PBS -M {user_settings.data['EMAIL']}\n")
                f.write("#PBS -m abe\n")
        f.write("\n")
        f.write("\n")

    def _write_change_to_job_directory(self, f):
        """
        Write PBS-specific directory change command.

        Uses PBS_O_WORKDIR environment variable to change to the
        job submission directory.

        Args:
            f: File handle for writing directory change command.
        """
        f.write("cd $PBS_O_WORKDIR\n\n")


class SLURMSubmitter(Submitter):
    """
    SLURM (Simple Linux Utility for Resource Management) job submitter.

    Handles job submission to SLURM cluster management systems.
    Creates SLURM-specific submission scripts with appropriate resource
    requests and scheduler directives.

    Attributes:
        NAME (str): Identifier for SLURM scheduler type ('SLURM').
        name (str): Inherited; instance identifier (often 'SLURM').
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission parameters passed to the base class.
    """

    NAME = "SLURM"

    def __init__(self, name="SLURM", job=None, server=None, **kwargs):
        """
        Initialize SLURM submitter.

        Args:
            name (str): Name identifier for this submitter. Defaults to "SLURM".
            job: Job instance to be submitted.
            server: Server configuration for SLURM submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write SLURM-specific scheduler directives.

        Writes SLURM directives for job name, output files, resource
        requests, partition selection, time limits, and user notifications.

        Args:
            f: File handle for writing SLURM directives.
        """
        f.write(f"#SBATCH --job-name={self.job.label}\n")
        f.write(f"#SBATCH --output={self.job.label}.slurmout\n")
        f.write(f"#SBATCH --error={self.job.label}.slurmerr\n")
        if self.server.num_gpus:
            f.write(f"#SBATCH --gres=gpu:{self.server.num_gpus}\n")
        f.write(
            f"#SBATCH --nodes=1 --ntasks-per-node={self.server.num_cores} --mem={self.server.mem_gb}G\n"
        )
        if self.server.queue_name:
            f.write(f"#SBATCH --partition={self.server.queue_name}\n")
        if self.server.num_hours:
            f.write(f"#SBATCH --time={self.server.num_hours}:00:00\n")
        if user_settings is not None:
            if user_settings.data.get("PROJECT"):
                f.write(f"#SBATCH --account={user_settings.data['PROJECT']}\n")
            if user_settings.data.get("EMAIL"):
                f.write(f"#SBATCH --mail-user={user_settings.data['EMAIL']}\n")
                f.write("#SBATCH --mail-type=END,FAIL\n")
        f.write("\n")
        f.write("\n")

    def _write_change_to_job_directory(self, f):
        """
        Write SLURM-specific directory change command.

        Uses SLURM_SUBMIT_DIR environment variable to change to the
        job submission directory.

        Args:
            f: File handle for writing directory change command.
        """
        f.write("cd $SLURM_SUBMIT_DIR\n\n")


class SLFSubmitter(Submitter):
    """
    LSF (Load Sharing Facility) job submitter.

    Handles job submission to IBM LSF cluster management systems.
    Creates LSF-specific submission scripts with appropriate resource
    requests and scheduler directives.

    Note: The class name 'SLFSubmitter' appears to be a typo for 'LSFSubmitter'
    but is maintained for compatibility.

    Attributes:
        NAME (str): Identifier for LSF scheduler type ('SLF').
        name (str): Inherited; instance identifier (often 'SLF').
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission parameters passed to the base class.
    """

    NAME = "SLF"

    def __init__(self, name="SLF", job=None, server=None, **kwargs):
        """
        Initialize LSF submitter.

        Args:
            name (str): Name identifier for this submitter. Defaults to "SLF".
            job: Job instance to be submitted.
            server: Server configuration for LSF submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write LSF-specific scheduler directives.

        Writes LSF directives for job name, output files, project assignment,
        node requests, GPU allocation, and walltime limits.

        Args:
            f: File handle for writing LSF directives.
        """
        f.write(f"#BSUB -J {self.job.label}\n")
        f.write(f"#BSUB -o {self.job.label}.bsubout\n")
        f.write(f"#BSUB -e {self.job.label}.bsuberr\n")
        if user_settings is not None:
            project_number = user_settings.data.get("PROJECT")
        if project_number is not None:
            f.write(f"#BSUB -P {project_number}\n")
        f.write(f"#BSUB -nnodes {self.server.num_nodes}\n")
        if self.server.num_gpus:
            f.write(f"#BSUB -gpu num={self.server.num_gpus}\n")
        f.write(f"#BSUB -W {self.server.num_hours}\n")
        f.write("#BSUB -alloc_flags gpumps\n")
        f.write("\n")
        f.write("\n")

    def _write_change_to_job_directory(self, f):
        """
        Write LSF-specific directory change command.

        Uses LS_SUBCWD environment variable to change to the
        job submission directory.

        Args:
            f: File handle for writing directory change command.
        """
        f.write("cd $LS_SUBCWD\n\n")


class FUGAKUSubmitter(Submitter):
    """
    FUGAKU supercomputer job submitter.

    Handles job submission to the FUGAKU supercomputer system using
    the Fujitsu Job Operation and Management (PJM) scheduler.
    Creates PJM-specific submission scripts with appropriate resource
    requests and scheduler directives.

    Attributes:
        NAME (str): Identifier for FUGAKU scheduler type ('FUGAKU').
        name (str): Inherited; instance identifier (often 'FUGAKU').
        job (Job): Job instance to be submitted.
        server (Server): Server configuration used for submission.
        kwargs (dict): Additional submission parameters passed to the base class.
    """

    NAME = "FUGAKU"

    def __init__(self, name="FUGAKU", job=None, server=None, **kwargs):
        """
        Initialize FUGAKU submitter.

        Args:
            name (str): Name identifier for this submitter. Defaults to "FUGAKU".
            job: Job instance to be submitted.
            server: Server configuration for FUGAKU submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write FUGAKU PJM-specific scheduler directives.

        Writes PJM directives for resource group, node allocation,
        elapsed time, MPI processes, project assignment, and output files.
        Includes FUGAKU-specific optimizations like LLIO cache settings.

        Args:
            f: File handle for writing PJM directives.
        """
        if user_settings is not None:
            f.write(f'#PJM -L rscgrp={user_settings.data["RSCGRP"]}\n')
        f.write("#PJM -L node=1\n")  # using one node here
        f.write(f"#PJM -L elapse={self.server.num_hours}\n")
        f.write(f"#PJM --mpi proc={self.server.num_cores}\n")
        f.write(f"#PJM -g {self.project}\n")
        f.write("#PJM -o pjm.%j.out\n")
        f.write("#PJM -e pjm.%j.err\n")
        f.write("#PJM -x PJM_LLIO_GFSCACHE=/vol0005:/vol0004\n")
        f.write("#PJM -S\n")
        f.write("\n")
        f.write("\n")

    def _write_change_to_job_directory(self, f):
        """
        Write FUGAKU PJM-specific directory change command.

        Uses PJM_O_WORKDIR environment variable to change to the
        job submission directory.

        Args:
            f: File handle for writing directory change command.
        """
        f.write("cd $PJM_O_WORKDIR\n\n")


class AmberMCPBSLURMSubmitter(SLURMSubmitter):
    """
    Specialized SLURM submitter for Amber MCPB workflows.

    This submitter creates sequential MCPB workflow scripts that execute
    multiple steps in order: preparation, small model generation, large model
    generation, QM calculations, and finalization.
    """

    NAME = "AMBER_MCPB_SLURM"

    def __init__(
        self, name="AMBER_MCPB_SLURM", job=None, server=None, **kwargs
    ):
        """
        Initialize Amber MCPB SLURM submitter.

        Args:
            name (str): Name identifier for this submitter.
            job: Amber MCPB job instance to be submitted.
            server: Server configuration for SLURM submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write SLURM directives optimized for MCPB workflows.

        Uses longer time limits and appropriate resource allocation
        for sequential MCPB processing.

        Args:
            f: File handle for writing SLURM directives.
        """
        # Use mcpb_workflow as default job name for MCPB jobs
        job_name = getattr(self.job, "mcpb_job_name", "mcpb_workflow")
        f.write(f"#SBATCH --job-name={job_name}\n")
        f.write(f"#SBATCH --output={self.job.label}.slurmout\n")
        f.write(f"#SBATCH --error={self.job.label}.slurmerr\n")

        # MCPB workflows typically need longer time limits
        time_hours = self.server.num_hours if self.server.num_hours else 48
        f.write(f"#SBATCH --time={time_hours}:00:00\n")

        # Use cpus-per-task instead of ntasks-per-node for MCPB
        f.write(f"#SBATCH --cpus-per-task={self.server.num_cores}\n")
        f.write(f"#SBATCH --mem={self.server.mem_gb}G\n")

        if self.server.queue_name:
            f.write(f"#SBATCH --partition={self.server.queue_name}\n")

        if user_settings is not None:
            if user_settings.data.get("PROJECT"):
                f.write(f"#SBATCH --account={user_settings.data['PROJECT']}\n")
            if user_settings.data.get("EMAIL"):
                f.write(f"#SBATCH --mail-user={user_settings.data['EMAIL']}\n")
                f.write("#SBATCH --mail-type=END,FAIL\n")
        f.write("\n")

    def _write_job_command(self, f):
        """
        Write the sequential MCPB workflow commands instead of the standard run script.

        Creates a sequential workflow that executes MCPB steps in order,
        running QM calculations between steps as needed. The filenames and
        commands are dynamically generated based on the job's settings.

        Args:
            f: File handle for writing job commands.
        """
        settings = self.job.settings

        # Get system-specific names from job settings
        reduce_command = getattr(settings, "reduce_command", None)
        antechamber_mol2_command = getattr(
            settings, "antechamber_mol2_command", None
        )
        rename_antechamber_output_command = getattr(
            settings, "rename_antechamber_output_command", None
        )
        parmchk2_command = getattr(settings, "parmchk2_command", None)
        metalpdb2mol_command = getattr(settings, "metalpdb2mol_command", None)
        ambpdb_command = getattr(settings, "ambpdb_command", None)
        combine_command = getattr(settings, "combine_command", None)
        pdb4amber_command = getattr(settings, "pdb4amber_command", None)
        mcpb_step1_command = getattr(settings, "mcpb_step1_command", None)
        mcpb_step2_command = getattr(settings, "mcpb_step2_command", None)

        group_name = getattr(settings, "group_name", "system")
        gaussian_exe = getattr(settings, "gaussian_exe", "g16")

        # Define QM input and output files based on group_name
        opt_com = f"{group_name}_small_opt.com"
        opt_log = f"{group_name}_small_opt.log"
        opt_chk = f"{group_name}_small_opt.chk"
        opt_fchk = f"{group_name}_small_opt.fchk"

        fc_com = f"{group_name}_small_fc.com"
        fc_log = f"{group_name}_small_fc.log"

        mk_com = f"{group_name}_large_mk.com"
        mk_log = f"{group_name}_large_mk.log"

        f.write("# Amber Force Field Parameterization Sequential Workflow\n")
        f.write("# Prepare PDB and mol2 files for the ligand\n")
        f.write(f"{reduce_command}\n ")
        f.write(f"{antechamber_mol2_command}\n")
        f.write(f"{rename_antechamber_output_command}\n")
        f.write(f"{parmchk2_command}\n")

        f.write("# Prepare PDB and mol2 files for the metal ion\n")
        f.write(f"{metalpdb2mol_command}\n")
        # todo: processing for water molecules
        # f.write("# Prepare PDB and mol2 files for the water\n")

        f.write("# Standard pdb file preparation\n")
        if ambpdb_command:
            f.write(f"{antechamber_mol2_command}\n")
            f.write(f"{combine_command}\n")
        # write pdb4chamber command directly if ambpdb is not used
        f.write(f"{pdb4amber_command}\n")

        f.write(
            "# Amber MCPB Sequential Workflow\n"
            "# Generate the PDB, Gaussian, GAMESS-US and fingerprint modeling files.\n"
        )
        # Step 1: Preparation
        f.write('echo "Step 1: MCPB Preparation"\n')
        f.write(f"{mcpb_step1_command}\n")
        f.write("set -e  # Exit on any error\n\n")

        # Step2 running QM calculations
        # QM calculations
        if getattr(settings, "auto_submit_qm", False):
            f.write(
                'echo "Running Gaussian calculations for small model..."\n'
            )

            # Run optimization
            f.write(f"if [ -f {opt_com} ]; then\n")
            f.write(f"    {gaussian_exe} < {opt_com} > {opt_log}\n")
            f.write(
                f"    if [ -f {opt_chk} ]; then formchk {opt_chk} {opt_fchk}; fi\n"
            )
            f.write("fi\n\n")

            # Run frequency calculation
            f.write(f"if [ -f {fc_com} ]; then\n")
            f.write(f"    {gaussian_exe} < {fc_com} > {fc_log}\n")
            f.write("fi\n\n")

            # Run RESP calculation
            f.write(f"if [ -f {mk_com} ]; then\n")
            f.write(f"    {gaussian_exe} < {mk_com} > {mk_log}\n")
            f.write("fi\n\n")
        else:
            f.write(
                'echo "QM calculations must be run manually before continuing to step 4"\n'
            )
            f.write(f'echo "Files to run: {opt_com}, {fc_com}, {mk_com}"\n\n')

        # Step 3: final modelling
        f.write('echo "Step 2: Small Model Generation"\n')
        f.write(f"{mcpb_step2_command}\n")

        # Step 3: Large model
        # f.write('echo "Step 3: Large Model Generation"\n')
        # f.write(f"MCPB.py -i {mcpb_input} -s 3\n\n")
        #
        # # Step 4: Finalization
        # if getattr(settings, "auto_submit_qm", False):
        #     f.write('echo "Step 4: MCPB Finalization"\n')
        #     f.write(f"MCPB.py -i {mcpb_input} -s 4\n\n")
        #
        #     # Run tleap to generate topology files
        #     tleap_in = f"{group_name}_tleap.in"
        #     tleap_out = f"{group_name}_tleap.out"
        #     f.write('echo "Running tleap to generate topology files..."\n')
        #     f.write(f"if [ -f {tleap_in} ]; then\n")
        #     f.write(f"    tleap -s -f {tleap_in} > {tleap_out}\n")
        #     f.write("fi\n\n")
        #
        #     f.write('echo "MCPB workflow completed successfully!"\n')
        #     f.write(
        #         f'echo "Generated files: {group_name}_solv.prmtop, {group_name}_solv.inpcrd"\n'
        #     )
        #
        # f.write("\n")


class AmberMCPBPBSSubmitter(PBSSubmitter):
    """
    Specialized PBS submitter for Amber MCPB workflows.

    Similar to AmberMCPBSLURMSubmitter but for PBS scheduler systems.
    """

    NAME = "AMBER_MCPB_PBS"

    def __init__(self, name="AMBER_MCPB_PBS", job=None, server=None, **kwargs):
        """
        Initialize Amber MCPB PBS submitter.

        Args:
            name (str): Name identifier for this submitter.
            job: Amber MCPB job instance to be submitted.
            server: Server configuration for PBS submission.
            **kwargs: Additional submission parameters.
        """
        super().__init__(name=name, job=job, server=server, **kwargs)

    def _write_scheduler_options(self, f):
        """
        Write PBS directives optimized for MCPB workflows.

        Args:
            f: File handle for writing PBS directives.
        """
        job_name = getattr(self.job, "mcpb_job_name", "mcpb_workflow")
        f.write(f"#PBS -N {job_name}\n")
        f.write(f"#PBS -o {self.job.label}.pbsout\n")
        f.write(f"#PBS -e {self.job.label}.pbserr\n")

        # MCPB workflows typically need longer time limits
        time_hours = self.server.num_hours if self.server.num_hours else 48
        f.write(f"#PBS -l walltime={time_hours}:00:00\n")

        f.write(
            f"#PBS -l select=1:ncpus={self.server.num_cores}:mem={self.server.mem_gb}GB\n"
        )

        if self.server.queue_name:
            f.write(f"#PBS -q {self.server.queue_name}\n")

        if user_settings is not None:
            if user_settings.data.get("PROJECT"):
                f.write(f"#PBS -P {user_settings.data['PROJECT']}\n")
            if user_settings.data.get("EMAIL"):
                f.write(f"#PBS -M {user_settings.data['EMAIL']}\n")
                f.write("#PBS -m abe\n")
        f.write("\n")

    def _write_job_command(self, f):
        """
        Write the sequential MCPB workflow commands for PBS.

        Uses the same workflow logic as the SLURM version, generating
        dynamic, system-specific filenames for all commands.

        Args:
            f: File handle for writing job commands.
        """
        # This method now directly implements the workflow logic
        # instead of calling the SLURM version's method.
        settings = self.job.settings

        mcpb_input = getattr(
            settings,
            "mcpb_input",
            f"{getattr(settings, 'group_name', 'system')}.in",
        )
        group_name = getattr(settings, "group_name", "system")
        gaussian_exe = getattr(settings, "gaussian_exe", "g16")

        opt_com = f"{group_name}_small_opt.com"
        opt_log = f"{group_name}_small_opt.log"
        opt_chk = f"{group_name}_small_opt.chk"
        opt_fchk = f"{group_name}_small_opt.fchk"

        fc_com = f"{group_name}_small_fc.com"
        fc_log = f"{group_name}_small_fc.log"
        fc_chk = f"{group_name}_small_fc.chk"
        fc_fchk = f"{group_name}_small_fc.fchk"

        mk_com = f"{group_name}_small_mk.com"
        mk_log = f"{group_name}_small_mk.log"
        mk_chk = f"{group_name}_small_mk.chk"
        mk_fchk = f"{group_name}_small_mk.fchk"

        f.write("# Amber MCPB Sequential Workflow (PBS)\n")
        f.write("set -e  # Exit on any error\n\n")

        f.write('echo "Step 1: MCPB Preparation"\n')
        f.write(f"MCPB.py -i {mcpb_input} -s 1\n\n")

        f.write('echo "Step 2: Small Model Generation"\n')
        f.write(f"MCPB.py -i {mcpb_input} -s 2\n\n")

        f.write('echo "Step 3: Large Model Generation"\n')
        f.write(f"MCPB.py -i {mcpb_input} -s 3\n\n")

        if getattr(settings, "auto_submit_qm", False):
            f.write(
                'echo "Running Gaussian calculations for small model..."\n'
            )

            f.write(f"if [ -f {opt_com} ]; then\n")
            f.write(f"    {gaussian_exe} < {opt_com} > {opt_log}\n")
            f.write(
                f"    if [ -f {opt_chk} ]; then formchk {opt_chk} {opt_fchk}; fi\n"
            )
            f.write("fi\n\n")

            f.write(f"if [ -f {fc_com} ]; then\n")
            f.write(f"    {gaussian_exe} < {fc_com} > {fc_log}\n")
            f.write(
                f"    if [ -f {fc_chk} ]; then formchk {fc_chk} {fc_fchk}; fi\n"
            )
            f.write("fi\n\n")

            f.write(f"if [ -f {mk_com} ]; then\n")
            f.write(f"    {gaussian_exe} < {mk_com} > {mk_log}\n")
            f.write(
                f"    if [ -f {mk_chk} ]; then formchk {mk_chk} {mk_fchk}; fi\n"
            )
            f.write("fi\n\n")
        else:
            f.write(
                'echo "QM calculations must be run manually before continuing to step 4"\n'
            )
            f.write(f'echo "Files to run: {opt_com}, {fc_com}, {mk_com}"\n\n')

        if getattr(settings, "auto_submit_qm", False):
            f.write('echo "Step 4: MCPB Finalization"\n')
            f.write(f"MCPB.py -i {mcpb_input} -s 4\n\n")

            tleap_in = f"{group_name}_tleap.in"
            tleap_out = f"{group_name}_tleap.out"
            f.write('echo "Running tleap to generate topology files..."\n')
            f.write(f"if [ -f {tleap_in} ]; then\n")
            f.write(f"    tleap -s -f {tleap_in} > {tleap_out}\n")
            f.write("fi\n\n")

            f.write('echo "MCPB workflow completed successfully!"\n')
            f.write(
                f'echo "Generated files: {group_name}_solv.prmtop, {group_name}_solv.inpcrd"\n'
            )

        f.write("\n")
