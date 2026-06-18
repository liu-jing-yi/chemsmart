import os
from io import StringIO

from chemsmart.settings.executable import GaussianExecutable, ORCAExecutable
from chemsmart.settings.server import Server
from chemsmart.settings.submitters import (
    PBSSubmitter,
    SLURMSubmitter,
    render_run_script_contents,
)


def _pbs_submitter(label="job1", program="orca"):
    server = Server(
        "custom-pbs",
        SCHEDULER="PBS",
        NUM_CORES=8,
        MEM_GB=24,
        NUM_GPUS=0,
    )
    job = type(
        "DummyJob",
        (),
        {"label": label, "folder": ".", "PROGRAM": program},
    )()
    return PBSSubmitter(job=job, server=server)


class TestServer:
    def test_server_yaml(self, server_yaml_file):
        assert os.path.exists(server_yaml_file)
        assert os.path.isfile(server_yaml_file)
        server = Server.from_yaml(name=server_yaml_file)
        assert server.scheduler.lower() == "pbs"
        assert server.queue_name == "normal"
        assert server.num_hours == 24
        assert server.mem_gb == 375
        assert server.num_cores == 64
        assert server.num_gpus == 0
        assert server.num_threads == 64
        assert server.submit_command == "qsub"
        assert server.scratch_dir is None
        assert server.use_hosts is True
        assert (
            server.extra_commands == """export PATH=$HOME/bin/chemsmart:$PATH
export PATH=$HOME/bin/chemsmart/chemsmart/cli:$PATH
export PATH=$HOME/bin/chemsmart/chemsmart/scripts:$PATH
export PYTHONPATH=$HOME/bin/chemsmart:$PYTHONPATH
"""
        )
        assert server.extra_scheduler_directives == "#PBS -m abe\n"

    def test_gaussian_executable(self, server_yaml_file):
        gaussian_executable = GaussianExecutable.from_servername(
            server_yaml_file
        )
        assert gaussian_executable.executable_folder == os.path.expanduser(
            "~/programs/g16"
        )
        assert gaussian_executable.local_run is True

        gaussian_conda_env = """source ~/anaconda3/etc/profile.d/conda.sh
conda activate ~/anaconda3/envs/chemsmart
"""
        assert gaussian_executable.conda_env == gaussian_conda_env

        gaussian_modules = """module purge
module load craype-x86-rome
module load libfabric/1.11.0.4.125
"""
        assert gaussian_executable.modules == gaussian_modules

        assert (
            gaussian_executable.scripts
            == 'tcsh -c "source ~/programs/g16/bsd/g16.login"\n'
        )

        gassian_envars = """export SCRATCH=~/scratch
export GAUSS_EXEDIR=~/programs/g16
export g16root=~/programs/g16

"""
        assert gaussian_executable.envars == gassian_envars

    def test_orca_executable(self, server_yaml_file):
        orca_executable = ORCAExecutable.from_servername(server_yaml_file)
        assert orca_executable.executable_folder == os.path.expanduser(
            "~/programs/orca_6_0_0"
        )
        assert orca_executable.local_run is False

        assert orca_executable.conda_env is None

        assert orca_executable.modules is None

        assert orca_executable.scripts is None

        orca_envars = """export PATH=~/programs/openmpi-4.1.6/build/bin:$PATH
export LD_LIBRARY_PATH=~/programs/openmpi-4.1.6/build/lib:$LD_LIBRARY_PATH
"""
        assert orca_executable.envars == orca_envars

    def test_slurm_submitter_writes_extra_scheduler_directives(self):
        server = Server(
            "custom-slurm",
            SCHEDULER="SLURM",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_GPUS=0,
            EXTRA_SCHEDULER_DIRECTIVES="#SBATCH --reservation=xlzhang_1\n",
        )
        job = type("DummyJob", (), {"label": "job1"})()
        submitter = SLURMSubmitter(job=job, server=server)

        buffer = StringIO()
        submitter._write_scheduler_options(buffer)
        assert "#SBATCH --reservation=xlzhang_1\n" in buffer.getvalue()

    def test_pbs_submitter_writes_extra_scheduler_directives(self):
        server = Server(
            "custom-pbs",
            SCHEDULER="PBS",
            NUM_CORES=8,
            MEM_GB=24,
            NUM_GPUS=0,
            EXTRA_SCHEDULER_DIRECTIVES="#PBS -m abe\n",
        )
        job = type("DummyJob", (), {"label": "job1"})()
        submitter = PBSSubmitter(job=job, server=server)

        buffer = StringIO()
        submitter._write_scheduler_options(buffer)
        assert "#PBS -m abe\n" in buffer.getvalue()


class TestRunScriptGeneration:
    def test_run_script_includes_prejob_hooks(self):
        cli_args = ["orca", "-f", "test.xyz", "sp"]
        contents = render_run_script_contents(cli_args)

        assert "from chemsmart.cli.prejob import run_prejob_hooks" in contents
        assert "CLI_ARGS = ['orca', '-f', 'test.xyz', 'sp']" in contents
        assert "run_prejob_hooks(CLI_ARGS)" in contents
        assert "run(CLI_ARGS)" in contents
        assert contents.index("run_prejob_hooks(CLI_ARGS)") < contents.index(
            "run(CLI_ARGS)"
        )

    def test_submit_script_only_launches_python_run_script(self):
        submitter = _pbs_submitter(label="ionic_crystal_qmmm")
        buffer = StringIO()
        submitter._write_job_command(buffer)
        contents = buffer.getvalue()

        assert "chmod +x ./chemsmart_run_ionic_crystal_qmmm.py\n" in contents
        assert "./chemsmart_run_ionic_crystal_qmmm.py &\n" in contents
        assert "wait\n" in contents
        assert "run_prejob_hooks" not in contents
        assert "ORCA_crystalprep" not in contents
