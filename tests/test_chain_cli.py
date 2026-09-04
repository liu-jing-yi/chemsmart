from unittest.mock import MagicMock, patch

import pytest
from click.testing import CliRunner

from chemsmart.cli.chain.chain import chain
from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.run import run
from chemsmart.cli.sub import sub
from chemsmart.jobs.chain import ChainJob
from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.user import CHEMSMARTUserSettings


@pytest.fixture()
def isolated_config_dir(tmp_path, monkeypatch):
    config_dir = tmp_path / "chemsmart_config"
    config_dir.mkdir()
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_dir))
    monkeypatch.setattr(
        CHEMSMARTUserSettings,
        "USER_CONFIG_DIR",
        str(config_dir),
    )
    return config_dir


def _write_chain_yaml(isolated_config_dir, name, content):
    chain_dir = isolated_config_dir / "chain"
    chain_dir.mkdir(exist_ok=True)
    (chain_dir / f"{name}.yaml").write_text(content)
    return name


def _invoke_chain(args, obj=None, **extra):
    runner = CliRunner()
    if obj is None:
        obj = {"jobrunner": MagicMock()}
    return runner.invoke(chain, args, obj=obj, catch_exceptions=False, **extra)


class TestChainHelp:
    def test_run_help_lists_chain(self):
        result = CliRunner().invoke(run, ["--help"])
        assert result.exit_code == 0, result.output
        assert "chain" in result.output

    def test_sub_help_lists_chain(self):
        result = CliRunner().invoke(sub, ["--help"])
        assert result.exit_code == 0, result.output
        assert "chain" in result.output

    def test_chain_help_lists_custom_subcommand(self):
        result = CliRunner().invoke(run, ["chain", "--help"])
        assert result.exit_code == 0, result.output
        assert "\n  custom" in result.output
        assert "\n  pka" not in result.output
        assert "\n  fukui" not in result.output
        assert "\n  reaction" not in result.output
        assert "\n  gaussian" not in result.output
        assert "\n  orca" not in result.output
        assert "\n  xtb" not in result.output
        assert "\n  crest" not in result.output


class TestChainPipeline:
    def test_pipeline_builds_chain_job(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            """
crest: test
xtb: test
gaussian: gas_solv
orca: gas_solv
""",
        )
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-s",
                "crest:conformers",
                "-s",
                "xtb:opt",
                "-s",
                "gaussian:opt",
                "-s",
                "orca:sp",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        job = result.return_value
        assert isinstance(job, ChainJob)
        assert job.label == "crest_best"
        assert [phase.name for phase in job.phases] == [
            "00_crest_conformers",
            "01_xtb_opt",
            "02_gaussian_opt",
            "03_orca_sp",
        ]
        children = job.phases[0].resolve_jobs()
        assert len(children) == 1
        assert isinstance(children[0], CRESTConformerSearchJob)
        assert children[0].label == "crest_best_00_crest_conformers"

    def test_pipeline_without_steps_is_usage_error(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "requires at least one -s/--steps" in result.output

    def test_pipeline_without_filename_is_usage_error(
        self, isolated_config_dir
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            ["-p", "combined", "-s", "gaussian:opt"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "requires a structure file" in result.output

    def test_unknown_nested_command_is_usage_error(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "pka",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "No such command 'pka'" in result.output

    def test_invalid_step_token_is_usage_error(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-s",
                "gaussian",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "expected PROGRAM:JOB" in result.output

    def test_nested_only_step_is_usage_error(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = CliRunner().invoke(
            chain,
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-s",
                "gaussian:pka",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "do not support gaussian 'pka'" in result.output


class TestChainPipelineExtras:
    def test_pipeline_quoted_step_applies_opt_options(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-s",
                "gaussian: -o maxstep=8,maxsize=12 opt",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        job = result.return_value
        assert isinstance(job, ChainJob)
        child = job.phases[0].resolve_jobs()[0]
        assert (
            child.settings.additional_opt_options_in_route
            == "maxstep=8,maxsize=12"
        )

    def test_custom_subcommand_builds_pipeline(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        result = _invoke_chain(
            [
                "-p",
                "combined",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-s",
                "gaussian:opt",
                "custom",
            ],
            standalone_mode=False,
        )
        assert result.exit_code == 0, result.output
        assert isinstance(result.return_value, ChainJob)
        assert [phase.name for phase in result.return_value.phases] == [
            "00_gaussian_opt"
        ]

    def test_missing_chain_project_is_usage_error(self, isolated_config_dir):
        result = CliRunner().invoke(
            chain,
            ["-p", "missing", "custom"],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output

    def test_run_gaussian_without_chain_uses_own_project(
        self, single_molecule_xyz_file
    ):
        args = [
            "-p",
            "test",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "opt",
        ]
        with (
            patch(
                "chemsmart.settings.gaussian.GaussianProjectSettings."
                "from_project",
                wraps=GaussianProjectSettings.from_project,
            ) as gaussian_from_project,
            patch(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
                return_value=MagicMock(),
            ),
        ):
            result = CliRunner().invoke(
                gaussian,
                args,
                obj={"jobrunner": MagicMock()},
                catch_exceptions=False,
            )
        assert result.exit_code == 0, result.output
        gaussian_from_project.assert_called_once_with("test")
