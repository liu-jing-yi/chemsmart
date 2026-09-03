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

    def test_chain_help_lists_nested_programs(self):
        result = CliRunner().invoke(run, ["chain", "--help"])
        assert result.exit_code == 0, result.output
        assert "gaussian" in result.output
        assert "orca" in result.output
        assert "xtb" in result.output
        assert "crest" in result.output

    def test_chain_gaussian_pka_help(
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
                "gaussian",
                "-f",
                single_molecule_xyz_file,
                "pka",
                "--help",
            ]
        )
        assert result.exit_code == 0, result.output
        assert "submit" in result.output
        assert "batch" in result.output


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

    def test_steps_with_nested_command_is_usage_error(
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
                "gaussian:opt",
                "gaussian",
                "opt",
            ],
            obj={"jobrunner": MagicMock()},
        )
        assert result.exit_code == 2, result.output
        assert "Cannot combine -s/--steps" in result.output

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


class TestChainNestedSlice:
    def test_gaussian_opt_uses_chain_alias_and_skips_other_programs(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\norca: gas_solv\n",
        )
        args = [
            "-p",
            "combined",
            "gaussian",
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
                "chemsmart.settings.orca.ORCAProjectSettings.from_project"
            ) as orca_from_project,
            patch(
                "chemsmart.settings.xtb.XTBProjectSettings.from_project"
            ) as xtb_from_project,
            patch(
                "chemsmart.settings.crest.CRESTProjectSettings.from_project"
            ) as crest_from_project,
            patch(
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
                return_value=MagicMock(),
            ) as mock_job,
        ):
            result = _invoke_chain(args)
        assert result.exit_code == 0, result.output
        gaussian_from_project.assert_called_once_with("gas_solv")
        orca_from_project.assert_not_called()
        xtb_from_project.assert_not_called()
        crest_from_project.assert_not_called()
        assert mock_job.called

    def test_nested_gaussian_uses_chain_level_file_charge_multiplicity(
        self, isolated_config_dir, single_molecule_xyz_file
    ):
        _write_chain_yaml(
            isolated_config_dir,
            "combined",
            "gaussian: gas_solv\n",
        )
        args = [
            "-p",
            "combined",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "gaussian",
            "opt",
        ]
        with patch(
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            return_value=MagicMock(),
        ) as mock_job:
            result = _invoke_chain(args)
        assert result.exit_code == 0, result.output
        assert mock_job.called
        call_kwargs = mock_job.call_args.kwargs
        assert call_kwargs["settings"].charge == 0
        assert call_kwargs["settings"].multiplicity == 1

    def test_missing_chain_project_is_usage_error(self, isolated_config_dir):
        result = CliRunner().invoke(
            chain,
            ["-p", "missing", "gaussian", "opt"],
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
