from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import click
import pytest
from click.testing import CliRunner

from chemsmart.cli.crest.crest import crest
from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.orca.orca import orca
from chemsmart.cli.utils import (
    CHAIN_PROJECT_SETTINGS_KEY,
    resolve_program_project,
)
from chemsmart.cli.xtb.xtb import xtb
from chemsmart.settings.chain import ChainProjectSettings
from chemsmart.settings.crest import CRESTProjectSettings
from chemsmart.settings.gaussian import GaussianProjectSettings
from chemsmart.settings.orca import ORCAProjectSettings
from chemsmart.settings.xtb import XTBProjectSettings


def _ctx(obj):
    return SimpleNamespace(obj=obj)


def _chain_settings(**aliases):
    return ChainProjectSettings(
        aliases=aliases,
        steps=(),
        project_name="combined",
    )


class TestResolveProgramProject:
    def test_explicit_project_wins_over_chain_alias(self):
        settings = _chain_settings(gaussian="gas_solv")
        ctx = _ctx({CHAIN_PROJECT_SETTINGS_KEY: settings})
        assert (
            resolve_program_project(ctx, "gaussian", "defaults") == "defaults"
        )

    def test_uses_chain_alias_when_project_omitted(self):
        settings = _chain_settings(gaussian="gas_solv")
        ctx = _ctx({CHAIN_PROJECT_SETTINGS_KEY: settings})
        assert resolve_program_project(ctx, "gaussian", None) == "gas_solv"

    def test_missing_program_section_is_usage_error(self):
        settings = _chain_settings(crest="test")
        ctx = _ctx({CHAIN_PROJECT_SETTINGS_KEY: settings})
        with pytest.raises(click.UsageError, match="has no gaussian section"):
            resolve_program_project(ctx, "gaussian", None)

    def test_without_chain_settings_returns_project_unchanged(self):
        ctx = _ctx({})
        assert (
            resolve_program_project(ctx, "gaussian", "gas_solv") == "gas_solv"
        )
        assert resolve_program_project(ctx, "gaussian", None) is None

    def test_none_ctx_obj_returns_project_unchanged(self):
        ctx = _ctx(None)
        assert resolve_program_project(ctx, "gaussian", None) is None
        assert resolve_program_project(ctx, "gaussian", "test") == "test"


def _cli_obj(chain_settings=None):
    obj = {"jobrunner": MagicMock()}
    if chain_settings is not None:
        obj[CHAIN_PROJECT_SETTINGS_KEY] = chain_settings
    return obj


def _invoke(cli, args, obj, job_class_path):
    runner = CliRunner()
    with patch(job_class_path, return_value=MagicMock()):
        return runner.invoke(cli, args, obj=obj, catch_exceptions=False)


class TestProgramCLIProjectInject:
    def test_gaussian_uses_chain_alias_and_skips_other_programs(
        self, single_molecule_xyz_file
    ):
        settings = _chain_settings(gaussian="gas_solv", orca="gas_solv")
        args = [
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
        ):
            result = _invoke(
                gaussian,
                args,
                _cli_obj(settings),
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            )
        assert result.exit_code == 0, result.output
        gaussian_from_project.assert_called_once_with("gas_solv")
        orca_from_project.assert_not_called()
        xtb_from_project.assert_not_called()
        crest_from_project.assert_not_called()

    def test_gaussian_explicit_p_overrides_chain_alias(
        self, single_molecule_xyz_file
    ):
        settings = _chain_settings(gaussian="defaults")
        args = [
            "-p",
            "gas_solv",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "opt",
        ]
        with patch(
            "chemsmart.settings.gaussian.GaussianProjectSettings."
            "from_project",
            wraps=GaussianProjectSettings.from_project,
        ) as gaussian_from_project:
            result = _invoke(
                gaussian,
                args,
                _cli_obj(settings),
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            )
        assert result.exit_code == 0, result.output
        gaussian_from_project.assert_called_once_with("gas_solv")

    def test_gaussian_without_chain_settings_uses_own_p(
        self, single_molecule_xyz_file
    ):
        args = [
            "-p",
            "gas_solv",
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "opt",
        ]
        with patch(
            "chemsmart.settings.gaussian.GaussianProjectSettings."
            "from_project",
            wraps=GaussianProjectSettings.from_project,
        ) as gaussian_from_project:
            result = _invoke(
                gaussian,
                args,
                _cli_obj(),
                "chemsmart.jobs.gaussian.opt.GaussianOptJob",
            )
        assert result.exit_code == 0, result.output
        gaussian_from_project.assert_called_once_with("gas_solv")

    def test_gaussian_missing_chain_section_is_usage_error(
        self, single_molecule_xyz_file
    ):
        settings = _chain_settings(crest="test")
        args = [
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "opt",
        ]
        result = _invoke(
            gaussian,
            args,
            _cli_obj(settings),
            "chemsmart.jobs.gaussian.opt.GaussianOptJob",
        )
        assert result.exit_code == 2, result.output
        assert "has no gaussian section" in result.output

    def test_orca_uses_chain_alias(self, single_molecule_xyz_file):
        settings = _chain_settings(orca="gas_solv")
        args = [
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "sp",
        ]
        with patch(
            "chemsmart.settings.orca.ORCAProjectSettings.from_project",
            wraps=ORCAProjectSettings.from_project,
        ) as orca_from_project:
            result = _invoke(
                orca,
                args,
                _cli_obj(settings),
                "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            )
        assert result.exit_code == 0, result.output
        orca_from_project.assert_called_once_with("gas_solv")

    def test_xtb_uses_chain_alias(self, single_molecule_xyz_file):
        settings = _chain_settings(xtb="test")
        args = [
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "opt",
        ]
        with patch(
            "chemsmart.settings.xtb.XTBProjectSettings.from_project",
            wraps=XTBProjectSettings.from_project,
        ) as xtb_from_project:
            result = _invoke(
                xtb,
                args,
                _cli_obj(settings),
                "chemsmart.jobs.xtb.opt.XTBOptJob",
            )
        assert result.exit_code == 0, result.output
        xtb_from_project.assert_called_once_with("test")

    def test_crest_uses_chain_alias(self, single_molecule_xyz_file):
        settings = _chain_settings(crest="test")
        args = [
            "-f",
            single_molecule_xyz_file,
            "-c",
            "0",
            "-m",
            "1",
            "conformers",
        ]
        with patch(
            "chemsmart.settings.crest.CRESTProjectSettings.from_project",
            wraps=CRESTProjectSettings.from_project,
        ) as crest_from_project:
            result = _invoke(
                crest,
                args,
                _cli_obj(settings),
                "chemsmart.jobs.crest.conformers.CRESTConformerSearchJob",
            )
        assert result.exit_code == 0, result.output
        crest_from_project.assert_called_once_with("test")
