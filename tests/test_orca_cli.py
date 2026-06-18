"""
Tests for ORCA CLI option propagation and subcommand behaviour.

This module verifies that solvent-related options
 (``-sm``/``--solvent-model``, ``-si``/``--solvent-id``,
 ``-so``/``--solvent-options``, and ``--remove-solvent``) on the
 ``orca`` CLI *group* and on individual subcommands (``sp``, ``opt``,
 ``ts``) are correctly propagated to the merged
 :class:`~chemsmart.jobs.orca.settings.ORCAJobSettings`.

Each test uses :class:`click.testing.CliRunner` to invoke the ``orca``
group and :mod:`unittest.mock` to intercept the job constructor so that
the merged settings can be inspected without running an actual calculation.
"""


class TestORCASolventCLISpCommand:
    """CLI solvent options propagated to the ``sp`` subcommand."""

    def test_solvent_model_and_id_injected_into_sp_settings_group_level(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water`` at group level sets solvent on sp settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cpcm",
                "-si",
                "water",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert (
            settings is not None
        ), "ORCASinglePointJob was never instantiated"
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_solvent_model_and_id_injected_into_sp_settings_subcommand_level(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``sp -sm cpcm -si water`` at subcommand level sets solvent on sp settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "sp",
                "-sm",
                "cpcm",
                "-si",
                "water",
            ],
        )

        assert result.exit_code == 0, result.output
        assert (
            settings is not None
        ), "ORCASinglePointJob was never instantiated"
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_solvent_options_injected_into_sp_settings(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water -so 'Epsilon 78.36'`` sets additional options on sp."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cpcm",
                "-si",
                "water",
                "-so",
                "Epsilon 78.36",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "Epsilon 78.36"

    def test_remove_solvent_clears_solvent_from_sp(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``--remove-solvent`` nulls the solvent on a project that has one."""
        # The ``solv`` project sets solvent_model=smd and solvent_id=cyclohexane
        # for every job type.  ``--remove-solvent`` must strip these from the
        # merged settings.
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_no_solvent_options_leaves_project_sp_settings_unchanged(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """No solvent CLI flags leave the project sp solvent settings intact."""
        # ``gas_solv`` project sp has smd/cyclohexane; no CLI flags → preserved.
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "cyclohexane"


class TestORCASolventCLIOptCommand:
    """CLI solvent options propagated to the ``opt`` subcommand."""

    def test_solvent_model_and_id_injected_into_opt_settings(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water`` sets solvent on the opt job settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cpcm",
                "-si",
                "water",
                "opt",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCAOptJob was never instantiated"
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_solvent_options_injected_into_opt_settings(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water -so 'Epsilon 78.36'`` propagates to opt."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cpcm",
                "-si",
                "water",
                "-so",
                "Epsilon 78.36",
                "opt",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "Epsilon 78.36"

    def test_remove_solvent_clears_solvent_from_opt(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``--remove-solvent`` nulls the solvent on a project that has one."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "opt",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model is None
        assert settings.solvent_id is None

    def test_subcommand_level_solvent_overrides_group_level(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """Subcommand-level ``-sm``/``-si`` overrides group-level solvent."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.opt.ORCAOptJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cpcm",
                "-si",
                "toluene",
                "opt",
                "-sm",
                "smd",
                "-si",
                "water",
            ],
        )

        assert result.exit_code == 0, result.output
        # Subcommand-level smd/water overrides group-level cpcm/toluene
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"


class TestORCASolventCLITsCommand:
    """CLI solvent options propagated to the ``ts`` subcommand."""

    def test_solvent_model_and_id_injected_into_ts_settings(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm cpcm -si water`` sets solvent on the ts job settings."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.ts.ORCATSJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cpcm",
                "-si",
                "water",
                "ts",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCATSJob was never instantiated"
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id == "water"

    def test_remove_solvent_clears_solvent_from_ts(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``--remove-solvent`` removes solvent from ts job."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.ts.ORCATSJob",
            [
                "-p",
                "solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "ts",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings.solvent_model is None
        assert settings.solvent_id is None


class TestORCACpcmBlockOptions:
    """Tests for the ORCA-specific ``%cpcm`` block options via CLI ``-so``."""

    def test_custom_epsilon_no_solvent_id(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """Custom dielectric via ``--remove-solvent`` + ``sp -sm cpcm -so 'Epsilon 78.36'``.

        The project-level solvent (cyclohexane) is cleared by ``--remove-solvent``
        at the group level; the subcommand-level flags then set the custom
        dielectric without a named solvent.
        """
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "sp",
                "-sm",
                "cpcm",
                "-so",
                "Epsilon 78.36",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id is None
        assert settings.additional_solvent_options == "Epsilon 78.36"

    def test_custom_epsilon_and_refrac(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """Custom Epsilon + Refrac via ``--remove-solvent`` + ``sp -sm cpcm -so '...'``.

        Both ``Epsilon`` and ``Refrac`` are passed as a newline-separated
        string to ``-so``; each should appear in ``additional_solvent_options``.
        """
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "--remove-solvent",
                "sp",
                "-sm",
                "cpcm",
                "-so",
                "Epsilon 78.36\nRefrac 1.33",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "cpcm"
        assert settings.solvent_id is None
        assert "Epsilon 78.36" in settings.additional_solvent_options
        assert "Refrac 1.33" in settings.additional_solvent_options

    def test_smd_with_surface_type(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm smd -si water -so 'SurfaceType gepol_ses'`` stores all options."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "water",
                "-so",
                "SurfaceType gepol_ses",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "SurfaceType gepol_ses"

    def test_smd_with_rsolv(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """``-sm smd -si water -so 'Rsolv 1.30'`` stores Rsolv option."""
        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "smd",
                "-si",
                "water",
                "-so",
                "Rsolv 1.30",
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "smd"
        assert settings.solvent_id == "water"
        assert settings.additional_solvent_options == "Rsolv 1.30"

    def test_solventfilename_injected_into_sp_settings(
        self,
        tmp_path,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """-sf /path/water.cosmorsxyz stores solventfilename on sp settings."""
        # Create a dummy .cosmorsxyz file so click.Path(exists=True) is satisfied
        sf = tmp_path / "water.cosmorsxyz"
        sf.write_text("")

        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cosmors",
                "-si",
                "water",
                "-sf",
                str(sf),
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solvent_model == "cosmors"
        assert settings.solvent_id == "water"
        assert settings.solventfilename == str(sf)

    def test_solventfilename_group_level_injected_into_sp_settings(
        self,
        tmp_path,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
    ):
        """-sf at the orca group level propagates solventfilename to sp settings."""
        sf = tmp_path / "custom.cosmorsxyz"
        sf.write_text("")

        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_molecule_xyz_file,
                "-c",
                "0",
                "-m",
                "1",
                "-sm",
                "cosmors",
                "-sf",
                str(sf),
                "sp",
            ],
        )

        assert result.exit_code == 0, result.output
        assert settings is not None
        assert settings.solventfilename == str(sf)


class TestORCALabelAndAuxBasisOptions:
    def test_short_a_sets_append_label(self, single_molecule_xyz_file):
        from os.path import basename, splitext
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        runner = CliRunner()
        with patch(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob"
        ) as mock:
            mock.return_value = MagicMock()
            result = runner.invoke(
                orca_cli,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "-a",
                    "tag",
                    "sp",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock.call_args is not None
        base = splitext(basename(single_molecule_xyz_file))[0]
        assert mock.call_args.kwargs["label"].startswith(f"{base}_tag")
        assert mock.call_args.kwargs["settings"].aux_basis is None

    def test_short_B_sets_aux_basis(self, single_molecule_xyz_file):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        runner = CliRunner()
        with patch(
            "chemsmart.jobs.orca.singlepoint.ORCASinglePointJob"
        ) as mock:
            mock.return_value = MagicMock()
            result = runner.invoke(
                orca_cli,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "-c",
                    "0",
                    "-m",
                    "1",
                    "-B",
                    "def2/J",
                    "sp",
                ],
                obj={},
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert mock.call_args is not None
        assert mock.call_args.kwargs["settings"].aux_basis == "def2/J"


class TestORCAQMMMCLISpCommand:
    """CLI accepts structure files for ORCA QMMM submission."""

    def test_pdb_structure_file_accepted_for_ionic_crystal_qmmm(
        self,
        single_model_pdb_file,
        run_orca_and_capture_settings,
        tmpdir,
    ):
        """``.pdb`` inputs should not raise an unrecognised filetype error."""
        from unittest.mock import MagicMock

        prms_file = tmpdir.join("NaCl.ORCAFF.prms")
        prms_file.write("dummy force field parameters")

        result, settings = run_orca_and_capture_settings(
            "chemsmart.jobs.orca.qmmm.ORCAQMMMJob",
            [
                "-p",
                "gas_solv",
                "-f",
                single_model_pdb_file,
                "sp",
                "qmmm",
                "-j",
                "IONIC-CRYSTAL-QMMM",
                "-hx",
                "PBE",
                "-hb",
                "def2-SVP",
                "-lm",
                str(prms_file),
                "-ct",
                "0",
                "-ch",
                "19",
                "-mh",
                "1",
            ],
            ctx_obj={"jobrunner": MagicMock()},
        )

        assert result.exit_code == 0, result.output
        assert settings is not None, "ORCAQMMMJob was never instantiated"
        assert settings.jobtype == "IONIC-CRYSTAL-QMMM"
        assert settings.high_level_functional == "PBE"
        assert settings.low_level_method == str(prms_file)

    def test_conv_charges_false_survives_cli_reconstruction(
        self,
        single_molecule_xyz_file,
        run_orca_and_capture_settings,
        tmpdir,
    ):
        """``-cc false`` must round-trip when ``sub`` rebuilds the CLI."""
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli
        from chemsmart.utils.cli import CtxObjArguments

        prms_file = tmpdir.join("NaCl.ORCAFF.prms")
        prms_file.write("dummy force field parameters")

        ctx_obj = {"jobrunner": MagicMock()}
        runner = CliRunner()
        cli_args = [
            "-p",
            "gas_solv",
            "-f",
            single_molecule_xyz_file,
            "sp",
            "qmmm",
            "-j",
            "IONIC-CRYSTAL-QMMM",
            "-hx",
            "PBE",
            "-hb",
            "def2-SVP",
            "-lm",
            str(prms_file),
            "-cc",
            "false",
            "-ecp",
            "SDD",
            "-ct",
            "0",
            "-ch",
            "19",
            "-mh",
            "1",
        ]

        with patch("chemsmart.jobs.orca.qmmm.ORCAQMMMJob") as mock_job:
            mock_job.return_value = MagicMock()
            first = runner.invoke(
                orca_cli, cli_args, obj=ctx_obj, catch_exceptions=False
            )
            assert first.exit_code == 0, first.output

            reconstructed = CtxObjArguments(
                ctx_obj["subcommand"], entry_point="orca"
            ).reconstruct_command_line()[1:]
            cc_idx = reconstructed.index("--conv-charges")
            assert reconstructed[cc_idx + 1] == "false"

            second = runner.invoke(
                orca_cli,
                reconstructed,
                obj={"jobrunner": MagicMock()},
                catch_exceptions=False,
            )
            assert second.exit_code == 0, second.output
            settings = mock_job.call_args.kwargs["settings"]
            assert settings.conv_charges is False
            assert settings.ecp_layer_ecp == "SDD"


class TestORCAQMMMGeninputCLI:
    """``--geninput`` and CrystalPrep CLI flag parsing for ionic-crystal QMMM."""

    @staticmethod
    def _ionic_crystal_qmmm_args(
        structure_file,
        prms_file,
        *,
        extra=None,
    ):
        args = [
            "-p",
            "gas_solv",
            "-f",
            structure_file,
            "sp",
            "qmmm",
            "-j",
            "IONIC-CRYSTAL-QMMM",
            "-hx",
            "PBE",
            "-hb",
            "def2-SVP",
            "-lm",
            str(prms_file),
            "-ct",
            "0",
            "-ch",
            "19",
            "-mh",
            "1",
        ]
        if extra:
            args.extend(extra)
        return args

    @staticmethod
    def _geninput_flags(cif_file=None, scdimension="15x15x15"):
        flags = [
            "--geninput",
            "--cp-scdimension",
            scdimension,
            "--cp-atomtype",
            "Na",
            "0",
            "1.0",
            "0.0",
            "--cp-atomtype",
            "Cl",
            "1",
            "-1.0",
            "0.0",
        ]
        if cif_file is not None:
            flags.extend(["--cp-input-cif", str(cif_file)])
        return flags

    def test_geninput_rejected_without_ionic_crystal_jobtype(
        self,
        single_molecule_xyz_file,
        orca_ionic_crystal_nacl_cif_file,
        tmpdir,
    ):
        from unittest.mock import MagicMock

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        prms_file = tmpdir.join("NaCl.ORCAFF.prms")
        prms_file.write("dummy force field parameters")

        runner = CliRunner()
        args = self._ionic_crystal_qmmm_args(
            single_molecule_xyz_file,
            prms_file,
            extra=self._geninput_flags(
                cif_file=orca_ionic_crystal_nacl_cif_file
            ),
        )
        args[args.index("-j") + 1] = "QMMM"

        result = runner.invoke(
            orca_cli,
            args,
            obj={"jobrunner": MagicMock()},
            catch_exceptions=True,
        )

        assert result.exit_code != 0
        assert "--geninput is only supported" in result.output

    def test_geninput_requires_cif_input(
        self,
        orca_ionic_crystal_nacl_supercell_pdb_file,
        tmpdir,
    ):
        from unittest.mock import MagicMock

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        prms_file = tmpdir.join("NaCl.ORCAFF.prms")
        prms_file.write("dummy force field parameters")

        runner = CliRunner()
        args = self._ionic_crystal_qmmm_args(
            orca_ionic_crystal_nacl_supercell_pdb_file,
            prms_file,
            extra=self._geninput_flags(),
        )

        result = runner.invoke(
            orca_cli,
            args,
            obj={"jobrunner": MagicMock()},
            catch_exceptions=True,
        )

        assert result.exit_code != 0
        assert "requires a .cif input" in result.output

    def test_geninput_requires_scdimension(
        self,
        orca_ionic_crystal_nacl_supercell_pdb_file,
        orca_ionic_crystal_nacl_cif_file,
        tmpdir,
    ):
        from unittest.mock import MagicMock

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        prms_file = tmpdir.join("NaCl.ORCAFF.prms")
        prms_file.write("dummy force field parameters")

        runner = CliRunner()
        args = self._ionic_crystal_qmmm_args(
            orca_ionic_crystal_nacl_supercell_pdb_file,
            prms_file,
            extra=[
                "--geninput",
                "--cp-input-cif",
                orca_ionic_crystal_nacl_cif_file,
                "--cp-atomtype",
                "Na",
                "0",
                "1.0",
                "0.0",
            ],
        )

        result = runner.invoke(
            orca_cli,
            args,
            obj={"jobrunner": MagicMock()},
            catch_exceptions=True,
        )

        assert result.exit_code != 0
        assert "requires --cp-scdimension" in result.output

    def test_geninput_requires_atomtype(
        self,
        orca_ionic_crystal_nacl_supercell_pdb_file,
        orca_ionic_crystal_nacl_cif_file,
        tmpdir,
    ):
        from unittest.mock import MagicMock

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        prms_file = tmpdir.join("NaCl.ORCAFF.prms")
        prms_file.write("dummy force field parameters")

        runner = CliRunner()
        args = self._ionic_crystal_qmmm_args(
            orca_ionic_crystal_nacl_supercell_pdb_file,
            prms_file,
            extra=[
                "--geninput",
                "--cp-input-cif",
                orca_ionic_crystal_nacl_cif_file,
                "--cp-scdimension",
                "15x15x15",
            ],
        )

        result = runner.invoke(
            orca_cli,
            args,
            obj={"jobrunner": MagicMock()},
            catch_exceptions=True,
        )

        assert result.exit_code != 0
        assert "requires at least one --cp-atomtype" in result.output

    def test_geninput_parses_crystalprep_options_from_cli(
        self,
        orca_ionic_crystal_nacl_supercell_pdb_file,
        orca_ionic_crystal_nacl_cif_file,
        tmpdir,
    ):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        prms_file = tmpdir.join("NaCl.ORCAFF.prms")
        prms_file.write("dummy force field parameters")

        ctx_obj = {"jobrunner": MagicMock()}
        runner = CliRunner()
        args = self._ionic_crystal_qmmm_args(
            orca_ionic_crystal_nacl_supercell_pdb_file,
            prms_file,
            extra=self._geninput_flags(
                cif_file=orca_ionic_crystal_nacl_cif_file
            )
            + [
                "--cp-neutralize",
                "true",
                "--cp-qccharge",
                "0",
                "--cp-qcmult",
                "1",
                "--cp-template-out",
                "custom.cp.inp",
            ],
        )

        with patch("chemsmart.jobs.orca.qmmm.ORCAQMMMJob") as mock_job:
            mock_job.return_value = MagicMock()
            result = runner.invoke(
                orca_cli,
                args,
                obj=ctx_obj,
                catch_exceptions=False,
            )

        assert result.exit_code == 0, result.output
        assert ctx_obj["geninput"] is True
        options = ctx_obj["crystalprep_options"]
        assert options.input_cif == orca_ionic_crystal_nacl_cif_file
        assert options.sc_dimension == "15x15x15"
        assert len(options.atom_types) == 2
        assert options.atom_types[0].symbol == "Na"
        assert options.atom_types[1].symbol == "Cl"
        assert options.neutralize is True
        assert options.qc_charge == 0
        assert options.qc_mult == 1
        assert options.template_out == "custom.cp.inp"

    def test_resolve_crystalprep_input_cif_prefers_explicit_flag(
        self, orca_ionic_crystal_nacl_cif_file, tmpdir
    ):
        from chemsmart.cli.orca.qmmm import resolve_crystalprep_input_cif

        explicit = str(tmpdir.join("explicit.cif"))
        assert (
            resolve_crystalprep_input_cif(
                explicit, orca_ionic_crystal_nacl_cif_file
            )
            == explicit
        )

    def test_resolve_crystalprep_input_cif_uses_f_when_cif(
        self, orca_ionic_crystal_nacl_cif_file
    ):
        from chemsmart.cli.orca.qmmm import resolve_crystalprep_input_cif

        assert (
            resolve_crystalprep_input_cif(
                None, orca_ionic_crystal_nacl_cif_file
            )
            == orca_ionic_crystal_nacl_cif_file
        )

    def test_non_ionic_qmmm_unchanged_without_geninput(
        self,
        single_molecule_xyz_file,
        tmpdir,
    ):
        from unittest.mock import MagicMock, patch

        from click.testing import CliRunner

        from chemsmart.cli.orca.orca import orca as orca_cli

        prms_file = tmpdir.join("test.ORCAFF.prms")
        prms_file.write("dummy force field parameters")
        ctx_obj = {"jobrunner": MagicMock()}
        runner = CliRunner()

        with patch("chemsmart.jobs.orca.qmmm.ORCAQMMMJob") as mock_job:
            mock_job.return_value = MagicMock()
            result = runner.invoke(
                orca_cli,
                [
                    "-p",
                    "gas_solv",
                    "-f",
                    single_molecule_xyz_file,
                    "sp",
                    "qmmm",
                    "-j",
                    "QMMM",
                    "-hx",
                    "PBE",
                    "-hb",
                    "def2-SVP",
                    "-lm",
                    str(prms_file),
                ],
                obj=ctx_obj,
                catch_exceptions=False,
            )
            settings = mock_job.call_args.kwargs["settings"]

        assert result.exit_code == 0, result.output
        assert settings.jobtype == "QMMM"
        assert ctx_obj.get("geninput") is False
        assert ctx_obj.get("crystalprep_options") is None

    def test_build_crystalprep_options_from_cli_unit(
        self, orca_ionic_crystal_nacl_cif_file
    ):
        from chemsmart.cli.orca.qmmm import build_crystalprep_options_from_cli

        options = build_crystalprep_options_from_cli(
            cp_input_cif=None,
            structure_filename=orca_ionic_crystal_nacl_cif_file,
            cp_scdimension="2x2x2",
            cp_atomtype=[("Na", 0, 1.0, 0.0)],
            cp_docif=False,
            cp_dosupercell=False,
        )

        assert options.input_cif == orca_ionic_crystal_nacl_cif_file
        assert options.sc_dimension == "2x2x2"
        assert options.do_cif is False
        assert options.do_supercell is False
        assert options.do_embedding is True
        assert options.do_layers is True
