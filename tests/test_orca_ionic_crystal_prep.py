"""Tests for ionic-crystal QM/MM prep filename inference and command generation."""

import pytest

from chemsmart.jobs.orca.qmmm import (
    IonicCrystalPrepFilenames,
    ORCAQMMMJob,
    build_ionic_crystal_prep_commands,
    build_orca_crystalprep_commands,
    build_orca_crystalprep_geninput_command,
    build_orca_crystalprep_run_command,
    build_orca_mm_makeff_command,
)
from chemsmart.jobs.orca.writer import (
    CrystalPrepAtomType,
    CrystalPrepOptions,
)


class TestIonicCrystalPrepFilenames:
    def test_infer_default_filenames_from_cif_and_dimension(self):
        names = IonicCrystalPrepFilenames.infer("NaCl.cif", "15x15x15")

        assert names.stem == "NaCl.cif_15x15x15"
        assert names.supercell_xyz == "NaCl.cif_15x15x15.xyz"
        assert names.supercell_pdb == "NaCl.cif_15x15x15.pdb"
        assert names.orcaff_prms == "NaCl.cif_15x15x15.ORCAFF.prms"
        assert names.icqmmm_inp == "NaCl.cif_15x15x15.xyz.ICQMMM.inp"

    def test_infer_uses_cif_basename_from_path(self):
        names = IonicCrystalPrepFilenames.infer(
            "/data/crystals/NaCl.cif", "2x2x2"
        )

        assert names.stem == "NaCl.cif_2x2x2"
        assert names.supercell_xyz == "NaCl.cif_2x2x2.xyz"

    def test_infer_respects_individual_filename_overrides(self):
        names = IonicCrystalPrepFilenames.infer(
            "NaCl.cif",
            "15x15x15",
            supercell_xyz="custom.xyz",
            orcaff_prms="custom.prms",
        )

        assert names.stem == "NaCl.cif_15x15x15"
        assert names.supercell_xyz == "custom.xyz"
        assert names.supercell_pdb == "NaCl.cif_15x15x15.pdb"
        assert names.orcaff_prms == "custom.prms"

    def test_infer_respects_stem_override(self):
        names = IonicCrystalPrepFilenames.infer(
            "NaCl.cif",
            "15x15x15",
            stem="NaCl_custom_stem",
        )

        assert names.stem == "NaCl_custom_stem"
        assert names.supercell_xyz == "NaCl_custom_stem.xyz"
        assert names.icqmmm_inp == "NaCl_custom_stem.xyz.ICQMMM.inp"


class TestIonicCrystalPrepCommands:
    def test_build_orca_crystalprep_commands(self):
        geninput_cmd = build_orca_crystalprep_geninput_command("NaCl.cp.inp")
        run_cmd = build_orca_crystalprep_run_command("NaCl.cp.inp")

        assert geninput_cmd == 'ORCA_crystalprep "NaCl.cp.inp" -geninput'
        assert run_cmd == 'ORCA_crystalprep "NaCl.cp.inp"'
        assert build_orca_crystalprep_commands("NaCl.cp.inp") == (
            geninput_cmd,
            run_cmd,
        )

    def test_build_orca_mm_makeff_from_atom_types(self):
        atom_types = [
            CrystalPrepAtomType("Na", 0, 1.0, 0.0),
            CrystalPrepAtomType("Cl", 1, -1.0, 0.0),
        ]
        command = build_orca_mm_makeff_command(
            "NaCl.cif_15x15x15.xyz", atom_types
        )

        assert (
            command == 'ORCA_mm -makeff "NaCl.cif_15x15x15.xyz" '
            "-CEL Na 1.0 -CEL Cl -1.0"
        )

    def test_build_orca_mm_makeff_requires_atom_types(self):
        with pytest.raises(
            ValueError, match="At least one CrystalPrepAtomType"
        ):
            build_orca_mm_makeff_command("NaCl.cif_15x15x15.xyz", [])

    def test_build_ionic_crystal_prep_commands_in_order(self):
        atom_types = [
            CrystalPrepAtomType("Na", 0, 1.0, 0.0),
            CrystalPrepAtomType("Cl", 1, -1.0, 0.0),
        ]
        commands = build_ionic_crystal_prep_commands(
            "NaCl.cp.inp",
            "NaCl.cif_15x15x15.xyz",
            atom_types,
        )

        assert commands == (
            'ORCA_crystalprep "NaCl.cp.inp" -geninput',
            'ORCA_crystalprep "NaCl.cp.inp"',
            'ORCA_mm -makeff "NaCl.cif_15x15x15.xyz" '
            "-CEL Na 1.0 -CEL Cl -1.0",
        )

    def test_orca_qmmm_job_build_prep_commands_from_options(self):
        options = CrystalPrepOptions(
            input_cif="NaCl.cif",
            sc_dimension="15x15x15",
            atom_types=[
                CrystalPrepAtomType("Na", 0, 1.0, 0.0),
                CrystalPrepAtomType("Cl", 1, -1.0, 0.0),
            ],
        )
        filenames = ORCAQMMMJob.infer_prep_filenames_from_options(options)
        commands = ORCAQMMMJob.build_prep_commands_from_options(
            "NaCl.cp.inp", options
        )

        assert filenames.supercell_xyz == "NaCl.cif_15x15x15.xyz"
        assert commands[2].startswith(
            'ORCA_mm -makeff "NaCl.cif_15x15x15.xyz"'
        )
