import importlib
from pathlib import Path

import click
import pytest
from click.testing import CliRunner

from chemsmart.cli.run import run
from chemsmart.cli.sub import sub

_ENSEMBLE_SOLVENT_CONTRADICTIONS = (
    "one job per species",
    "one solvent SP per species",
    "lowest-energy optimized conformer",
    "use the lowest conformer",
    "solvent SPs stay",
    "solvent SPs remain",
    "Solvent single-points remain one",
)


def _assert_no_ensemble_solvent_contradiction(text):
    text = " ".join(text.split())
    for phrase in _ENSEMBLE_SOLVENT_CONTRADICTIONS:
        assert phrase not in text, phrase


def _assert_ensemble_solvent_help(text):
    normalized = " ".join(text.split())
    _assert_no_ensemble_solvent_contradiction(normalized)
    assert "matching solvent single-point" in normalized
    assert "num-cores" in normalized
    assert "_c1" in normalized


def _molecule_from_smiles(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from chemsmart.io.molecules.structure import Molecule

    rdkit_mol = Chem.MolFromSmiles(smiles)
    rdkit_mol = Chem.AddHs(rdkit_mol)
    AllChem.EmbedMolecule(rdkit_mol, randomSeed=0xC0FFEE)
    AllChem.UFFOptimizeMolecule(rdkit_mol)
    return Molecule.from_rdkit_mol(rdkit_mol)


def _write_signature_file(path: Path, program: str):
    signatures = {
        "gaussian": "Gaussian, Inc.\n",
        "orca": "* O   R   C   A *\n",
        "unknown": "Some random text\n",
    }
    path.write_text(signatures[program])


def _build_outputs(tmp_path: Path, program: str):
    names = [
        "ha.log",
        "a.log",
        "hb.log",
        "b.log",
        "has.log",
        "as.log",
        "hbs.log",
        "bs.log",
    ]
    files = {}
    for n in names:
        p = tmp_path / n
        _write_signature_file(p, program)
        files[n] = str(p)
    return files


def _invoke_pka_direct(runner, files, delta_g_proton=None):
    args = ["pka", "-s", "direct"]
    if delta_g_proton is not None:
        args.extend(["-dG", str(delta_g_proton)])
    args.extend(
        [
            "analyze",
            "-ha",
            files["ha.log"],
            "-a",
            files["a.log"],
            "-has",
            files["has.log"],
            "-as",
            files["as.log"],
        ]
    )
    return runner.invoke(run, args)


def _invoke_pka(runner, files):
    return runner.invoke(
        run,
        [
            "pka",
            "analyze",
            "-ha",
            files["ha.log"],
            "-a",
            files["a.log"],
            "-hr",
            files["hb.log"],
            "-r",
            files["b.log"],
            "-has",
            files["has.log"],
            "-as",
            files["as.log"],
            "--href-solv",
            files["hbs.log"],
            "--ref-solv",
            files["bs.log"],
            "-rp",
            "6.75",
        ],
    )


def _require_backend_pka_subcommand(command_group, backend):
    runner = CliRunner()
    result = runner.invoke(command_group, [backend, "--help"])
    assert result.exit_code == 0, result.output
    if "\n  pka" not in result.output:
        pytest.skip(
            f"{backend} backend pka subcommand is not registered in this build."
        )


class _FakeThermochemistry:
    def __init__(self, filename, **kwargs):
        self.filename = filename
        self.electronic_energy = -627.0
        self.qrrho_gibbs_free_energy = -628.0


def _install_fake_thermochemistry(monkeypatch, constructed=None):
    constructed = [] if constructed is None else constructed

    class _TrackingFakeThermochemistry(_FakeThermochemistry):
        def __init__(self, filename, **kwargs):
            constructed.append(Path(filename).name)
            super().__init__(filename, **kwargs)

    monkeypatch.setattr(
        "chemsmart.cli.pka.Thermochemistry",
        _TrackingFakeThermochemistry,
    )
    return constructed


def _write_test_backend_project(tmp_path, backend):
    config_root = tmp_path / "chemsmart_cfg"
    backend_cfg_dir = config_root / backend
    backend_cfg_dir.mkdir(parents=True)
    (backend_cfg_dir / "test.yaml").write_text(
        "gas:\n"
        "  functional: B3LYP\n"
        "  basis: def2-SVP\n"
        "solv:\n"
        "  functional: B3LYP\n"
        "  basis: def2-SVP\n"
        "  freq: false\n"
        "  solvent_model: smd\n"
        "  solvent_id: water\n"
    )
    return config_root


def _setup_sub_pka_batch_test(tmp_path, monkeypatch, backend, captured):
    """Shared fixtures for sub ... pka batch submission tests."""
    acid1 = tmp_path / "acid1.xyz"
    acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    acid2 = tmp_path / "acid2.xyz"
    acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

    table = tmp_path / "pka_scale.csv"
    table.write_text(
        "structure,filepath,proton_index,charge,multiplicity\n"
        f"acid1,{acid1},2,0,1\n"
        f"acid2,{acid2},2,1,2\n"
    )

    config_root = _write_test_backend_project(tmp_path, backend)
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

    from chemsmart.settings.server import Server

    fake_server = Server(name="dummy")
    captured["submissions"] = []
    fake_server.submit = lambda job, test=False, cli_args=None, **kw: captured[
        "submissions"
    ].append((job, test, cli_args))
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )
    return table, captured


def _capture_sub_pka_jobs(tmp_path, monkeypatch, backend, extra_sub_args):
    """Run ``chemsmart sub ... pka`` in test mode and return created jobs."""
    _require_backend_pka_subcommand(sub, backend)
    config_root = _write_test_backend_project(tmp_path, backend)
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

    from chemsmart.settings.server import Server

    fake_server = Server(name="dummy")
    captured = {"jobs": []}
    fake_server.submit = lambda job, test=False, cli_args=None, **kw: captured[
        "jobs"
    ].append(job)
    monkeypatch.setattr(
        "chemsmart.settings.server.Server.from_servername",
        lambda _name: fake_server,
    )

    runner = CliRunner()
    result = runner.invoke(
        sub,
        [
            "--test",
            "--server",
            "dummy",
            "--no-scratch",
            backend,
            "-p",
            "test",
            *extra_sub_args,
        ],
    )
    return result, captured["jobs"]


def _build_pka_batch_table(tmp_path):
    acid1 = tmp_path / "acid1.xyz"
    acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
    acid2 = tmp_path / "acid2.xyz"
    acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

    table = tmp_path / "pka_scale.csv"
    table.write_text(
        "filepath,proton_index,charge,multiplicity\n"
        f"{acid1},2,0,1\n"
        f"{acid2},2,1,2\n"
    )
    return table


def _first_hydrogen_index(molecule):
    return next(
        i + 1 for i, symbol in enumerate(molecule.symbols) if symbol == "H"
    )


def _write_crest_main_out(path, normal=True):
    text = (
        "Conformer-Rotamer Ensemble Sampling Tool\n"
        "$ crest mol.xyz --chrg 0 --uhf 0\n"
    )
    if normal:
        text += "CREST terminated normally.\n"
    else:
        text += "CREST failed.\n"
    Path(path).write_text(text)


def _write_translated_conformers(path, molecule, n, shift=0.5):
    for index in range(n):
        frame = molecule.copy()
        frame.positions = frame.positions + [index * shift, 0.0, 0.0]
        frame.write_xyz(str(path), mode="w" if index == 0 else "a")


def _direct_pka_job(
    backend,
    molecule,
    jobrunner,
    *,
    sampling=False,
    num_conformers=1,
    label="mol_pka",
    scheme="direct",
    reference_file=None,
    reference_proton_index=None,
):
    proton_index = _first_hydrogen_index(molecule)
    molecule.charge = 0
    molecule.multiplicity = 1
    reference_kwargs = {}
    if reference_file is not None:
        reference_kwargs = {
            "reference_file": reference_file,
            "reference_proton_index": reference_proton_index,
            "reference_charge": 0,
            "reference_multiplicity": 1,
        }
    if backend == "gaussian":
        from chemsmart.jobs.gaussian.pka import GaussianpKaJob
        from chemsmart.jobs.gaussian.settings import GaussianpKaJobSettings

        settings = GaussianpKaJobSettings(
            proton_index=proton_index,
            scheme=scheme,
            functional="B3LYP",
            basis="6-31G*",
            sampling=sampling,
            num_conformers=num_conformers,
            **reference_kwargs,
        )
        job_cls = GaussianpKaJob
    else:
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme=scheme,
            functional="B3LYP",
            basis="def2-SVP",
            sampling=sampling,
            num_conformers=num_conformers,
            **reference_kwargs,
        )
        job_cls = ORCApKaJob
    return job_cls(
        molecule=molecule,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
    )


def _make_crest_search_job(molecule, jobrunner, label="ha_crest"):
    from chemsmart.jobs.crest.conformers import CRESTConformerSearchJob
    from chemsmart.jobs.crest.settings import CRESTJobSettings

    settings = CRESTJobSettings.default()
    settings.jobtype = "conformers"
    settings.charge = 0 if molecule.charge is None else molecule.charge
    settings.multiplicity = (
        1 if molecule.multiplicity is None else molecule.multiplicity
    )
    return CRESTConformerSearchJob(
        molecule=molecule,
        settings=settings,
        label=label,
        jobrunner=jobrunner,
    )


class TestAqueousProtonSolutionFreeEnergy:
    def test_value_at_298_15_k(self):
        from chemsmart.cli.pka import (
            aqueous_proton_solution_free_energy_kcal_mol,
        )

        g_soln = aqueous_proton_solution_free_energy_kcal_mol(298.15)
        assert g_soln == pytest.approx(-270.3, abs=0.05)

    def test_temperature_dependence_of_gas_and_standard_state_terms(self):
        import math

        from chemsmart.cli.pka import (
            aqueous_proton_solution_free_energy_kcal_mol,
        )
        from chemsmart.utils.constants import R, atm_to_pa, energy_conversion

        delta_g_solv = -265.9
        entropy_kcal_mol_k = 26.016 / 1000.0
        r_kcal_mol_k = energy_conversion("j/mol", "kcal/mol", R)
        r_liter_atm_mol_k = R / atm_to_pa * 1000.0

        def gas_plus_standard_state(temperature):
            g_gas = (
                2.5 * r_kcal_mol_k * temperature
                - temperature * entropy_kcal_mol_k
            )
            g_std = (
                r_kcal_mol_k
                * temperature
                * math.log(r_liter_atm_mol_k * temperature)
            )
            return g_gas + g_std

        g_298 = aqueous_proton_solution_free_energy_kcal_mol(
            298.15, delta_g_solv=delta_g_solv
        )
        g_373 = aqueous_proton_solution_free_energy_kcal_mol(
            373.15, delta_g_solv=delta_g_solv
        )
        assert g_298 - delta_g_solv == pytest.approx(
            gas_plus_standard_state(298.15)
        )
        assert g_373 - delta_g_solv == pytest.approx(
            gas_plus_standard_state(373.15)
        )
        assert g_373 != pytest.approx(g_298)

    def test_compute_pka_direct_uses_computed_default(
        self, tmp_path, monkeypatch
    ):
        files = _build_outputs(tmp_path, "gaussian")
        _install_fake_thermochemistry(monkeypatch)
        from chemsmart.cli.pka import (
            aqueous_proton_solution_free_energy_kcal_mol,
            compute_pka,
        )

        result = compute_pka(
            ha_gas_file=files["ha.log"],
            a_gas_file=files["a.log"],
            ha_solv_file=files["has.log"],
            a_solv_file=files["as.log"],
            scheme="direct",
            temperature=298.15,
        )
        expected = aqueous_proton_solution_free_energy_kcal_mol(298.15)
        assert result["delta_G_proton_kcal_mol"] == pytest.approx(expected)
        assert result["delta_G_proton_user_supplied"] is False

    def test_compute_pka_direct_honors_user_override(
        self, tmp_path, monkeypatch
    ):
        files = _build_outputs(tmp_path, "gaussian")
        _install_fake_thermochemistry(monkeypatch)
        from chemsmart.cli.pka import compute_pka

        result = compute_pka(
            ha_gas_file=files["ha.log"],
            a_gas_file=files["a.log"],
            ha_solv_file=files["has.log"],
            a_solv_file=files["as.log"],
            scheme="direct",
            temperature=298.15,
            delta_G_proton=-270.0,
        )
        assert result["delta_G_proton_kcal_mol"] == -270.0
        assert result["delta_G_proton_user_supplied"] is True


class TestPkaEnsembleAnalysis:
    """Ensemble G_eff analysis for multi-conformer pKa outputs."""

    def test_ensemble_effective_free_energy_two_equal_g(self):
        import math

        from chemsmart.cli.pka import ensemble_effective_free_energy
        from chemsmart.utils.constants import R, energy_conversion

        g = -1.0
        temperature = 298.15
        g_eff = ensemble_effective_free_energy([g, g], temperature)
        rt_hartree = energy_conversion("j/mol", "hartree", R * temperature)
        assert g_eff == pytest.approx(g - rt_hartree * math.log(2))

    def test_ensemble_effective_free_energy_single_value_unchanged(self):
        from chemsmart.cli.pka import ensemble_effective_free_energy

        assert ensemble_effective_free_energy([-0.5], 298.15) == pytest.approx(
            -0.5
        )

    def test_compute_pka_two_file_pairs_uses_g_eff(
        self, tmp_path, monkeypatch
    ):
        import math

        files = _build_outputs(tmp_path, "gaussian")
        _install_fake_thermochemistry(monkeypatch)
        ha2 = tmp_path / "ha2.log"
        a2 = tmp_path / "a2.log"
        has2 = tmp_path / "has2.log"
        as2 = tmp_path / "as2.log"
        for path in (ha2, a2, has2, as2):
            _write_signature_file(path, "gaussian")

        from chemsmart.cli.pka import compute_pka
        from chemsmart.utils.constants import R, energy_conversion

        temperature = 298.15
        single = compute_pka(
            ha_gas_file=files["ha.log"],
            a_gas_file=files["a.log"],
            ha_solv_file=files["has.log"],
            a_solv_file=files["as.log"],
            scheme="direct",
            temperature=temperature,
        )
        result = compute_pka(
            ha_gas_file=[files["ha.log"], str(ha2)],
            a_gas_file=[files["a.log"], str(a2)],
            ha_solv_file=[files["has.log"], str(has2)],
            a_solv_file=[files["as.log"], str(as2)],
            scheme="direct",
            temperature=temperature,
        )
        rt_hartree = energy_conversion("j/mol", "hartree", R * temperature)
        expected_g_ha = single["G_soln_HA_au"] - rt_hartree * math.log(2)
        expected_g_a = single["G_soln_A_au"] - rt_hartree * math.log(2)
        assert result["G_soln_HA_au"] == pytest.approx(expected_g_ha)
        assert result["G_soln_A_au"] == pytest.approx(expected_g_a)
        assert result["num_conformers_HA"] == 2
        assert result["num_conformers_A"] == 2

    def test_print_pka_summary_mentions_ensemble_g_eff(
        self, tmp_path, monkeypatch, capsys
    ):
        files = _build_outputs(tmp_path, "gaussian")
        _install_fake_thermochemistry(monkeypatch)
        ha2 = tmp_path / "ha2.log"
        a2 = tmp_path / "a2.log"
        has2 = tmp_path / "has2.log"
        as2 = tmp_path / "as2.log"
        for path in (ha2, a2, has2, as2):
            _write_signature_file(path, "gaussian")

        from chemsmart.cli.pka import print_pka_summary

        print_pka_summary(
            ha_gas_file=[files["ha.log"], str(ha2)],
            a_gas_file=[files["a.log"], str(a2)],
            ha_solv_file=[files["has.log"], str(has2)],
            a_solv_file=[files["as.log"], str(as2)],
            scheme="direct",
            temperature=298.15,
        )
        output = capsys.readouterr().out
        assert "G_eff" in output
        assert "2 conformers" in output

    def test_analyze_auto_discovers_ensemble(self, tmp_path, monkeypatch):
        for name in (
            "acid1_pka_HA_opt_c1.log",
            "acid1_pka_HA_opt_c2.log",
            "acid1_pka_A_opt_c1.log",
            "acid1_pka_A_opt_c2.log",
            "acid1_pka_HA_sp_c1.log",
            "acid1_pka_HA_sp_c2.log",
            "acid1_pka_A_sp_c1.log",
            "acid1_pka_A_sp_c2.log",
        ):
            _write_signature_file(tmp_path / name, "gaussian")

        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)
        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "-s",
                "direct",
                "analyze",
                "-ha",
                str(tmp_path / "acid1_pka_HA_opt_c1.log"),
            ],
        )
        assert result.exit_code == 0, result.output
        assert called["kwargs"]["ha_gas_file"] == [
            str(tmp_path / "acid1_pka_HA_opt_c1.log"),
            str(tmp_path / "acid1_pka_HA_opt_c2.log"),
        ]
        assert called["kwargs"]["a_gas_file"] == [
            str(tmp_path / "acid1_pka_A_opt_c1.log"),
            str(tmp_path / "acid1_pka_A_opt_c2.log"),
        ]
        assert called["kwargs"]["ha_solv_file"] == [
            str(tmp_path / "acid1_pka_HA_sp_c1.log"),
            str(tmp_path / "acid1_pka_HA_sp_c2.log"),
        ]
        assert called["kwargs"]["a_solv_file"] == [
            str(tmp_path / "acid1_pka_A_sp_c1.log"),
            str(tmp_path / "acid1_pka_A_sp_c2.log"),
        ]


class TestPkaCrestSampling:
    """CREST sampling helpers, pKa job labels, and conformer extract path."""

    def test_pka_subjob_label_n1_and_n_greater_than_1(self):
        from chemsmart.utils.datasets import pka_subjob_label

        assert (
            pka_subjob_label("mol_pka", "HA", "opt", 1, 1) == "mol_pka_HA_opt"
        )
        assert pka_subjob_label("mol_pka", "A", "sp", 2, 1) == "mol_pka_A_sp"
        assert (
            pka_subjob_label("mol_pka", "HA", "opt", 1, 3)
            == "mol_pka_HA_opt_c1"
        )
        assert (
            pka_subjob_label("mol_pka", "A", "sp", 2, 3) == "mol_pka_A_sp_c2"
        )

    def test_select_crest_conformers_waits_when_output_missing(
        self, temporary_working_dir, water_molecule, crest_jobrunner_no_scratch
    ):
        from chemsmart.cli.pka import select_crest_conformers

        crest_job = _make_crest_search_job(
            water_molecule, crest_jobrunner_no_scratch
        )
        assert select_crest_conformers(crest_job, 1, water_molecule) is None

    def test_select_crest_conformers_missing_xyz_falls_back(
        self,
        temporary_working_dir,
        water_molecule,
        crest_jobrunner_no_scratch,
        caplog,
    ):
        import logging

        from chemsmart.cli.pka import select_crest_conformers

        water_molecule.charge = 0
        water_molecule.multiplicity = 1
        crest_job = _make_crest_search_job(
            water_molecule, crest_jobrunner_no_scratch
        )
        _write_crest_main_out(crest_job.outputfile, normal=True)
        with caplog.at_level(logging.WARNING):
            selected = select_crest_conformers(crest_job, 1, water_molecule)
        assert len(selected) == 1
        assert selected[0].positions == pytest.approx(water_molecule.positions)
        assert selected[0].charge == 0
        assert "input" in caplog.text.lower()

    def test_select_crest_conformers_abnormal_termination_falls_back(
        self,
        temporary_working_dir,
        water_molecule,
        crest_jobrunner_no_scratch,
        caplog,
    ):
        import logging

        from chemsmart.cli.pka import select_crest_conformers

        water_molecule.charge = 1
        water_molecule.multiplicity = 1
        crest_job = _make_crest_search_job(
            water_molecule, crest_jobrunner_no_scratch, label="failed_crest"
        )
        _write_crest_main_out(crest_job.outputfile, normal=False)
        with caplog.at_level(logging.WARNING):
            selected = select_crest_conformers(crest_job, 1, water_molecule)
        assert len(selected) == 1
        assert selected[0].positions == pytest.approx(water_molecule.positions)
        assert selected[0].charge == 1
        assert "did not terminate normally" in caplog.text

    def test_select_crest_conformers_fewer_than_n_uses_available(
        self,
        temporary_working_dir,
        water_molecule,
        crest_jobrunner_no_scratch,
        caplog,
    ):
        import logging

        from chemsmart.cli.pka import select_crest_conformers

        water_molecule.charge = 0
        water_molecule.multiplicity = 1
        crest_job = _make_crest_search_job(
            water_molecule, crest_jobrunner_no_scratch, label="short_crest"
        )
        _write_crest_main_out(crest_job.outputfile, normal=True)
        _write_translated_conformers(
            Path(crest_job.folder) / "crest_conformers.xyz",
            water_molecule,
            2,
        )
        with caplog.at_level(logging.WARNING):
            selected = select_crest_conformers(crest_job, 3, water_molecule)
        assert len(selected) == 2
        assert selected[0].charge == 0
        assert "requested 3" in caplog.text

    def test_select_crest_conformers_slices_fixture_ensemble(
        self,
        temporary_working_dir,
        water_molecule,
        crest_jobrunner_no_scratch,
        multiple_molecules_xyz_file,
    ):
        import shutil

        from chemsmart.cli.pka import select_crest_conformers
        from chemsmart.io.molecules.structure import Molecule

        water_molecule.charge = 0
        water_molecule.multiplicity = 1
        crest_job = _make_crest_search_job(
            water_molecule, crest_jobrunner_no_scratch, label="ens_crest"
        )
        _write_crest_main_out(crest_job.outputfile, normal=True)
        shutil.copy(
            multiple_molecules_xyz_file,
            Path(crest_job.folder) / "crest_conformers.xyz",
        )
        selected = select_crest_conformers(crest_job, 3, water_molecule)
        expected = Molecule.from_filepath(
            multiple_molecules_xyz_file, index=":", return_list=True
        )[:3]
        assert len(selected) == 3
        for got, want in zip(selected, expected):
            assert got.positions == pytest.approx(want.positions)
            assert got.charge == 0

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_default_path_has_no_crest_jobs_and_keeps_labels(
        self,
        backend,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        orca_jobrunner_no_scratch,
    ):
        from chemsmart.io.molecules.structure import Molecule

        jobrunner = {
            "gaussian": gaussian_jobrunner_no_scratch,
            "orca": orca_jobrunner_no_scratch,
        }[backend]
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        job = _direct_pka_job(backend, mol, jobrunner, label="1a_pka")
        assert job.crest_jobs == []
        assert job.protonated_crest_job is None
        assert job.conjugate_base_crest_job is None
        assert job.protonated_job.label == "1a_pka_HA_opt"
        assert job.conjugate_base_job.label == "1a_pka_A_opt"
        if job.sp_jobs is None:
            job._create_sp_jobs()
        assert job.protonated_sp_job.label == "1a_pka_HA_sp"
        assert job.conjugate_base_sp_job.label == "1a_pka_A_sp"

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sampling_n1_creates_crest_children_and_keeps_dft_labels(
        self,
        backend,
        temporary_working_dir,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        orca_jobrunner_no_scratch,
    ):
        from chemsmart.io.molecules.structure import Molecule

        jobrunner = {
            "gaussian": gaussian_jobrunner_no_scratch,
            "orca": orca_jobrunner_no_scratch,
        }[backend]
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        job = _direct_pka_job(
            backend, mol, jobrunner, sampling=True, label="1a_pka"
        )
        assert len(job.crest_jobs) == 2
        assert job.protonated_crest_job.label == "1a_pka_HA_crest"
        assert job.conjugate_base_crest_job.label == "1a_pka_A_crest"
        assert job.protonated_job.label == "1a_pka_HA_opt"
        assert job.conjugate_base_job.label == "1a_pka_A_opt"
        if job.sp_jobs is None:
            job._create_sp_jobs()
        assert job.protonated_sp_job.label == "1a_pka_HA_sp"
        assert job.conjugate_base_sp_job.label == "1a_pka_A_sp"

    def test_sampling_n1_uses_extracted_crest_geometry(
        self,
        temporary_working_dir,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
    ):
        import numpy as np

        from chemsmart.io.molecules.structure import Molecule

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        job = _direct_pka_job(
            "gaussian",
            mol,
            gaussian_jobrunner_no_scratch,
            sampling=True,
            label="1a_pka",
        )
        original = np.array(mol.positions, copy=True)
        shifted = mol.copy()
        shifted.positions = original + np.array([1.0, 0.0, 0.0])
        ha_crest = job.protonated_crest_job
        a_crest = job.conjugate_base_crest_job
        _write_crest_main_out(ha_crest.outputfile, normal=True)
        _write_crest_main_out(a_crest.outputfile, normal=True)
        best_path = Path(ha_crest.folder) / "crest_best.xyz"
        shifted.write_xyz(str(best_path), mode="w")

        selected = job._selected_crest_conformers()
        assert selected is not None
        ha_confs, a_confs, *_ = selected
        assert ha_confs[0].positions == pytest.approx(shifted.positions)
        assert not np.allclose(ha_confs[0].positions, original)
        job._prepare_target_opt_jobs(ha_confs, a_confs)
        assert job.protonated_job.molecule.positions == pytest.approx(
            shifted.positions
        )
        assert job.protonated_job.label == "1a_pka_HA_opt"

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sampling_n3_labels_opt_c1_to_c3(
        self,
        backend,
        temporary_working_dir,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        orca_jobrunner_no_scratch,
    ):
        from chemsmart.io.molecules.structure import Molecule

        jobrunner = {
            "gaussian": gaussian_jobrunner_no_scratch,
            "orca": orca_jobrunner_no_scratch,
        }[backend]
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        job = _direct_pka_job(
            backend,
            mol,
            jobrunner,
            sampling=True,
            num_conformers=3,
            label="1a_pka",
        )
        assert [child.label for child in job.protonated_opt_jobs] == [
            "1a_pka_HA_opt_c1",
            "1a_pka_HA_opt_c2",
            "1a_pka_HA_opt_c3",
        ]
        assert [child.label for child in job.conjugate_base_opt_jobs] == [
            "1a_pka_A_opt_c1",
            "1a_pka_A_opt_c2",
            "1a_pka_A_opt_c3",
        ]
        assert job.protonated_job.label == "1a_pka_HA_opt_c1"
        if job.sp_jobs is None:
            job._create_sp_jobs()
        assert [child.label for child in job.protonated_sp_jobs] == [
            "1a_pka_HA_sp_c1",
            "1a_pka_HA_sp_c2",
            "1a_pka_HA_sp_c3",
        ]
        assert [child.label for child in job.conjugate_base_sp_jobs] == [
            "1a_pka_A_sp_c1",
            "1a_pka_A_sp_c2",
            "1a_pka_A_sp_c3",
        ]

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_submit_rejects_num_conformers_without_sampling(
        self, tmp_path, monkeypatch, backend
    ):
        _require_backend_pka_subcommand(run, backend)
        acid = tmp_path / "acid.xyz"
        acid.write_text("2\nacid\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "--no-scratch",
                "--fake",
                backend,
                "-p",
                "test",
                "-f",
                str(acid),
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-N",
                "3",
                "-s",
                "direct",
                "submit",
            ],
        )
        assert result.exit_code != 0
        assert "-N/--num-conformers requires --sampling" in result.output


def _conformer_index(filename):
    suffix = Path(filename).stem.rsplit("_c", 1)[-1]
    if suffix.isdigit() and suffix != Path(filename).stem:
        return int(suffix)
    return 1


def _conformer_mock_energies(filename):
    """Return distinct J/mol energies so each conformer changes G_soln.

    HA and A- use different conformer slopes so the ensemble ΔG is not the
    same as the conformer-1 ΔG.
    """
    name = Path(filename).name
    index = _conformer_index(name)
    slope = 8000.0 if "_HA_" in name else 400.0
    electronic = -1.0e6 - index * slope
    if "_sp" in name:
        electronic -= 5.0e4
    return electronic, electronic + 1000.0


def _install_conformer_thermochemistry(monkeypatch):
    class _FakeThermochemistry:
        def __init__(self, filename, **kwargs):
            electronic, qh = _conformer_mock_energies(filename)
            self.electronic_energy = electronic
            self.qrrho_gibbs_free_energy = qh
            self.zero_point_energy = electronic
            self.enthalpy = electronic
            self.qrrho_enthalpy = qh
            self.gibbs_free_energy = qh

    monkeypatch.setattr(
        "chemsmart.cli.pka.Thermochemistry",
        _FakeThermochemistry,
    )


def _touch_output(path, program):
    signature = (
        "Gaussian, Inc.\n" if program == "gaussian" else "* O   R   C   A *\n"
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(signature)


class TestPkaEnsembleJobOutputs:
    """Job API, analyze, and batch-analyze share one ensemble file set."""

    def test_conformer_paths_sort_c10_after_c9(self):
        from chemsmart.utils.datasets import output_paths_in_conformer_order

        class _Job:
            def __init__(self, outputfile):
                self.outputfile = outputfile

        ordered = output_paths_in_conformer_order(
            [
                _Job("mol_HA_opt_c10.log"),
                _Job("mol_HA_opt_c2.log"),
                _Job("mol_HA_opt_c9.log"),
                _Job("mol_HA_opt_c1.log"),
            ]
        )
        assert [Path(path).name for path in ordered] == [
            "mol_HA_opt_c1.log",
            "mol_HA_opt_c2.log",
            "mol_HA_opt_c9.log",
            "mol_HA_opt_c10.log",
        ]

    def test_mismatched_gas_and_sp_counts_name_the_species(self):
        from chemsmart.utils.datasets import collect_pka_species_outputs

        class _Job:
            def __init__(self, outputfile):
                self.outputfile = outputfile

        with pytest.raises(
            ValueError,
            match=r"HA gas-phase and solvent file counts must match \(2 vs 1\)",
        ):
            collect_pka_species_outputs(
                ha_opt_jobs=[_Job("ha_c1.log"), _Job("ha_c2.log")],
                a_opt_jobs=[_Job("a_c1.log")],
                ha_sp_jobs=[_Job("ha_sp_c1.log")],
                a_sp_jobs=[_Job("a_sp_c1.log")],
            )

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_one_conformer_paths_stay_legacy(
        self,
        backend,
        temporary_working_dir,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        orca_jobrunner_no_scratch,
        monkeypatch,
    ):
        from chemsmart.cli.pka import compute_pka
        from chemsmart.io.molecules.structure import Molecule

        _install_conformer_thermochemistry(monkeypatch)
        jobrunner = {
            "gaussian": gaussian_jobrunner_no_scratch,
            "orca": orca_jobrunner_no_scratch,
        }[backend]
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        job = _direct_pka_job(backend, mol, jobrunner, label="legacy_pka")
        monkeypatch.setattr(job, "_opt_jobs_are_complete", lambda: True)
        files = job._pka_output_files()
        assert list(files) == ["HA", "A-"]
        for species in ("HA", "A-"):
            assert len(files[species]["gas"]) == 1
            assert len(files[species]["solv"]) == 1
            assert "_c" not in Path(files[species]["gas"][0]).name
            assert "_c" not in Path(files[species]["solv"][0]).name

        job_result = job.print_thermochemistry()
        analyzed = compute_pka(
            ha_gas_file=files["HA"]["gas"][0],
            a_gas_file=files["A-"]["gas"][0],
            ha_solv_file=files["HA"]["solv"][0],
            a_solv_file=files["A-"]["solv"][0],
            scheme="direct",
        )
        assert job_result["pKa"] == pytest.approx(analyzed["pKa"])
        assert job_result["num_conformers_HA"] == 1
        thermo = job.compute_thermochemistry()
        assert isinstance(thermo["HA"], dict)
        assert thermo["HA"]["name"] == "HA"

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_three_conformer_job_passes_three_gas_and_sp_paths(
        self,
        backend,
        temporary_working_dir,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        orca_jobrunner_no_scratch,
    ):
        from chemsmart.io.molecules.structure import Molecule

        jobrunner = {
            "gaussian": gaussian_jobrunner_no_scratch,
            "orca": orca_jobrunner_no_scratch,
        }[backend]
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        job = _direct_pka_job(
            backend,
            mol,
            jobrunner,
            sampling=True,
            num_conformers=3,
            label="acid_pka",
        )
        files = job._pka_output_files()
        for species, gas_tag, solv_tag in (
            ("HA", "HA_opt", "HA_sp"),
            ("A-", "A_opt", "A_sp"),
        ):
            gas_names = [Path(path).name for path in files[species]["gas"]]
            solv_names = [Path(path).name for path in files[species]["solv"]]
            assert gas_names == [
                f"acid_pka_{gas_tag}_c1.{_output_ext(backend)}",
                f"acid_pka_{gas_tag}_c2.{_output_ext(backend)}",
                f"acid_pka_{gas_tag}_c3.{_output_ext(backend)}",
            ]
            assert solv_names == [
                f"acid_pka_{solv_tag}_c1.{_output_ext(backend)}",
                f"acid_pka_{solv_tag}_c2.{_output_ext(backend)}",
                f"acid_pka_{solv_tag}_c3.{_output_ext(backend)}",
            ]

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_reference_ensembles_include_all_four_species(
        self,
        backend,
        temporary_working_dir,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        orca_jobrunner_no_scratch,
    ):
        from chemsmart.io.molecules.structure import Molecule

        jobrunner = {
            "gaussian": gaussian_jobrunner_no_scratch,
            "orca": orca_jobrunner_no_scratch,
        }[backend]
        reference = temporary_working_dir / "ref_acid.xyz"
        reference.write_text("2\nref\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        job = _direct_pka_job(
            backend,
            mol,
            jobrunner,
            sampling=True,
            num_conformers=3,
            label="acid_pka",
            scheme="proton exchange",
            reference_file=str(reference),
            reference_proton_index=2,
        )
        files = job._pka_output_files()
        assert list(files) == ["HA", "A-", "HRef", "Ref-"]
        for species in files:
            assert len(files[species]["gas"]) == 3
            assert len(files[species]["solv"]) == 3
            gas_indexes = [
                Path(path).stem.rsplit("_c", 1)[1]
                for path in files[species]["gas"]
            ]
            solv_indexes = [
                Path(path).stem.rsplit("_c", 1)[1]
                for path in files[species]["solv"]
            ]
            assert gas_indexes == ["1", "2", "3"]
            assert solv_indexes == ["1", "2", "3"]

    def test_gaussian_sp_phase_rebuilds_after_early_output_collection(
        self,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        monkeypatch,
    ):
        """An early path collection must not freeze solvent jobs to the input geometry."""
        import numpy as np

        from chemsmart.io.molecules.structure import Molecule

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        job = _direct_pka_job(
            "gaussian", mol, gaussian_jobrunner_no_scratch, label="acid_pka"
        )
        input_positions = np.array(job.protonated_sp_job.molecule.positions)
        job._pka_output_files()
        assert job.protonated_sp_job.molecule.positions == pytest.approx(
            input_positions
        )

        shifted = mol.copy()
        shifted.positions = input_positions + np.array([1.0, 0.0, 0.0])

        class _Optimized:
            normal_termination = True
            molecule = shifted

        for opt_job in job.opt_jobs:
            opt_job.is_complete = lambda: True
            opt_job._output = lambda: _Optimized()

        monkeypatch.setattr(
            "chemsmart.jobs.chain.pka.run_phase_jobs", lambda **kwargs: None
        )
        job._run_sp_jobs()
        assert job.protonated_sp_job.molecule.positions == pytest.approx(
            shifted.positions
        )
        assert not np.allclose(
            job.protonated_sp_job.molecule.positions, input_positions
        )

    def test_job_analyze_and_batch_analyze_share_ensemble_pka(
        self,
        temporary_working_dir,
        single_molecule_xyz_file,
        gaussian_jobrunner_no_scratch,
        monkeypatch,
    ):
        from chemsmart.cli.pka import compute_pka
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.utils.datasets import PKaOutputTable

        _install_conformer_thermochemistry(monkeypatch)
        mol = Molecule.from_filepath(single_molecule_xyz_file)
        job = _direct_pka_job(
            "gaussian",
            mol,
            gaussian_jobrunner_no_scratch,
            sampling=True,
            num_conformers=3,
            label="acid_pka",
        )
        monkeypatch.setattr(job, "_opt_jobs_are_complete", lambda: True)
        files = job._pka_output_files()
        for species in files:
            for path in files[species]["gas"] + files[species]["solv"]:
                _touch_output(path, "gaussian")

        job_result = job.print_thermochemistry()
        analyzed = compute_pka(
            ha_gas_file=files["HA"]["gas"],
            a_gas_file=files["A-"]["gas"],
            ha_solv_file=files["HA"]["solv"],
            a_solv_file=files["A-"]["solv"],
            scheme="direct",
        )
        first_only = compute_pka(
            ha_gas_file=files["HA"]["gas"][0],
            a_gas_file=files["A-"]["gas"][0],
            ha_solv_file=files["HA"]["solv"][0],
            a_solv_file=files["A-"]["solv"][0],
            scheme="direct",
        )

        table_path = temporary_working_dir / "outputs.csv"
        table_path.write_text("basename\nacid\n")
        table = PKaOutputTable.from_file(str(table_path))
        table.prepare(check_file_exists=True, scheme="direct")
        batch_results = table.run_pka(
            output_cls=compute_pka,
            scheme="direct",
        )

        assert job_result["num_conformers_HA"] == 3
        assert job_result["num_conformers_A"] == 3
        assert job_result["pKa"] == pytest.approx(analyzed["pKa"])
        assert batch_results[0]["pKa"] == pytest.approx(job_result["pKa"])
        assert job_result["pKa"] != pytest.approx(first_only["pKa"])

        thermo = job.compute_thermochemistry()
        assert len(thermo["HA"]) == 3
        assert len(thermo["A"]) == 3

        runner = CliRunner()
        cli_result = runner.invoke(
            run,
            [
                "pka",
                "-s",
                "direct",
                "analyze",
                "-ha",
                files["HA"]["gas"][0],
            ],
        )
        assert cli_result.exit_code == 0, cli_result.output
        assert f"{job_result['pKa']:.2f}" in cli_result.output


def _output_ext(backend):
    return "log" if backend == "gaussian" else "out"


class TestPkbConversion:
    """pKb = pKs − pKa conversion for analysis and summaries."""

    def test_pks_to_pkb_arithmetic(self):
        from chemsmart.cli.pka import pks_to_pkb

        assert pks_to_pkb(4.5, 14.0) == pytest.approx(9.5)
        assert pks_to_pkb(6.75, 16.7) == pytest.approx(9.95)
        assert pks_to_pkb(10.0, 14.0) == pytest.approx(4.0)

    def test_resolve_pkb_reporting_defaults(self):
        from chemsmart.cli.pka import DEFAULT_PKS, resolve_pkb_reporting

        assert resolve_pkb_reporting() == (False, None, False)
        assert resolve_pkb_reporting(pkb=False, pks=None) == (
            False,
            None,
            False,
        )
        assert resolve_pkb_reporting(pkb=True) == (True, DEFAULT_PKS, True)
        assert resolve_pkb_reporting(pkb=True, pks=16.7) == (True, 16.7, False)
        assert resolve_pkb_reporting(pkb=False, pks=16.7) == (
            True,
            16.7,
            False,
        )

    def test_compute_pka_does_not_return_pkb(self, tmp_path, monkeypatch):
        files = _build_outputs(tmp_path, "gaussian")
        _install_fake_thermochemistry(monkeypatch)
        from chemsmart.cli.pka import compute_pka

        result = compute_pka(
            ha_gas_file=files["ha.log"],
            a_gas_file=files["a.log"],
            href_gas_file=files["hb.log"],
            ref_gas_file=files["b.log"],
            ha_solv_file=files["has.log"],
            a_solv_file=files["as.log"],
            href_solv_file=files["hbs.log"],
            ref_solv_file=files["bs.log"],
            pka_reference=6.75,
        )
        assert "pKa" in result
        assert "pKb" not in result
        assert "pKs" not in result

    def test_warn_if_default_pks_non_aqueous(self, caplog):
        import logging

        from chemsmart.cli.pka import warn_if_default_pks_non_aqueous

        with caplog.at_level(logging.WARNING):
            warn_if_default_pks_non_aqueous(True, "acetonitrile")
        assert any(
            "pKs = 14.0" in rec.message and "acetonitrile" in rec.message
            for rec in caplog.records
        )

        caplog.clear()
        with caplog.at_level(logging.WARNING):
            warn_if_default_pks_non_aqueous(True, "water")
            warn_if_default_pks_non_aqueous(False, "acetonitrile")
            warn_if_default_pks_non_aqueous(True, None)
        assert not caplog.records

    def test_analyze_pkb_defaults_pks_to_14(self, tmp_path, monkeypatch):
        files = _build_outputs(tmp_path, "gaussian")
        _install_fake_thermochemistry(monkeypatch)
        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "--pkb",
                "analyze",
                "-ha",
                files["ha.log"],
                "-a",
                files["a.log"],
                "-hr",
                files["hb.log"],
                "-r",
                files["b.log"],
                "-has",
                files["has.log"],
                "-as",
                files["as.log"],
                "--href-solv",
                files["hbs.log"],
                "--ref-solv",
                files["bs.log"],
                "-rp",
                "6.75",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "pKs = 14.00 (default aqueous)" in result.output
        assert "Computed pKa(HA) = 6.75" in result.output
        assert "Computed pKb(B)  = 7.25" in result.output

    def test_analyze_pks_without_pkb_enables_reporting(
        self, tmp_path, monkeypatch
    ):
        files = _build_outputs(tmp_path, "gaussian")
        _install_fake_thermochemistry(monkeypatch)
        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "--pks",
                "16.7",
                "analyze",
                "-ha",
                files["ha.log"],
                "-a",
                files["a.log"],
                "-hr",
                files["hb.log"],
                "-r",
                files["b.log"],
                "-has",
                files["has.log"],
                "-as",
                files["as.log"],
                "--href-solv",
                files["hbs.log"],
                "--ref-solv",
                files["bs.log"],
                "-rp",
                "6.75",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "pKs = 16.70 (user-supplied)" in result.output
        assert "Computed pKa(HA) = 6.75" in result.output
        assert "Computed pKb(B)  = 9.95" in result.output

    def test_analyze_pkb_with_custom_pks(self, tmp_path, monkeypatch):
        files = _build_outputs(tmp_path, "gaussian")
        _install_fake_thermochemistry(monkeypatch)
        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "--pkb",
                "--pks",
                "16.7",
                "analyze",
                "-ha",
                files["ha.log"],
                "-a",
                files["a.log"],
                "-hr",
                files["hb.log"],
                "-r",
                files["b.log"],
                "-has",
                files["has.log"],
                "-as",
                files["as.log"],
                "--href-solv",
                files["hbs.log"],
                "--ref-solv",
                files["bs.log"],
                "-rp",
                "6.75",
            ],
        )
        assert result.exit_code == 0, result.output
        assert "pKs = 16.70 (user-supplied)" in result.output
        assert "Computed pKb(B)  = 9.95" in result.output

    def test_batch_analyze_pkb_adds_table_columns(self, tmp_path, monkeypatch):
        monkeypatch.chdir(tmp_path)
        basename = "target"
        for suffix in ("_pka_HA_opt", "_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"{basename}{suffix}.log").write_text(
                "Gaussian, Inc.\n"
            )
        for name in ("ref_HA_opt", "ref_A_opt", "ref_HA_sp", "ref_A_sp"):
            (tmp_path / f"{name}.log").write_text("Gaussian, Inc.\n")

        table = tmp_path / "pka_output.csv"
        table.write_text(
            "basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref\n"
            f"{basename},,,,,ref_HA_opt.log,ref_A_opt.log,ref_HA_sp.log,ref_A_sp.log,6.75\n"
        )
        _install_fake_thermochemistry(monkeypatch)

        runner = CliRunner()
        result = runner.invoke(
            run,
            ["pka", "--pkb", "batch-analyze", "-o", str(table)],
        )
        assert result.exit_code == 0, result.output
        assert "pKb" in result.output
        assert "pKs" in result.output
        assert "7.25" in result.output
        assert "14.00" in result.output

    def test_print_pka_summary_warns_default_pks_non_water(
        self, tmp_path, monkeypatch, caplog
    ):
        import logging

        files = _build_outputs(tmp_path, "gaussian")
        _install_fake_thermochemistry(monkeypatch)
        from chemsmart.cli.pka import print_pka_summary

        with caplog.at_level(logging.WARNING):
            print_pka_summary(
                ha_gas_file=files["ha.log"],
                a_gas_file=files["a.log"],
                href_gas_file=files["hb.log"],
                ref_gas_file=files["b.log"],
                ha_solv_file=files["has.log"],
                a_solv_file=files["as.log"],
                href_solv_file=files["hbs.log"],
                ref_solv_file=files["bs.log"],
                pka_reference=6.75,
                pkb=True,
                solvent_id="acetonitrile",
            )
        assert any(
            "pKs = 14.0" in rec.message and "acetonitrile" in rec.message
            for rec in caplog.records
        )


class TestPKa:
    """pKa CLI, batch submission, and job workflow tests."""

    def test_run_pka_detects_gaussian_and_dispatches(
        self, tmp_path, monkeypatch
    ):
        files = _build_outputs(tmp_path, "gaussian")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = _invoke_pka(runner, files)

        assert result.exit_code == 0
        assert "kwargs" in called
        assert called["kwargs"]["ha_gas_file"] == files["ha.log"]
        assert called["kwargs"]["a_solv_file"] == files["as.log"]
        assert called["kwargs"]["pka_reference"] == 6.75

    def test_run_pka_direct_analyze_dispatches(self, tmp_path, monkeypatch):
        files = _build_outputs(tmp_path, "gaussian")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = _invoke_pka_direct(runner, files, delta_g_proton=-270.0)

        assert result.exit_code == 0
        assert called["kwargs"]["ha_gas_file"] == files["ha.log"]
        assert called["kwargs"]["delta_G_proton"] == -270.0
        assert called["kwargs"]["scheme"] == "direct"

    def test_run_pka_direct_omitted_delta_g_computes_default(
        self, tmp_path, monkeypatch
    ):
        files = _build_outputs(tmp_path, "gaussian")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = _invoke_pka_direct(runner, files)

        assert result.exit_code == 0, result.output
        assert called["kwargs"]["delta_G_proton"] is None
        assert called["kwargs"]["scheme"] == "direct"

    def test_run_pka_detects_orca_and_dispatches(self, tmp_path, monkeypatch):
        files = _build_outputs(tmp_path, "orca")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = _invoke_pka(runner, files)

        assert result.exit_code == 0
        assert "kwargs" in called
        assert called["kwargs"]["href_gas_file"] == files["hb.log"]

    def test_run_pka_mixed_programs_analyze(self, tmp_path, monkeypatch):
        files = _build_outputs(tmp_path, "gaussian")
        _write_signature_file(Path(files["bs.log"]), "orca")
        called = {}

        def _fake_print(*args, **kwargs):
            called["kwargs"] = kwargs

        import chemsmart.cli.pka as pka_cli

        monkeypatch.setattr(pka_cli, "print_pka_summary", _fake_print)

        runner = CliRunner()
        result = _invoke_pka(runner, files)

        assert result.exit_code == 0, result.output
        assert called["kwargs"]["href_solv_file"] == files["hbs.log"]
        assert called["kwargs"]["ref_solv_file"] == files["bs.log"]

    def test_run_pka_batch_analyze_orca_outputs(self, tmp_path, monkeypatch):
        """batch-analyze should build Thermochemistry objects for ORCA files."""
        monkeypatch.chdir(tmp_path)
        basename = "target"
        for suffix in ("_pka_HA_opt", "_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"{basename}{suffix}.out").write_text(
                "* O   R   C   A *\n"
            )
        for name in ("ref_HA_opt", "ref_A_opt", "ref_HA_sp", "ref_A_sp"):
            (tmp_path / f"{name}.out").write_text("* O   R   C   A *\n")

        table = tmp_path / "pka_output.csv"
        table.write_text(
            "basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref\n"
            f"{basename},,,,,ref_HA_opt.out,ref_A_opt.out,ref_HA_sp.out,ref_A_sp.out,10.6\n"
        )

        constructed = _install_fake_thermochemistry(monkeypatch)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "-T",
                "333.15",
                "-csg",
                "100",
                "-ch",
                "100",
                "batch-analyze",
                "-o",
                str(table),
            ],
        )

        assert result.exit_code == 0, result.output
        assert "target_pka_HA_sp.out" in constructed
        assert "ref_HA_sp.out" in constructed
        assert "pKa" in result.output

    def test_run_pka_batch_analyze_mixed_gaussian_orca(
        self, tmp_path, monkeypatch
    ):
        """batch-analyze stays program-agnostic via Thermochemistry(filename=...)."""
        monkeypatch.chdir(tmp_path)
        basename = "target"
        for suffix in ("_pka_HA_opt", "_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"{basename}{suffix}.out").write_text(
                "* O   R   C   A *\n"
            )
        for name in ("ref_HA_opt", "ref_A_opt", "ref_HA_sp", "ref_A_sp"):
            (tmp_path / f"{name}.log").write_text("Gaussian, Inc.\n")

        table = tmp_path / "pka_output.csv"
        table.write_text(
            "basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref\n"
            f"{basename},,,,,ref_HA_opt.log,ref_A_opt.log,ref_HA_sp.log,ref_A_sp.log,10.6\n"
        )

        constructed = _install_fake_thermochemistry(monkeypatch)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "pka",
                "batch-analyze",
                "-o",
                str(table),
            ],
        )

        assert result.exit_code == 0, result.output
        assert "target_pka_HA_sp.out" in constructed
        assert "ref_HA_sp.log" in constructed
        assert "pKa" in result.output

    def test_run_pka_batch_analyze_gaussian_outputs(
        self, tmp_path, monkeypatch
    ):
        """batch-analyze should preserve Gaussian table behavior."""
        monkeypatch.chdir(tmp_path)
        basename = "target"
        for suffix in ("_pka_HA_opt", "_pka_A_opt", "_pka_HA_sp", "_pka_A_sp"):
            (tmp_path / f"{basename}{suffix}.log").write_text(
                "Gaussian, Inc.\n"
            )
        for name in ("ref_HA_opt", "ref_A_opt", "ref_HA_sp", "ref_A_sp"):
            (tmp_path / f"{name}.log").write_text("Gaussian, Inc.\n")

        table = tmp_path / "pka_output.csv"
        table.write_text(
            "basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref\n"
            f"{basename},,,,,ref_HA_opt.log,ref_A_opt.log,ref_HA_sp.log,ref_A_sp.log,6.75\n"
        )

        constructed = _install_fake_thermochemistry(monkeypatch)

        runner = CliRunner()
        result = runner.invoke(
            run,
            ["pka", "batch-analyze", "-o", str(table)],
        )

        assert result.exit_code == 0, result.output
        assert "target_pka_HA_opt.log" in constructed
        assert "ref_HA_sp.log" in constructed
        assert "pKa" in result.output

    def test_pka_thermochemistry_missing_scf_energy(
        self, tmp_path, monkeypatch
    ):
        class _MissingScfThermochemistry:
            electronic_energy = None
            qrrho_gibbs_free_energy = -1.0

            def __init__(self, filename, **kwargs):
                pass

        monkeypatch.setattr(
            "chemsmart.cli.pka.Thermochemistry",
            _MissingScfThermochemistry,
        )

        from chemsmart.cli.pka import pka_solvent_scf_energy

        with pytest.raises(ValueError, match="Could not extract SCF energy"):
            pka_solvent_scf_energy(str(tmp_path / "missing.out"))

    def test_pka_thermochemistry_missing_qh_gibbs(self, tmp_path, monkeypatch):
        class _MissingQhThermochemistry:
            electronic_energy = -1.0
            qrrho_gibbs_free_energy = None

            def __init__(self, filename, **kwargs):
                pass

        monkeypatch.setattr(
            "chemsmart.cli.pka.Thermochemistry",
            _MissingQhThermochemistry,
        )

        from chemsmart.cli.pka import pka_gas_phase_data

        with pytest.raises(
            ValueError,
            match="Could not extract quasi-harmonic Gibbs free energy",
        ):
            pka_gas_phase_data(str(tmp_path / "gas.out"))

    def test_run_pka_unparseable_output_raises(self, tmp_path):
        """analyze no longer pre-detects program type; parsing fails on bad files."""
        files = _build_outputs(tmp_path, "unknown")

        runner = CliRunner()
        result = _invoke_pka(runner, files)

        assert result.exit_code != 0

    def test_sub_orca_pka_batch_reconstructs_per_job_cli_args(
        self, tmp_path, monkeypatch, captured
    ):
        _require_backend_pka_subcommand(sub, "orca")
        orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.settings.server import Server

        acid1 = tmp_path / "acid1.xyz"
        acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        acid2 = tmp_path / "acid2.xyz"
        acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        table = tmp_path / "batch.xyz"
        table.write_text(
            "filepath proton_index charge multiplicity\n"
            f"{acid1} 2 0 1\n"
            f"{acid2} 2 1 2\n"
        )

        config_root = tmp_path / "chemsmart_cfg"
        orca_cfg_dir = config_root / "orca"
        orca_cfg_dir.mkdir(parents=True)
        (orca_cfg_dir / "test.yaml").write_text(
            "gas:\n"
            "  functional: B3LYP\n"
            "  basis: def2-SVP\n"
            "solv:\n"
            "  functional: B3LYP\n"
            "  basis: def2-SVP\n"
            "  freq: false\n"
            "  solvent_model: smd\n"
            "  solvent_id: water\n"
        )
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured["submissions"] = []

        fake_server = Server(name="dummy")
        real_from_filepath = Molecule.from_filepath

        def _fake_from_filepath(filepath, *args, **kwargs):
            if str(filepath) == str(table):
                placeholder = Molecule(
                    symbols=["C", "H"],
                    positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                    charge=0,
                    multiplicity=1,
                )
                if kwargs.get("return_list"):
                    return [placeholder]
                return placeholder
            return real_from_filepath(filepath, *args, **kwargs)

        def _fake_submit(job, test=False, cli_args=None, **kwargs):
            captured["submissions"].append((job, test, cli_args))

        monkeypatch.setattr(fake_server, "submit", _fake_submit)
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )
        monkeypatch.setattr(
            orca_cli.Molecule, "from_filepath", _fake_from_filepath
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--server",
                "dummy",
                "--test",
                "orca",
                "--project",
                "test",
                "--filename",
                str(table),
                "pka",
                "--scheme",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 2
        first_job, first_test, first_args = captured["submissions"][0]
        second_job, second_test, second_args = captured["submissions"][1]
        assert first_test is True
        assert second_test is True
        assert isinstance(first_args, list)
        assert isinstance(second_args, list)
        # Per-entry submit scripts should execute a single-row submission.
        assert "submit" in first_args
        assert "batch" not in first_args
        assert str(table) not in first_args

    def test_sub_orca_pka_batch_rewrites_per_entry_file_args(
        self, tmp_path, monkeypatch, captured
    ):
        _require_backend_pka_subcommand(sub, "orca")
        orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.settings.server import Server

        acid1 = tmp_path / "acid1.xyz"
        acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        acid2 = tmp_path / "acid2.xyz"
        acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        table = tmp_path / "batch.xyz"
        table.write_text(
            "filepath proton_index charge multiplicity\n"
            f"{acid1} 2 0 1\n"
            f"{acid2} 2 1 2\n"
        )

        config_root = tmp_path / "chemsmart_cfg"
        orca_cfg_dir = config_root / "orca"
        orca_cfg_dir.mkdir(parents=True)
        (orca_cfg_dir / "test.yaml").write_text(
            "gas:\n"
            "  functional: B3LYP\n"
            "  basis: def2-SVP\n"
            "solv:\n"
            "  functional: B3LYP\n"
            "  basis: def2-SVP\n"
            "  freq: false\n"
            "  solvent_model: smd\n"
            "  solvent_id: water\n"
        )
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured["submissions"] = []

        fake_server = Server(name="dummy")
        real_from_filepath = Molecule.from_filepath

        def _fake_from_filepath(filepath, *args, **kwargs):
            if str(filepath) == str(table):
                placeholder = Molecule(
                    symbols=["C", "H"],
                    positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                    charge=0,
                    multiplicity=1,
                )
                if kwargs.get("return_list"):
                    return [placeholder]
                return placeholder
            return real_from_filepath(filepath, *args, **kwargs)

        def _fake_submit(job, test=False, cli_args=None, **kwargs):
            captured["submissions"].append((job, test, cli_args))

        monkeypatch.setattr(fake_server, "submit", _fake_submit)
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )
        monkeypatch.setattr(
            orca_cli.Molecule, "from_filepath", _fake_from_filepath
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--server",
                "dummy",
                "--test",
                "orca",
                "--project",
                "test",
                "--filename",
                str(table),
                "pka",
                "--scheme",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 2

        first_args = captured["submissions"][0][2]
        second_args = captured["submissions"][1][2]

        assert str(table) not in first_args
        assert str(table) not in second_args
        assert str(acid1) in first_args
        assert str(acid2) in second_args
        # Rewritten row-level options must be placed in the correct command scope:
        # --charge/--multiplicity on backend command and --proton-index on pka.
        assert "--charge" in first_args
        assert "--multiplicity" in first_args
        assert "--proton-index" in first_args
        assert "submit" in first_args
        assert "batch" not in first_args
        assert first_args.index("--charge") < first_args.index("pka")
        assert first_args.index("--multiplicity") < first_args.index("pka")
        assert first_args.index("submit") < first_args.index("--proton-index")

    def test_sub_orca_pka_batch_shared_reference_loaded_once(
        self, tmp_path, monkeypatch, captured
    ):
        _require_backend_pka_subcommand(sub, "orca")
        orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings
        from chemsmart.settings.server import Server

        acid1 = tmp_path / "acid1.xyz"
        acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        acid2 = tmp_path / "acid2.xyz"
        acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        reference = tmp_path / "ref.xyz"
        reference.write_text("2\nref\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        table = tmp_path / "batch.xyz"
        table.write_text(
            "filepath proton_index charge multiplicity\n"
            f"{acid1} 2 0 1\n"
            f"{acid2} 2 1 2\n"
        )

        config_root = tmp_path / "chemsmart_cfg"
        orca_cfg_dir = config_root / "orca"
        orca_cfg_dir.mkdir(parents=True)
        (orca_cfg_dir / "test.yaml").write_text(
            "gas:\n"
            "  functional: B3LYP\n"
            "  basis: def2-SVP\n"
            "solv:\n"
            "  functional: B3LYP\n"
            "  basis: def2-SVP\n"
            "  freq: false\n"
            "  solvent_model: smd\n"
            "  solvent_id: water\n"
        )
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured["submissions"] = []
        reference_pair_call_count = {"count": 0}

        fake_server = Server(name="dummy")
        real_from_filepath = Molecule.from_filepath
        real_reference_pair_molecules = (
            ORCApKaJobSettings.reference_pair_molecules
        )

        def _fake_from_filepath(filepath, *args, **kwargs):
            if str(filepath) == str(table):
                placeholder = Molecule(
                    symbols=["C", "H"],
                    positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                    charge=0,
                    multiplicity=1,
                )
                if kwargs.get("return_list"):
                    return [placeholder]
                return placeholder
            return real_from_filepath(filepath, *args, **kwargs)

        def _counting_reference_pair(self):
            reference_pair_call_count["count"] += 1
            return real_reference_pair_molecules(self)

        def _fake_submit(job, test=False, cli_args=None, **kwargs):
            captured["submissions"].append((job, test, cli_args))

        monkeypatch.setattr(fake_server, "submit", _fake_submit)
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )
        monkeypatch.setattr(
            orca_cli.Molecule, "from_filepath", _fake_from_filepath
        )
        monkeypatch.setattr(
            ORCApKaJobSettings,
            "reference_pair_molecules",
            _counting_reference_pair,
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--server",
                "dummy",
                "--test",
                "orca",
                "--project",
                "test",
                "--filename",
                str(table),
                "pka",
                "--scheme",
                "proton exchange",
                "--reference",
                str(reference),
                "--reference-proton-index",
                "2",
                "--reference-charge",
                "0",
                "--reference-multiplicity",
                "1",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 2
        # Two batch jobs share one cached HRef/Ref- pair built at construction.
        assert reference_pair_call_count["count"] == 1

    def test_sub_orca_pka_batch_first_exchange_rest_direct(
        self, tmp_path, monkeypatch, captured
    ):
        _require_backend_pka_subcommand(sub, "orca")
        orca_cli = importlib.import_module("chemsmart.cli.orca.orca")

        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.settings.server import Server

        acid1 = tmp_path / "acid1.xyz"
        acid1.write_text("2\nacid1\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        acid2 = tmp_path / "acid2.xyz"
        acid2.write_text("2\nacid2\nN 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        reference = tmp_path / "ref.xyz"
        reference.write_text("2\nref\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")

        table = tmp_path / "batch.xyz"
        table.write_text(
            "filepath proton_index charge multiplicity\n"
            f"{acid1} 2 0 1\n"
            f"{acid2} 2 1 2\n"
        )

        config_root = tmp_path / "chemsmart_cfg"
        orca_cfg_dir = config_root / "orca"
        orca_cfg_dir.mkdir(parents=True)
        (orca_cfg_dir / "test.yaml").write_text(
            "gas:\n"
            "  functional: B3LYP\n"
            "  basis: def2-SVP\n"
            "solv:\n"
            "  functional: B3LYP\n"
            "  basis: def2-SVP\n"
            "  freq: false\n"
            "  solvent_model: smd\n"
            "  solvent_id: water\n"
        )
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured["submissions"] = []

        fake_server = Server(name="dummy")
        real_from_filepath = Molecule.from_filepath

        def _fake_from_filepath(filepath, *args, **kwargs):
            if str(filepath) == str(table):
                placeholder = Molecule(
                    symbols=["C", "H"],
                    positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
                    charge=0,
                    multiplicity=1,
                )
                if kwargs.get("return_list"):
                    return [placeholder]
                return placeholder
            return real_from_filepath(filepath, *args, **kwargs)

        def _fake_submit(job, test=False, cli_args=None, **kwargs):
            captured["submissions"].append((job, test, cli_args))

        monkeypatch.setattr(fake_server, "submit", _fake_submit)
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )
        monkeypatch.setattr(
            orca_cli.Molecule, "from_filepath", _fake_from_filepath
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--server",
                "dummy",
                "--test",
                "orca",
                "--project",
                "test",
                "--filename",
                str(table),
                "pka",
                "--scheme",
                "proton exchange",
                "--reference",
                str(reference),
                "--reference-proton-index",
                "2",
                "--reference-charge",
                "0",
                "--reference-multiplicity",
                "1",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 2

        first_job, _, first_args = captured["submissions"][0]
        second_job, _, second_args = captured["submissions"][1]

        assert first_job.settings.scheme == "proton exchange"
        assert second_job.settings.scheme == "direct"

        assert "--scheme" in first_args
        assert (
            first_args[first_args.index("--scheme") + 1] == "proton exchange"
        )
        assert "--reference" in first_args

        assert "--scheme" in second_args
        assert second_args[second_args.index("--scheme") + 1] == "direct"
        assert "--reference" not in second_args
        assert "--reference-proton-index" not in second_args
        assert "--reference-charge" not in second_args
        assert "--reference-multiplicity" not in second_args

    def test_run_gaussian_pka_help_is_submission_only(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        _require_backend_pka_subcommand(run, "gaussian")
        config_root = _write_test_backend_project(tmp_path, "gaussian")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "--no-scratch",
                "--fake",
                "gaussian",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "pka",
                "--help",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "\n  submit" in result.output
        assert "\n  batch" in result.output
        assert "\n  analyze" not in result.output
        assert "\n  thermo" not in result.output
        assert "\n  batch-analyze" not in result.output
        assert "--sampling" in result.output
        assert "--no-sampling" in result.output
        assert "--num-conformers" in result.output
        assert "-N" in result.output
        assert "--preview" in result.output
        _assert_ensemble_solvent_help(result.output)

    def test_run_orca_pka_help_is_submission_only(
        self, tmp_path, monkeypatch, single_molecule_xyz_file
    ):
        _require_backend_pka_subcommand(run, "orca")
        config_root = _write_test_backend_project(tmp_path, "orca")
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "--no-scratch",
                "--fake",
                "orca",
                "-p",
                "test",
                "-f",
                single_molecule_xyz_file,
                "pka",
                "--help",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "\n  submit" in result.output
        assert "\n  batch" in result.output
        assert "\n  analyze" not in result.output
        assert "\n  thermo" not in result.output
        assert "\n  batch-analyze" not in result.output
        assert "--sampling" in result.output
        assert "--no-sampling" in result.output
        assert "--num-conformers" in result.output
        assert "-N" in result.output
        assert "--preview" in result.output
        _assert_ensemble_solvent_help(result.output)

    def test_pka_docs_describe_matching_conformer_solvent_jobs(self):
        root = Path(__file__).resolve().parents[1]
        required = (
            "N solvent single-point calculations",
            "matching",
            "_c1",
            "--num-cores",
        )
        for relative in (
            "docs/source/pka-calculations.rst",
            "docs/source/gaussian-pka-calculations.rst",
            "docs/source/orca-pka-calculations.rst",
        ):
            text = (root / relative).read_text()
            _assert_no_ensemble_solvent_contradiction(text)
            for phrase in required:
                assert phrase in text, f"{relative} missing {phrase!r}"

    def test_run_pka_help_keeps_output_analysis_commands(self):
        runner = CliRunner()
        result = runner.invoke(
            run, ["--no-scratch", "--fake", "pka", "--help"]
        )

        assert result.exit_code == 0, result.output
        assert "analyze" in result.output
        assert "batch-analyze" in result.output
        assert "--pkb" in result.output
        assert "--pks" in result.output
        assert "--sampling" not in result.output
        assert "--num-conformers" not in result.output

        analyze = runner.invoke(
            run,
            ["--no-scratch", "--fake", "pka", "analyze", "--help"],
        )
        assert analyze.exit_code == 0, analyze.output
        analyze_help = " ".join(analyze.output.split())
        _assert_no_ensemble_solvent_contradiction(analyze_help)
        assert "matching pair" in analyze_help
        assert "_c1" in analyze_help

    def test_resolve_pka_sampling_options_rejects_n_without_sampling(self):
        from chemsmart.cli.pka import resolve_pka_sampling_options

        assert resolve_pka_sampling_options(False, 1) == (False, 1)
        assert resolve_pka_sampling_options(True, 3) == (True, 3)

        with pytest.raises(
            click.UsageError, match="-N/--num-conformers requires --sampling"
        ):
            resolve_pka_sampling_options(False, 3)

    def test_pka_settings_sampling_defaults_and_builders(self):
        from chemsmart.jobs.gaussian.settings import (
            GaussianJobSettings,
            GaussianpKaJobSettings,
        )
        from chemsmart.jobs.orca.settings import (
            ORCAJobSettings,
            ORCApKaJobSettings,
        )

        gaussian_defaults = GaussianpKaJobSettings(proton_index=10)
        assert gaussian_defaults.sampling is False
        assert gaussian_defaults.num_conformers == 1

        orca_defaults = ORCApKaJobSettings(proton_index=10)
        assert orca_defaults.sampling is False
        assert orca_defaults.num_conformers == 1

        shared = {
            "scheme": "direct",
            "reference": None,
            "reference_proton_index": None,
            "reference_charge": None,
            "reference_multiplicity": None,
            "reference_conjugate_base_charge": None,
            "reference_conjugate_base_multiplicity": None,
            "delta_g_proton": None,
            "conjugate_base_charge": None,
            "conjugate_base_multiplicity": None,
            "solvent_model": None,
            "solvent_id": None,
            "sampling": True,
            "num_conformers": 3,
            "temperature": 298.15,
            "concentration": 1.0,
            "pressure": 1.0,
            "cutoff_entropy_grimme": 100.0,
            "cutoff_enthalpy": 100.0,
            "skip_completed": True,
            "reference_color_code": 4,
        }
        gaussian = GaussianpKaJobSettings.build_gaussian_pka_settings(
            10,
            shared,
            GaussianJobSettings(functional="B3LYP", basis="6-31G*"),
        )
        orca = ORCApKaJobSettings.build_orca_pka_settings(
            10,
            shared,
            ORCAJobSettings(functional="B3LYP", basis="def2-SVP"),
        )
        for settings in (gaussian, orca):
            assert settings.sampling is True
            assert settings.num_conformers == 3
            assert settings.proton_index == 10

        with pytest.raises(ValueError, match="num_conformers must be >= 1"):
            ORCApKaJobSettings(proton_index=10, num_conformers=0)

    def test_validate_reference_options_requires_reference_for_proton_exchange(
        self, tmp_path
    ):
        from chemsmart.cli.pka import validate_reference_options

        with pytest.raises(click.UsageError, match="-r/--reference"):
            validate_reference_options(
                {"scheme": "proton exchange", "reference": None}
            )
        validate_reference_options({"scheme": "direct", "reference": None})

        ref = tmp_path / "ref.xyz"
        ref.write_text("2\nref\nO 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        shared = {
            "scheme": "proton exchange",
            "reference": str(ref),
            "reference_proton_index": None,
            "reference_color_code": None,
            "reference_charge": None,
            "reference_multiplicity": None,
        }
        with pytest.raises(
            click.UsageError, match="When --reference is provided"
        ):
            validate_reference_options(shared)

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_submit_proton_exchange_without_reference_fails(
        self, tmp_path, monkeypatch, backend
    ):
        """Default proton-exchange submit requires -r/--reference."""
        _require_backend_pka_subcommand(run, backend)
        acid = tmp_path / "acid.xyz"
        acid.write_text("2\nacid\nC 0.0 0.0 0.0\nH 0.0 0.0 1.0\n")
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "--no-scratch",
                "--fake",
                backend,
                "-p",
                "test",
                "-f",
                str(acid),
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-pi",
                "2",
                "submit",
            ],
        )

        assert result.exit_code != 0
        assert "-r/--reference" in result.output

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_auto_routes_to_batch(
        self, tmp_path, monkeypatch, backend, captured
    ):
        """Table -f input should use batch workflow without requiring -pi."""
        _require_backend_pka_subcommand(sub, backend)
        table, captured = _setup_sub_pka_batch_test(
            tmp_path, monkeypatch, backend, captured
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "proton-index is required" not in result.output
        assert len(captured["submissions"]) == 2

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_without_batch_subcommand(
        self, tmp_path, monkeypatch, backend, captured
    ):
        """Omitting the batch subcommand still routes table input to batch."""
        _require_backend_pka_subcommand(sub, backend)
        table, captured = _setup_sub_pka_batch_test(
            tmp_path, monkeypatch, backend, captured
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "proton-index is required" not in result.output
        assert len(captured["submissions"]) == 2

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_submit_subcommand_routes_to_batch(
        self, tmp_path, monkeypatch, backend, captured
    ):
        """Explicit submit with table -f still uses row-wise batch processing."""
        _require_backend_pka_subcommand(sub, backend)
        table, captured = _setup_sub_pka_batch_test(
            tmp_path, monkeypatch, backend, captured
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
                "submit",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "proton-index is required" not in result.output
        assert len(captured["submissions"]) == 2

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_batch_reconstructed_run_args_accept_proton_index(
        self, tmp_path, monkeypatch, backend, captured
    ):
        """Per-row chemsmart_run_*.py args must parse --proton-index under run."""
        _require_backend_pka_subcommand(sub, backend)
        table, captured = _setup_sub_pka_batch_test(
            tmp_path, monkeypatch, backend, captured
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )
        assert result.exit_code == 0, result.output
        assert captured["submissions"]

        cli_args = captured["submissions"][0][2]
        assert "submit" in cli_args
        assert "--proton-index" in cli_args
        assert cli_args.index("submit") < cli_args.index("--proton-index")

        from chemsmart.jobs.job import Job

        def _fake_run(self):
            return None

        monkeypatch.setattr(Job, "run", _fake_run)

        run_result = runner.invoke(run, ["--no-scratch", "--fake"] + cli_args)
        assert run_result.exit_code == 0, run_result.output
        assert "proton-index is required" not in run_result.output

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_cdxml_batch_uses_coloured_proton_fragments(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
        captured,
    ):
        """CDXML batch should create one job per coloured-proton fragment."""
        _require_backend_pka_subcommand(sub, backend)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured["labels"] = []
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "labels"
            ].append(job.label)
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                colored_proton_two_molecule_cdxml_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "Expected 5 fields" not in result.output
        assert "proton-index is required" not in result.output
        assert len(captured["labels"]) == 2
        assert all("_frag" in label for label in captured["labels"])

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_cdxml_batch_uses_molecule_charge_without_parent_flags(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
        captured,
    ):
        """CDXML batch should infer charge/mult from parsed Molecule objects."""
        _require_backend_pka_subcommand(sub, backend)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured["jobs"] = []
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "jobs"
            ].append(job)
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                colored_proton_two_molecule_cdxml_file,
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["jobs"]) == 2
        assert "ChemDraw molecular fragment" in result.output
        assert "Parent -c/--charge" not in result.output
        for job in captured["jobs"]:
            assert job.settings.charge == 0
            assert job.settings.multiplicity == 1
            assert job._batch_entry["charge"] == 0
            assert job._batch_entry["multiplicity"] == 1
            assert job._batch_entry["label"] == job.label

    def test_get_pka_molecules_auto_assigns_charge_and_multiplicity(
        self, colored_proton_cdxml_file
    ):
        from chemsmart.io.file import PKaCDXFile

        pka_mol = PKaCDXFile(colored_proton_cdxml_file).get_pka_molecules(
            index="-1"
        )
        assert pka_mol.charge == 0
        assert pka_mol.multiplicity == 1

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_cdxml_batch_reconstructed_scripts_target_single_fragment(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
        captured,
    ):
        """Each CDXML fragment script must submit only that fragment, not re-batch all."""
        _require_backend_pka_subcommand(sub, backend)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured["submissions"] = []
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "submissions"
            ].append((job, test, cli_args))
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                colored_proton_two_molecule_cdxml_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 2

        fragment_indices = []
        for job, _test, cli_args in captured["submissions"]:
            assert "batch" not in cli_args
            assert "submit" in cli_args
            assert "--proton-index" in cli_args
            assert "--index" in cli_args
            assert "--label" in cli_args
            assert cli_args[cli_args.index("--label") + 1] == job.label
            fragment_indices.append(cli_args[cli_args.index("--index") + 1])

        assert fragment_indices == ["1", "2"]

        from chemsmart.jobs.job import Job

        def _fake_run(self):
            return None

        monkeypatch.setattr(Job, "run", _fake_run)

        for job, _test, cli_args in captured["submissions"]:
            run_labels = []

            def _fake_run(self):
                run_labels.append(self.label)
                return None

            monkeypatch.setattr(Job, "run", _fake_run)
            run_result = runner.invoke(
                run, ["--no-scratch", "--fake"] + cli_args
            )
            assert run_result.exit_code == 0, run_result.output
            assert "proton-index is required" not in run_result.output
            assert run_labels == [job.label]

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_cdxml_batch_ignores_sibling_csv(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
        captured,
    ):
        """CDXML batch must not fall back to a sibling CSV submission table."""
        _require_backend_pka_subcommand(sub, backend)
        sibling_csv = Path(colored_proton_two_molecule_cdxml_file).with_suffix(
            ".csv"
        )
        sibling_csv.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            "only_one_row.xyz,1,0,1\n"
        )

        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured["labels"] = []
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "labels"
            ].append(job.label)
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                colored_proton_two_molecule_cdxml_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["labels"]) == 2
        assert all("_frag" in label for label in captured["labels"])

    def test_pka_resolve_proton_index_accepts_explicit_index(self):
        from chemsmart.cli.pka import resolve_proton_index

        proton_index, molecules = resolve_proton_index("acid.xyz", 8, None)
        assert proton_index == 8
        assert molecules is None

    def test_resolve_proton_index_uses_smarts_for_xyz_without_pi(
        self, tmp_path
    ):
        from chemsmart.cli.pka import (
            resolve_ionizable_site,
            resolve_proton_index,
        )

        mol = _molecule_from_smiles("c1ccccc1O")
        path = tmp_path / "phenol.xyz"
        mol.write(str(path), format="xyz")

        proton_index, molecules = resolve_proton_index(str(path), None, None)
        assert molecules is None
        assert mol.chemical_symbols[proton_index - 1] == "H"
        assert proton_index == resolve_ionizable_site(mol, mode="acid")

    def test_resolve_proton_index_uses_base_smarts_without_pi(self, tmp_path):
        from chemsmart.cli.pka import (
            resolve_ionizable_site,
            resolve_proton_index,
        )

        mol = _molecule_from_smiles("c1ccncc1")
        path = tmp_path / "pyridine.xyz"
        mol.write(str(path), format="xyz")

        site, molecules = resolve_proton_index(
            str(path), None, None, mode="base"
        )
        assert molecules is None
        assert mol.chemical_symbols[site - 1] == "N"
        assert site == resolve_ionizable_site(mol, mode="base")

    def test_resolve_proton_index_cdxml_pkb_colored_n(
        self, colored_basic_atom_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_proton_index
        from chemsmart.io.molecules.structure import Molecule

        site, molecules = resolve_proton_index(
            colored_basic_atom_cdxml_file, None, None, mode="base"
        )
        mol = Molecule.from_filepath(colored_basic_atom_cdxml_file)
        assert molecules is None
        assert mol.chemical_symbols[site - 1] == "N"

    def test_resolve_proton_index_cdxml_pkb_color_code(
        self, two_color_basic_atom_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_proton_index
        from chemsmart.io.molecules.structure import Molecule

        with pytest.raises(ValueError, match="Multiple uniquely coloured"):
            resolve_proton_index(
                two_color_basic_atom_cdxml_file, None, None, mode="base"
            )
        site, molecules = resolve_proton_index(
            two_color_basic_atom_cdxml_file, None, 4, mode="base"
        )
        mol = Molecule.from_filepath(two_color_basic_atom_cdxml_file)
        assert molecules is None
        assert mol.chemical_symbols[site - 1] == "N"

    def test_resolve_proton_index_cdxml_pkb_colored_h_errors(
        self, colored_proton_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_proton_index

        with pytest.raises(ValueError, match="colour the basic"):
            resolve_proton_index(
                colored_proton_cdxml_file, None, None, mode="base"
            )

    def test_resolve_proton_index_uncolored_cdxml_falls_through_to_pka_smarts(
        self, single_molecule_cdxml_file_benzene
    ):
        from chemsmart.cli.pka import resolve_proton_index

        with pytest.raises(ValueError, match="SMARTS matches"):
            resolve_proton_index(
                single_molecule_cdxml_file_benzene, None, None, mode="acid"
            )

    def test_resolve_proton_index_uncolored_cdxml_falls_through_to_pkb_smarts(
        self, uncolored_pyridine_cdxml_file
    ):
        from chemsmart.cli.pka import (
            resolve_ionizable_site,
            resolve_proton_index,
        )
        from chemsmart.io.molecules.structure import Molecule

        site, molecules = resolve_proton_index(
            uncolored_pyridine_cdxml_file, None, None, mode="base"
        )
        mol = Molecule.from_filepath(uncolored_pyridine_cdxml_file)
        assert molecules is None
        assert mol.chemical_symbols[site - 1] == "N"
        assert site == resolve_ionizable_site(mol, mode="base")

    def test_resolve_proton_index_cdxml_pkb_two_fragments(
        self, colored_basic_atom_two_molecule_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_proton_index

        site, molecules = resolve_proton_index(
            colored_basic_atom_two_molecule_cdxml_file,
            None,
            None,
            mode="base",
        )
        assert site is None
        assert len(molecules) == 2
        for mol in molecules:
            assert mol.chemical_symbols[mol.proton_index - 1] == "N"

    def test_resolve_proton_index_pka_colour_still_wins_over_smarts(
        self, colored_proton_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_proton_index

        proton_index, molecules = resolve_proton_index(
            colored_proton_cdxml_file, None, None, mode="acid"
        )
        assert molecules is None
        assert proton_index == 8

    def test_prepare_pkb_submit_molecule_protonates_and_increments_charge(
        self,
    ):
        from chemsmart.cli.pka import prepare_pkb_submit_molecule

        mol = _molecule_from_smiles("N")
        n_index = next(
            i + 1
            for i, symbol in enumerate(mol.chemical_symbols)
            if symbol == "N"
        )
        mol.charge = 0
        mol.multiplicity = 1
        settings = type("Settings", (), {"charge": 0, "multiplicity": 1})()

        protonated, proton_index, updated = prepare_pkb_submit_molecule(
            mol, n_index, settings
        )
        assert protonated.num_atoms == mol.num_atoms + 1
        assert protonated.chemical_symbols[proton_index - 1] == "H"
        assert proton_index == protonated.num_atoms
        assert protonated.charge == 1
        assert updated.charge == 1
        assert updated.multiplicity == 1

    def test_prepare_pkb_submit_molecule_rejects_hydrogen(self):
        from chemsmart.cli.pka import prepare_pkb_submit_molecule

        mol = _molecule_from_smiles("N")
        h_index = next(
            i + 1
            for i, symbol in enumerate(mol.chemical_symbols)
            if symbol == "H"
        )
        settings = type("Settings", (), {"charge": 0, "multiplicity": 1})()
        with pytest.raises(ValueError, match="heavy atom"):
            prepare_pkb_submit_molecule(mol, h_index, settings)

    def test_resolve_pka_batch_row_auto_detects_coloured_proton(
        self, colored_proton_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_pka_batch_row
        from chemsmart.io.molecules.structure import PKaMolecule

        proton_index, molecule = resolve_pka_batch_row(
            colored_proton_cdxml_file, proton_index=None
        )
        assert proton_index == 8
        assert isinstance(molecule, PKaMolecule)
        assert molecule.proton_index == 8

    def test_resolve_pka_batch_row_pkb_auto_detects_coloured_n(
        self, colored_basic_atom_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_pka_batch_row

        site, molecule = resolve_pka_batch_row(
            colored_basic_atom_cdxml_file, proton_index=None, pkb=True
        )
        assert molecule.chemical_symbols[site - 1] == "N"

    def test_resolve_pka_batch_row_explicit_index_overrides_cdxml(
        self, colored_proton_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_pka_batch_row
        from chemsmart.io.molecules.structure import Molecule

        proton_index, molecule = resolve_pka_batch_row(
            colored_proton_cdxml_file, proton_index=8
        )
        assert proton_index == 8
        assert isinstance(molecule, Molecule)

    def test_resolve_pka_batch_row_rejects_multi_molecule_cdxml(
        self, colored_proton_two_molecule_cdxml_file
    ):
        from chemsmart.cli.pka import resolve_pka_batch_row

        with pytest.raises(ValueError, match="single-molecule CDXML"):
            resolve_pka_batch_row(
                colored_proton_two_molecule_cdxml_file, proton_index=None
            )

    def test_resolve_pka_batch_row_uses_smarts_for_xyz_without_pi(
        self, tmp_path
    ):
        from chemsmart.cli.pka import (
            resolve_ionizable_site,
            resolve_pka_batch_row,
        )

        mol = _molecule_from_smiles("c1ccccc1O")
        path = tmp_path / "phenol.xyz"
        mol.write(str(path), format="xyz")

        proton_index, molecule = resolve_pka_batch_row(
            str(path), proton_index=None
        )
        assert molecule.chemical_symbols[proton_index - 1] == "H"
        assert proton_index == resolve_ionizable_site(mol, mode="acid")

    def test_resolve_pka_batch_row_smarts_errors_without_unique_xyz_site(
        self, single_molecule_xyz_file
    ):
        from chemsmart.cli.pka import resolve_pka_batch_row

        with pytest.raises(ValueError, match="SMARTS matches"):
            resolve_pka_batch_row(single_molecule_xyz_file, proton_index=None)

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_xyz_without_pi_uses_smarts(
        self, tmp_path, monkeypatch, backend
    ):
        mol = _molecule_from_smiles("c1ccccc1O")
        path = tmp_path / "phenol.xyz"
        mol.write(str(path), format="xyz")
        from chemsmart.cli.pka import resolve_ionizable_site

        expected = resolve_ionizable_site(mol, mode="acid")
        result, jobs = _capture_sub_pka_jobs(
            tmp_path,
            monkeypatch,
            backend,
            [
                "-f",
                str(path),
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "-s",
                "direct",
                "submit",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(jobs) == 1
        assert jobs[0].settings.proton_index == expected
        assert jobs[0].settings.charge == 0
        assert jobs[0].molecule.num_atoms == mol.num_atoms

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pkb_xyz_with_pi_protonates_before_job(
        self, tmp_path, monkeypatch, backend
    ):
        mol = _molecule_from_smiles("N")
        path = tmp_path / "ammonia.xyz"
        mol.write(str(path), format="xyz")
        n_index = next(
            i + 1
            for i, symbol in enumerate(mol.chemical_symbols)
            if symbol == "N"
        )
        result, jobs = _capture_sub_pka_jobs(
            tmp_path,
            monkeypatch,
            backend,
            [
                "-f",
                str(path),
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "--pkb",
                "-s",
                "direct",
                "-pi",
                str(n_index),
                "submit",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(jobs) == 1
        job = jobs[0]
        assert job.settings.pkb is True
        assert job.settings.charge == 1
        assert job.molecule.num_atoms == mol.num_atoms + 1
        assert job.settings.proton_index == job.molecule.num_atoms
        assert job.molecule.chemical_symbols[
            job.settings.proton_index - 1
        ] == ("H")
        assert job.molecule.charge == 1
        assert job.protonated_job.settings.charge == 1
        assert job.conjugate_base_job.settings.charge == 0

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pkb_cdxml_colored_n_without_pi(
        self, tmp_path, monkeypatch, backend, colored_basic_atom_cdxml_file
    ):
        from chemsmart.io.molecules.structure import Molecule

        mol = Molecule.from_filepath(colored_basic_atom_cdxml_file)
        n_index = next(
            i + 1
            for i, symbol in enumerate(mol.chemical_symbols)
            if symbol == "N"
        )
        result, jobs = _capture_sub_pka_jobs(
            tmp_path,
            monkeypatch,
            backend,
            [
                "-f",
                colored_basic_atom_cdxml_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "--pkb",
                "-s",
                "direct",
                "submit",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(jobs) == 1
        job = jobs[0]
        assert job.settings.pkb is True
        assert job.settings.charge == 1
        assert job.molecule.num_atoms == mol.num_atoms + 1
        assert job.settings.proton_index == job.molecule.num_atoms
        assert job.molecule.chemical_symbols[
            job.settings.proton_index - 1
        ] == ("H")
        assert mol.chemical_symbols[n_index - 1] == "N"

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pkb_uncolored_cdxml_falls_through_to_smarts(
        self, tmp_path, monkeypatch, backend, uncolored_pyridine_cdxml_file
    ):
        from chemsmart.io.molecules.structure import Molecule

        mol = Molecule.from_filepath(uncolored_pyridine_cdxml_file)
        result, jobs = _capture_sub_pka_jobs(
            tmp_path,
            monkeypatch,
            backend,
            [
                "-f",
                uncolored_pyridine_cdxml_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "--pkb",
                "-s",
                "direct",
                "submit",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(jobs) == 1
        job = jobs[0]
        assert job.settings.pkb is True
        assert job.settings.charge == 1
        assert job.molecule.num_atoms == mol.num_atoms + 1
        assert job.molecule.chemical_symbols[
            job.settings.proton_index - 1
        ] == ("H")

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pkb_cdxml_batch_per_fragment(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_basic_atom_two_molecule_cdxml_file,
    ):
        result, jobs = _capture_sub_pka_jobs(
            tmp_path,
            monkeypatch,
            backend,
            [
                "-f",
                colored_basic_atom_two_molecule_cdxml_file,
                "-c",
                "0",
                "-m",
                "1",
                "pka",
                "--pkb",
                "-s",
                "direct",
                "batch",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(jobs) == 2
        for job in jobs:
            assert job.settings.pkb is True
            assert job.settings.charge == 1
            assert job.molecule.chemical_symbols[
                job.settings.proton_index - 1
            ] == ("H")
            assert job._batch_entry["proton_index"] is not None
            assert (
                job.molecule.chemical_symbols[
                    job._batch_entry["proton_index"] - 1
                ]
                == "N"
            )
            assert job._batch_entry["charge"] == 0
            assert "_frag" in job.label

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pkb_batch_table_uses_basic_atom_index(
        self, tmp_path, monkeypatch, backend
    ):
        mol = _molecule_from_smiles("N")
        path = tmp_path / "ammonia.xyz"
        mol.write(str(path), format="xyz")
        n_index = next(
            i + 1
            for i, symbol in enumerate(mol.chemical_symbols)
            if symbol == "N"
        )
        table = tmp_path / "pkb.csv"
        table.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            f"{path},{n_index},0,1\n"
        )
        result, jobs = _capture_sub_pka_jobs(
            tmp_path,
            monkeypatch,
            backend,
            [
                "-f",
                str(table),
                "pka",
                "--pkb",
                "-s",
                "direct",
                "batch",
            ],
        )
        assert result.exit_code == 0, result.output
        assert len(jobs) == 1
        job = jobs[0]
        assert job.settings.charge == 1
        assert job.settings.proton_index == mol.num_atoms + 1
        assert job.molecule.chemical_symbols[
            job.settings.proton_index - 1
        ] == ("H")
        assert job._batch_entry["proton_index"] == n_index
        assert job._batch_entry["charge"] == 0

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_cdxml_blank_proton_index_auto_detects(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_cdxml_file,
        captured,
    ):
        """CDXML rows with blank proton_index auto-detect the coloured proton."""
        _require_backend_pka_subcommand(sub, backend)

        table = tmp_path / "pka_cdxml.csv"
        table.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            f"{colored_proton_cdxml_file},,0,1\n"
        )

        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured["submissions"] = []
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "submissions"
            ].append((job, test, cli_args))
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 1
        job = captured["submissions"][0][0]
        assert job.settings.proton_index == 8
        assert job.settings.delta_G_proton is None

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_cdxml_explicit_proton_index_overrides(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_cdxml_file,
        captured,
    ):
        """Explicit table proton_index overrides CDXML coloured-proton detection."""
        _require_backend_pka_subcommand(sub, backend)

        table = tmp_path / "pka_cdxml.csv"
        table.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            f"{colored_proton_cdxml_file},8,0,1\n"
        )

        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        captured["submissions"] = []
        fake_server.submit = (
            lambda job, test=False, cli_args=None, **kw: captured[
                "submissions"
            ].append((job, test, cli_args))
        )
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert len(captured["submissions"]) == 1
        job = captured["submissions"][0][0]
        assert job.settings.proton_index == 8
        assert job.settings.delta_G_proton is None

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_sub_pka_csv_table_rejects_multi_molecule_cdxml_row(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
    ):
        """Multi-molecule CDXML paths in a table row must fail clearly."""
        _require_backend_pka_subcommand(sub, backend)

        table = tmp_path / "pka_cdxml.csv"
        table.write_text(
            "filepath,proton_index,charge,multiplicity\n"
            f"{colored_proton_two_molecule_cdxml_file},,0,1\n"
        )

        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        from chemsmart.settings.server import Server

        fake_server = Server(name="dummy")
        fake_server.submit = lambda job, test=False, cli_args=None, **kw: None
        monkeypatch.setattr(
            "chemsmart.settings.server.Server.from_servername",
            lambda _name: fake_server,
        )

        runner = CliRunner()
        result = runner.invoke(
            sub,
            [
                "--test",
                "--server",
                "dummy",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code != 0
        assert "single-molecule CDXML" in result.output

    def test_orca_pka_job_generates_ha_and_a_subjobs(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch
    ):
        """ORCA pKa should prepare HA/A opt and SP jobs with Gaussian-style labels."""
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.opt import ORCAOptJob
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings
        from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )

        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )
        assert settings.delta_G_proton is None

        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="1a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )

        assert len(job.opt_jobs) == 2
        assert isinstance(job.protonated_job, ORCAOptJob)
        assert isinstance(job.conjugate_base_job, ORCAOptJob)
        assert job.protonated_job.label == "1a_pka_HA_opt"
        assert job.conjugate_base_job.label == "1a_pka_A_opt"
        assert job.conjugate_base_job.settings.charge == -1

        assert len(job.sp_jobs) == 2
        assert isinstance(job.protonated_sp_job, ORCASinglePointJob)
        assert isinstance(job.conjugate_base_sp_job, ORCASinglePointJob)
        assert job.protonated_sp_job.label == "1a_pka_HA_sp"
        assert job.conjugate_base_sp_job.label == "1a_pka_A_sp"

    def test_orca_pka_subjob_is_complete_uses_parent_folder(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch, tmp_path
    ):
        """Sub-jobs should detect completed outputs in the parent pKa folder."""
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )
        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )
        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="5a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )
        job.folder = str(tmp_path)

        for name in ("5a_pka_HA_opt", "5a_pka_A_opt"):
            (tmp_path / f"{name}.out").write_text(
                "****ORCA TERMINATED NORMALLY****\n"
            )

        assert all(j.is_complete() for j in job.opt_jobs)

    def test_orca_pka_run_sp_jobs_after_completed_opt(
        self,
        single_molecule_xyz_file,
        orca_jobrunner_no_scratch,
        tmp_path,
        monkeypatch,
        captured,
    ):
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )
        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )
        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="5a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )
        job.folder = str(tmp_path)

        for name in ("5a_pka_HA_opt", "5a_pka_A_opt"):
            (tmp_path / f"{name}.out").write_text(
                "****ORCA TERMINATED NORMALLY****\n"
            )

        captured["sp_labels"] = []

        def _fake_run_phase_jobs(*, jobs=None, jobs_factory=None, **kwargs):
            phase_jobs = jobs_factory() if jobs_factory is not None else jobs
            for child_job in phase_jobs:
                captured["sp_labels"].append(child_job.label)

        monkeypatch.setattr(
            "chemsmart.jobs.chain.pka.run_phase_jobs",
            _fake_run_phase_jobs,
        )
        monkeypatch.setattr(job, "_run_opt_jobs", lambda: None)
        monkeypatch.setattr(
            job, "_subjob_output", lambda *args, **kwargs: None
        )

        job._run()
        assert all(j.is_complete() for j in job.opt_jobs)
        assert captured["sp_labels"] == ["5a_pka_HA_sp", "5a_pka_A_sp"]

    def test_orca_pka_subjob_is_complete_recognizes_legacy_output(
        self, single_molecule_xyz_file, orca_jobrunner_no_scratch, tmp_path
    ):
        """Pre-rename ORCA pKa outputs should still count as complete."""
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )
        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )
        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="1a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )
        job.folder = str(tmp_path)

        legacy_out = tmp_path / "1a_pka.out"
        legacy_out.write_text("****ORCA TERMINATED NORMALLY****\n")

        assert job._subjob_is_complete(
            job.protonated_job, legacy_label="1a_pka"
        )

    def test_orca_pka_run_executes_ha_and_a_opt_jobs(
        self,
        single_molecule_xyz_file,
        orca_jobrunner_no_scratch,
        monkeypatch,
        captured,
    ):
        """ORCA pKa opt phase should run both acid and conjugate-base jobs."""
        from chemsmart.io.molecules.structure import Molecule
        from chemsmart.jobs.orca.pka import ORCApKaJob
        from chemsmart.jobs.orca.settings import ORCApKaJobSettings

        mol = Molecule.from_filepath(single_molecule_xyz_file)
        mol.charge = 0
        mol.multiplicity = 1
        proton_index = next(
            i + 1 for i, symbol in enumerate(mol.symbols) if symbol == "H"
        )

        settings = ORCApKaJobSettings(
            proton_index=proton_index,
            scheme="direct",
            functional="B3LYP",
            basis="def2-SVP",
        )
        job = ORCApKaJob(
            molecule=mol,
            settings=settings,
            label="1a_pka",
            jobrunner=orca_jobrunner_no_scratch,
        )

        captured["labels"] = []

        def _fake_run_phase_jobs(*, jobs, **kwargs):
            for child_job in jobs:
                captured["labels"].append(child_job.label)

        monkeypatch.setattr(
            "chemsmart.jobs.chain.pka.run_phase_jobs",
            _fake_run_phase_jobs,
        )

        job._run_opt_jobs()
        assert captured["labels"] == ["1a_pka_HA_opt", "1a_pka_A_opt"]

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_run_pka_batch_table_processing(
        self, tmp_path, monkeypatch, backend, captured
    ):
        """pKa table batch returns multiple jobs; run executes each locally."""
        _require_backend_pka_subcommand(run, backend)
        table = _build_pka_batch_table(tmp_path)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        captured["runs"] = []

        from chemsmart.jobs.job import Job

        def _fake_run(self):
            captured["runs"].append(self.label)

        monkeypatch.setattr(Job, "run", _fake_run)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "--no-scratch",
                "--fake",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "Batch job submission is not supported" not in result.output
        assert len(captured["runs"]) == 2
        if backend == "gaussian":
            assert set(captured["runs"]) == {"acid1", "acid2"}
        else:
            assert set(captured["runs"]) == {"acid1_pka", "acid2_pka"}

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_run_pka_batch_with_no_scratch(
        self, tmp_path, monkeypatch, backend, captured
    ):
        """Explicit --no-scratch should not require a scratch directory."""
        _require_backend_pka_subcommand(run, backend)
        table = _build_pka_batch_table(tmp_path)
        config_root = _write_test_backend_project(tmp_path, backend)
        monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))

        missing_scratch = tmp_path / "missing_scratch"
        from chemsmart.jobs import runner as runner_module

        monkeypatch.setattr(
            runner_module.user_settings, "scratch", str(missing_scratch)
        )

        captured["runs"] = []

        from chemsmart.jobs.job import Job

        def _fake_run(self):
            captured["runs"].append(self.label)

        monkeypatch.setattr(Job, "run", _fake_run)

        runner = CliRunner()
        result = runner.invoke(
            run,
            [
                "--fake",
                "--no-scratch",
                backend,
                "-p",
                "test",
                "-f",
                str(table),
                "pka",
                "-s",
                "direct",
                "batch",
            ],
        )

        assert result.exit_code == 0, result.output
        assert "Specified scratch dir does not exist" not in result.output
        assert len(captured["runs"]) == 2

    def test_run_rejects_non_job_batch_payload(self, pbs_server):
        """Scheduler-style batch payloads that are not Job lists stay blocked."""

        from chemsmart.cli.run import process_pipeline
        from chemsmart.jobs.runner import JobRunner

        ctx = click.Context(run)
        ctx.ensure_object(dict)
        ctx.obj["jobrunner"] = JobRunner(server=pbs_server, fake=True)

        with pytest.raises(
            ValueError, match="Batch job submission is not supported"
        ):
            process_pipeline.__wrapped__(ctx, ["not-a-job", "also-not-a-job"])


class TestIonizableSiteSMARTS:
    def test_acid_and_base_smarts_are_the_resolver_patterns(self):
        from chemsmart.cli.pka import (
            _IONIZABLE_SITE_SMARTS,
            PKA_ACID_SMARTS,
            PKB_BASE_SMARTS,
        )

        assert _IONIZABLE_SITE_SMARTS["acid"] is PKA_ACID_SMARTS
        assert _IONIZABLE_SITE_SMARTS["base"] is PKB_BASE_SMARTS
        assert "[CX3](=O)[OX2H1][#1]" in PKA_ACID_SMARTS
        assert "[c][OX2H1][#1]" in PKA_ACID_SMARTS
        assert "[SX2H1][#1]" in PKA_ACID_SMARTS
        assert "[NX4;+1][#1]" in PKA_ACID_SMARTS
        assert "[nH;+1][#1]" in PKA_ACID_SMARTS
        assert "[nX2;H0]" in PKB_BASE_SMARTS

    @pytest.mark.parametrize("smiles", ["c1ccccc1O", "CC(=O)O"])
    def test_acid_smarts_unique_hydrogen(self, smiles):
        from chemsmart.cli.pka import resolve_ionizable_site

        mol = _molecule_from_smiles(smiles)
        proton_index = resolve_ionizable_site(mol, mode="acid")
        assert mol.chemical_symbols[proton_index - 1] == "H"

    def test_acid_smarts_errors_when_carboxylic_and_phenol_both_present(self):
        from chemsmart.cli.pka import resolve_ionizable_site

        mol = _molecule_from_smiles("O=C(O)c1ccc(O)cc1")
        with pytest.raises(ValueError, match="2 SMARTS matches"):
            resolve_ionizable_site(mol, mode="acid")

    def test_acid_smarts_does_not_select_generic_alcohol(self):
        from chemsmart.cli.pka import resolve_ionizable_site

        mol = _molecule_from_smiles("CCO")
        with pytest.raises(ValueError, match="0 SMARTS matches"):
            resolve_ionizable_site(mol, mode="acid")

    def test_base_smarts_unique_pyridine_nitrogen(self):
        from chemsmart.cli.pka import resolve_ionizable_site

        mol = _molecule_from_smiles("c1ccncc1")
        site = resolve_ionizable_site(mol, mode="base")
        assert mol.chemical_symbols[site - 1] == "N"

    def test_base_smarts_errors_when_two_amines_present(self):
        from chemsmart.cli.pka import resolve_ionizable_site

        mol = _molecule_from_smiles("NCCN")
        with pytest.raises(ValueError, match="2 SMARTS matches"):
            resolve_ionizable_site(mol, mode="base")

    def test_base_smarts_does_not_select_amide_nitrogen(self):
        from chemsmart.cli.pka import resolve_ionizable_site

        mol = _molecule_from_smiles("CC(=O)N")
        with pytest.raises(ValueError, match="0 SMARTS matches"):
            resolve_ionizable_site(mol, mode="base")


def _forbid_pka_job_construction(monkeypatch):
    import importlib

    def _reject(*args, **kwargs):
        raise AssertionError("pKa job constructor called during preview")

    for module_name, attribute in (
        ("chemsmart.cli.gaussian.pka", "GaussianpKaJob"),
        ("chemsmart.jobs.gaussian.pka", "GaussianpKaJob"),
        ("chemsmart.jobs.orca.pka", "ORCApKaJob"),
        ("chemsmart.cli.pka", "build_pka_crest_job"),
    ):
        module = importlib.import_module(module_name)
        monkeypatch.setattr(module, attribute, _reject)


def _preflight_table(output):
    lines = []
    started = False
    for line in output.splitlines():
        if line.startswith("fragment  "):
            started = True
        if not started:
            continue
        if not line.strip():
            break
        lines.append(line)
    return "\n".join(lines)


def _invoke_pka_preview(tmp_path, monkeypatch, backend, filename, *pka_args):
    _require_backend_pka_subcommand(run, backend)
    config_root = _write_test_backend_project(tmp_path, backend)
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
    _forbid_pka_job_construction(monkeypatch)
    runner = CliRunner()
    args = [
        "--no-scratch",
        "--fake",
        backend,
        "-p",
        "test",
        "-f",
        str(filename),
        "-c",
        "0",
        "-m",
        "1",
        "pka",
        "--preview",
        "-s",
        "direct",
        *pka_args,
    ]
    return runner.invoke(run, args)


class TestPkaPreview:
    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_single_coloured_acid_preview(
        self, tmp_path, monkeypatch, backend, colored_proton_cdxml_file
    ):
        result = _invoke_pka_preview(
            tmp_path,
            monkeypatch,
            backend,
            colored_proton_cdxml_file,
            "batch",
        )
        assert result.exit_code == 0, result.output
        table = _preflight_table(result.output)
        assert table.startswith("fragment  label")
        assert "phenol_pka" in table
        assert "pKa" in table
        assert "8 H" in table
        assert "ChemDraw colour" in table
        assert "explicit" not in table
        assert "SMARTS" not in table
        assert "0" in table
        assert "1" in table
        assert table.endswith("ok")
        assert table.count("\n") == 2
        assert "ChemDraw molecular fragment" not in result.output

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_two_coloured_acid_fragments_preview(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
    ):
        from chemsmart.cli.pka import list_pka_preflight_sites

        result = _invoke_pka_preview(
            tmp_path,
            monkeypatch,
            backend,
            colored_proton_two_molecule_cdxml_file,
            "batch",
        )
        assert result.exit_code == 0, result.output
        table = _preflight_table(result.output)
        sites = list_pka_preflight_sites(
            colored_proton_two_molecule_cdxml_file, None, "acid"
        )
        assert len(sites) == 2
        rows = [line for line in table.splitlines() if line[:1].isdigit()]
        assert len(rows) == 2
        assert rows[0].startswith("1  ")
        assert rows[1].startswith("2  ")
        for number, (molecule, site, source) in enumerate(sites, start=1):
            element = molecule.chemical_symbols[site - 1]
            site_text = f"{site} {element}"
            assert source == "ChemDraw colour"
            assert f"phenol_two_molecule_frag{number}_pka" in rows[number - 1]
            assert site_text in rows[number - 1]
            assert "ChemDraw colour" in rows[number - 1]
            assert "pKa" in rows[number - 1]
        assert rows[0] != rows[1]
        assert "frag1_pka" in rows[0]
        assert "frag2_pka" in rows[1]
        assert "ChemDraw molecular fragment" in result.output
        assert (
            "Parent -c/--charge 0 and -m/--multiplicity 1 apply "
            "to every ChemDraw molecular fragment."
        ) in result.output
        assert "Salts, counterions" in result.output

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_two_coloured_base_fragments_preview(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_basic_atom_two_molecule_cdxml_file,
    ):
        result = _invoke_pka_preview(
            tmp_path,
            monkeypatch,
            backend,
            colored_basic_atom_two_molecule_cdxml_file,
            "--pkb",
            "batch",
        )
        assert result.exit_code == 0, result.output
        table = _preflight_table(result.output)
        rows = [line for line in table.splitlines() if line[:1].isdigit()]
        assert len(rows) == 2
        for number, row in enumerate(rows, start=1):
            assert row.startswith(f"{number}  ")
            assert f"pyridine_two_molecule_frag{number}_pka" in row
            assert "pKb" in row
            assert " N" in row
            assert "ChemDraw colour" in row
            assert "added H" in row

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_uncoloured_smarts_preview(
        self, tmp_path, monkeypatch, backend, uncolored_pyridine_cdxml_file
    ):
        result = _invoke_pka_preview(
            tmp_path,
            monkeypatch,
            backend,
            uncolored_pyridine_cdxml_file,
            "--pkb",
            "batch",
        )
        assert result.exit_code == 0, result.output
        table = _preflight_table(result.output)
        assert "SMARTS" in table
        assert "ChemDraw colour" not in table
        assert "explicit" not in table
        assert "pKb" in table
        assert " N" in table
        assert "pyridine_uncolored_pka" in table
        assert "added H" in table

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_explicit_index_overrides_colour(
        self, tmp_path, monkeypatch, backend, colored_proton_cdxml_file
    ):
        from chemsmart.io.molecules.structure import Molecule

        molecule = Molecule.from_filepath(colored_proton_cdxml_file)
        explicit = next(
            index
            for index, symbol in enumerate(molecule.chemical_symbols, start=1)
            if symbol == "H" and index != 8
        )
        result = _invoke_pka_preview(
            tmp_path,
            monkeypatch,
            backend,
            colored_proton_cdxml_file,
            "-pi",
            str(explicit),
            "batch",
        )
        assert result.exit_code == 0, result.output
        table = _preflight_table(result.output)
        assert f"{explicit} H" in table
        assert "explicit" in table
        assert "ChemDraw colour" not in table
        assert "SMARTS" not in table
        assert "8 H" not in table

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_ambiguous_colour_fails_before_jobs(
        self, tmp_path, monkeypatch, backend, two_color_basic_atom_cdxml_file
    ):
        result = _invoke_pka_preview(
            tmp_path,
            monkeypatch,
            backend,
            two_color_basic_atom_cdxml_file,
            "--pkb",
            "batch",
        )
        assert result.exit_code != 0
        assert "Multiple uniquely coloured" in result.output
        assert "constructor called" not in result.output
        assert _preflight_table(result.output) == ""

    def test_multiple_smarts_sites_fail_before_jobs(
        self, tmp_path, monkeypatch
    ):
        molecule = _molecule_from_smiles("Oc1ccc(O)cc1")
        path = tmp_path / "hydroquinone.xyz"
        molecule.write(str(path), format="xyz")
        result = _invoke_pka_preview(
            tmp_path, monkeypatch, "gaussian", path, "batch"
        )
        assert result.exit_code != 0
        assert "SMARTS matches" in result.output
        assert "constructor called" not in result.output

    def test_preview_without_batch_does_not_create_jobs(
        self, tmp_path, monkeypatch, colored_proton_cdxml_file
    ):
        result = _invoke_pka_preview(
            tmp_path,
            monkeypatch,
            "gaussian",
            colored_proton_cdxml_file,
        )
        assert result.exit_code == 0, result.output
        table = _preflight_table(result.output)
        assert "phenol_pka" in table
        assert "8 H" in table
        assert "constructor called" not in result.output


def _invoke_pka_without_charge(
    tmp_path, monkeypatch, backend, filename, *pka_args
):
    _require_backend_pka_subcommand(run, backend)
    config_root = _write_test_backend_project(tmp_path, backend)
    monkeypatch.setenv("CHEMSMART_CONFIG_DIR", str(config_root))
    _forbid_pka_job_construction(monkeypatch)
    runner = CliRunner()
    args = [
        "--no-scratch",
        "--fake",
        backend,
        "-p",
        "test",
        "-f",
        str(filename),
        "pka",
        "-s",
        "direct",
        *pka_args,
    ]
    return runner.invoke(run, args)


class TestPkaSeminarValidation:
    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_missing_charge_and_multiplicity_fail_before_preview(
        self, tmp_path, monkeypatch, backend
    ):
        molecule = _molecule_from_smiles("Oc1ccccc1")
        path = tmp_path / "phenol.xyz"
        molecule.write(str(path), format="xyz")
        preview = _invoke_pka_without_charge(
            tmp_path, monkeypatch, backend, path, "--preview"
        )
        assert preview.exit_code != 0
        assert "Fragment 1" in preview.output
        assert "-c/--charge" in preview.output
        assert "-m/--multiplicity" in preview.output
        assert "constructor called" not in preview.output

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_missing_charge_and_multiplicity_fail_before_submit(
        self, tmp_path, monkeypatch, backend
    ):
        molecule = _molecule_from_smiles("Oc1ccccc1")
        path = tmp_path / "phenol.xyz"
        molecule.write(str(path), format="xyz")
        submitted = _invoke_pka_without_charge(
            tmp_path, monkeypatch, backend, path
        )
        assert submitted.exit_code != 0
        assert "Charge and multiplicity are required" in submitted.output
        assert "-c/--charge" in submitted.output
        assert "-m/--multiplicity" in submitted.output
        assert "constructor called" not in submitted.output

    @pytest.mark.parametrize("backend", ["gaussian", "orca"])
    def test_multi_fragment_preview_without_parent_flags_keeps_parsed_charge(
        self,
        tmp_path,
        monkeypatch,
        backend,
        colored_proton_two_molecule_cdxml_file,
    ):
        result = _invoke_pka_without_charge(
            tmp_path,
            monkeypatch,
            backend,
            colored_proton_two_molecule_cdxml_file,
            "--preview",
            "batch",
        )
        assert result.exit_code == 0, result.output
        assert "ChemDraw molecular fragment" in result.output
        assert "Parent -c/--charge" not in result.output
        assert "Parent -m/--multiplicity" not in result.output
        table = _preflight_table(result.output)
        rows = [line for line in table.splitlines() if line[:1].isdigit()]
        assert len(rows) == 2
        assert "constructor called" not in result.output

    def test_output_errors_name_species_and_conformer(self, monkeypatch):
        from chemsmart.cli.pka import _species_solution_free_energy

        def _fail_gas(filepath, **kwargs):
            if str(filepath).endswith("_c2.log"):
                raise ValueError(
                    f"File '{filepath}' did not terminate normally. "
                    "Skipping thermochemistry calculation for this file."
                )
            return -1.0, 0.01

        monkeypatch.setattr("chemsmart.cli.pka.pka_gas_phase_data", _fail_gas)
        monkeypatch.setattr(
            "chemsmart.cli.pka.pka_solvent_scf_energy",
            lambda filepath, **kwargs: -1.1,
        )
        with pytest.raises(ValueError, match=r"HA c2: File '.*_c2.log'"):
            _species_solution_free_energy(
                ["ha_opt_c1.log", "ha_opt_c2.log"],
                ["ha_sp_c1.log", "ha_sp_c2.log"],
                {},
                298.15,
                "HA",
            )

        def _imaginary(filepath, **kwargs):
            raise ValueError(
                f"Invalid geometry optimization for {filepath}. "
                "A valid optimized geometry should not contain "
                "imaginary frequencies."
            )

        monkeypatch.setattr("chemsmart.cli.pka.pka_gas_phase_data", _imaginary)
        with pytest.raises(
            ValueError, match=r"A- conformer 1: Invalid geometry"
        ):
            _species_solution_free_energy(
                ["a_opt.log"],
                ["a_sp.log"],
                {},
                298.15,
                "A-",
            )

    def test_inconsistent_conformer_temperatures_warn(self, caplog):
        import logging

        from chemsmart.cli.pka import (
            _warn_inconsistent_conformer_temperatures,
        )

        with caplog.at_level(logging.WARNING):
            _warn_inconsistent_conformer_temperatures(
                "HA",
                [
                    ("ha_opt_c1.log", 298.15),
                    ("ha_opt_c2.log", 310.0),
                ],
            )
        assert "HA conformer outputs use inconsistent temperatures" in (
            caplog.text
        )
        assert "ha_opt_c1.log=298.15 K" in caplog.text
        assert "ha_opt_c2.log=310 K" in caplog.text

        caplog.clear()
        with caplog.at_level(logging.WARNING):
            _warn_inconsistent_conformer_temperatures(
                "HA",
                [
                    ("ha_opt_c1.log", 298.15),
                    ("ha_opt_c2.log", 298.15),
                ],
            )
        assert "inconsistent temperatures" not in caplog.text
