import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.job import GaussianComJob
from chemsmart.jobs.gaussian.reaction import (
    build_qst_input_string,
    make_qst_com_job,
)
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.writer import GaussianInputWriter


def _h2(z_h2=0.74, charge=0, multiplicity=1):
    return Molecule(
        symbols=["H", "H"],
        positions=[[0.0, 0.0, 0.0], [0.0, 0.0, z_h2]],
        charge=charge,
        multiplicity=multiplicity,
    )


def _water():
    return Molecule(
        symbols=["O", "H", "H"],
        positions=[
            [0.0, 0.0, 0.0],
            [0.96, 0.0, 0.0],
            [-0.24, 0.93, 0.0],
        ],
        charge=0,
        multiplicity=1,
    )


def _qst_settings(**kwargs):
    values = dict(
        functional="b3lyp",
        basis="def2svp",
        charge=0,
        multiplicity=1,
        freq=False,
        chk=True,
        title="QST job",
    )
    values.update(kwargs)
    return GaussianJobSettings(**values)


class TestGaussianQSTInput:
    def test_qst2_builds_two_coordinate_blocks(
        self, gaussian_jobrunner_no_scratch
    ):
        text = build_qst_input_string(
            reactant=_h2(0.74),
            product=_h2(0.90),
            settings=_qst_settings(),
            label="rxn_qst",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        lowered = text.lower()
        assert "opt=qst2" in lowered
        assert "opt=qst3" not in lowered
        assert text.count("QST job: reactant") == 1
        assert text.count("QST job: product") == 1
        assert "QST job: ts guess" not in text
        assert text.count("\n0 1\n") == 2
        assert "0.7400000000" in text
        assert "0.9000000000" in text
        assert "%chk=rxn_qst.chk" in text
        assert (
            f"%nprocshared={gaussian_jobrunner_no_scratch.num_cores}" in text
        )
        assert f"%mem={gaussian_jobrunner_no_scratch.mem_gb}GB" in text

    def test_qst3_builds_three_coordinate_blocks(
        self, gaussian_jobrunner_no_scratch
    ):
        text = build_qst_input_string(
            reactant=_h2(0.74),
            product=_h2(0.90),
            settings=_qst_settings(),
            ts_guess=_h2(0.82),
            label="rxn_qst",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        lowered = text.lower()
        assert "opt=qst3" in lowered
        assert "opt=qst2" not in lowered
        assert text.count("QST job: reactant") == 1
        assert text.count("QST job: product") == 1
        assert text.count("QST job: ts guess") == 1
        assert text.count("\n0 1\n") == 3
        assert "0.7400000000" in text
        assert "0.9000000000" in text
        assert "0.8200000000" in text

    def test_atom_count_mismatch_errors(self):
        with pytest.raises(ValueError, match="same number of atoms"):
            build_qst_input_string(
                reactant=_h2(),
                product=_water(),
                settings=_qst_settings(),
            )

    def test_ts_guess_atom_count_mismatch_errors(self):
        with pytest.raises(ValueError, match="same number of atoms"):
            build_qst_input_string(
                reactant=_h2(),
                product=_h2(0.90),
                settings=_qst_settings(),
                ts_guess=_water(),
            )

    def test_atom_order_mismatch_errors(self):
        swapped = Molecule(
            symbols=["H", "C"],
            positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]],
            charge=0,
            multiplicity=1,
        )
        reactant = Molecule(
            symbols=["C", "H"],
            positions=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.1]],
            charge=0,
            multiplicity=1,
        )
        with pytest.raises(ValueError, match="same atom order"):
            build_qst_input_string(
                reactant=reactant,
                product=swapped,
                settings=_qst_settings(),
            )

    def test_qst_route_reuses_existing_ts_route(
        self, gaussian_jobrunner_no_scratch
    ):
        settings = _qst_settings(
            jobtype="ts",
            freq=True,
            additional_route_parameters="empiricaldispersion=gd3bj",
        )
        original = settings.route_string.lower()
        assert "opt=(ts,calcfc,noeigentest)" in original
        text = build_qst_input_string(
            reactant=_h2(),
            product=_h2(0.90),
            settings=settings,
            jobrunner=gaussian_jobrunner_no_scratch,
        ).lower()
        assert "opt=(qst2,calcfc,noeigentest)" in text
        assert "opt=(ts," not in text
        assert "empiricaldispersion=gd3bj" in text
        assert " freq " in text

    def test_additional_opt_options_go_in_qst_opt_keyword(
        self, gaussian_jobrunner_no_scratch
    ):
        text = build_qst_input_string(
            reactant=_h2(),
            product=_h2(0.90),
            settings=_qst_settings(
                jobtype="opt",
                additional_opt_options_in_route="calcfc,noeigentest",
            ),
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert "opt=(qst2,calcfc,noeigentest)" in text.lower()

    def test_make_qst_com_job_sets_input_string(
        self, gaussian_jobrunner_no_scratch
    ):
        job = make_qst_com_job(
            reactant=_h2(),
            product=_h2(0.90),
            settings=_qst_settings(),
            label="sn2_qst",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert isinstance(job, GaussianComJob)
        assert job.TYPE == "g16com"
        assert job.settings.input_string
        assert "opt=qst2" in job.settings.input_string.lower()
        assert "%chk=sn2_qst.chk" in job.settings.input_string

    def test_writer_writes_qst_input_string_as_com(
        self, tmpdir, gaussian_jobrunner_no_scratch
    ):
        job = make_qst_com_job(
            reactant=_h2(),
            product=_h2(0.90),
            ts_guess=_h2(0.82),
            settings=_qst_settings(),
            label="rxn_qst",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        GaussianInputWriter(job=job).write(target_directory=str(tmpdir))
        written = tmpdir.join("rxn_qst.com").read()
        assert written == job.settings.input_string
        assert "opt=qst3" in written.lower()
        assert written.count("\n0 1\n") == 3
