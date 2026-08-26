import pytest

from chemsmart.io.molecules.structure import Molecule
from chemsmart.jobs.gaussian.job import GaussianComJob, GaussianJob
from chemsmart.jobs.gaussian.opt import GaussianOptJob
from chemsmart.jobs.gaussian.reaction import (
    GaussianReactionJob,
    build_qst_input_string,
    make_qst_com_job,
)
from chemsmart.jobs.gaussian.settings import GaussianJobSettings
from chemsmart.jobs.gaussian.singlepoint import GaussianSinglePointJob
from chemsmart.jobs.gaussian.ts import GaussianTSJob
from chemsmart.jobs.gaussian.writer import GaussianInputWriter
from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.neb import ORCANEBJob
from chemsmart.jobs.orca.opt import ORCAOptJob
from chemsmart.jobs.orca.reaction import ORCAReactionJob
from chemsmart.jobs.orca.settings import ORCAJobSettings
from chemsmart.jobs.orca.singlepoint import ORCASinglePointJob
from chemsmart.jobs.orca.ts import ORCATSJob


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


def _opt_settings(**kwargs):
    values = dict(
        functional="b3lyp",
        basis="def2svp",
        charge=0,
        multiplicity=1,
        jobtype="opt",
        freq=True,
        title="opt",
    )
    values.update(kwargs)
    return GaussianJobSettings(**values)


def _ts_settings(**kwargs):
    values = dict(
        functional="b3lyp",
        basis="def2svp",
        charge=0,
        multiplicity=1,
        jobtype="ts",
        freq=True,
        title="ts",
    )
    values.update(kwargs)
    return GaussianJobSettings(**values)


def _sp_settings(**kwargs):
    values = dict(
        functional="wb97xd",
        basis="def2tzvp",
        charge=0,
        multiplicity=1,
        jobtype="sp",
        freq=False,
        solvent_model="smd",
        solvent_id="water",
        title="sp",
    )
    values.update(kwargs)
    return GaussianJobSettings(**values)


def _orca_settings(jobtype="opt", **kwargs):
    values = dict(
        functional="b3lyp",
        basis="def2-svp",
        charge=0,
        multiplicity=1,
        jobtype=jobtype,
        freq=jobtype != "sp",
        title=jobtype,
    )
    values.update(kwargs)
    return ORCAJobSettings(**values)


class TestGaussianReactionJob:
    def test_case1_skips_guess_and_builds_ts_opt(
        self, gaussian_jobrunner_no_scratch
    ):
        job = GaussianReactionJob(
            molecule=_h2(0.82),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
            opt_settings=_opt_settings(),
            ts_settings=_ts_settings(),
            sp_settings=_sp_settings(),
        )
        assert isinstance(job, GaussianJob)
        assert job.TYPE == "g16reaction"
        assert [phase.name for phase in job.phases] == ["Guess", "Opt", "SP"]
        assert job.phase_by_name("Guess").should_skip()
        assert job.guess_job is None
        assert isinstance(job.ts_opt_job, GaussianTSJob)
        assert job.ts_opt_job.label == "sn2_TS_opt"
        assert job.ts_opt_job.settings.jobtype == "ts"
        assert job.ts_opt_job.settings.freq is True
        assert job.reactant_opt_jobs == []
        assert job.product_opt_jobs == []
        assert job.opt_jobs == [job.ts_opt_job]

    def test_case1_optional_reactant_and_product_are_opt_jobs(
        self, gaussian_jobrunner_no_scratch
    ):
        job = GaussianReactionJob(
            molecule=_h2(0.82),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
            reactants=[_h2(0.74, charge=-1)],
            products=[_h2(0.90, charge=-1)],
            no_path_search=True,
            opt_settings=_opt_settings(),
            ts_settings=_ts_settings(),
            sp_settings=_sp_settings(),
        )
        assert job.phase_by_name("Guess").should_skip()
        assert job.guess_job is None
        assert isinstance(job.reactant_opt_jobs[0], GaussianOptJob)
        assert job.reactant_opt_jobs[0].label == "sn2_R_opt"
        assert job.reactant_opt_jobs[0].settings.charge == -1
        assert isinstance(job.ts_opt_job, GaussianTSJob)
        assert isinstance(job.product_opt_jobs[0], GaussianOptJob)
        assert job.product_opt_jobs[0].label == "sn2_P_opt"
        assert [child.label for child in job.opt_jobs] == [
            "sn2_R_opt",
            "sn2_TS_opt",
            "sn2_P_opt",
        ]

    def test_case2_uses_qst_com_job_then_ts_job(
        self, gaussian_jobrunner_no_scratch
    ):
        job = GaussianReactionJob(
            molecule=_h2(0.74),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
            products=[_h2(0.90)],
            opt_settings=_opt_settings(),
            ts_settings=_ts_settings(),
            sp_settings=_sp_settings(),
        )
        assert not job.phase_by_name("Guess").should_skip()
        assert isinstance(job.guess_job, GaussianComJob)
        assert job.guess_job.TYPE == "g16com"
        assert job.guess_job.label == "sn2_qst"
        qst_text = job.guess_job.settings.input_string.lower()
        assert "qst2" in qst_text
        assert "qst3" not in qst_text
        assert isinstance(job.ts_opt_job, GaussianTSJob)
        assert job.reactant_opt_jobs[0].label == "sn2_R_opt"
        assert job.product_opt_jobs[0].label == "sn2_P_opt"

    def test_case2_qst3_when_ts_guess_given(
        self, gaussian_jobrunner_no_scratch
    ):
        job = GaussianReactionJob(
            molecule=_h2(0.74),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
            products=[_h2(0.90)],
            ts_guess=_h2(0.82),
            opt_settings=_opt_settings(),
            ts_settings=_ts_settings(),
        )
        qst_text = job.guess_job.settings.input_string.lower()
        assert "qst3" in qst_text
        assert "qst2" not in qst_text

    def test_case2_reactant_and_product_use_parent_as_qst3_guess(
        self, gaussian_jobrunner_no_scratch
    ):
        job = GaussianReactionJob(
            molecule=_h2(0.82),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
            reactants=[_h2(0.74)],
            products=[_h2(0.90)],
        )
        qst_text = job.guess_job.settings.input_string.lower()
        assert "qst3" in qst_text
        assert job.reactant_opt_jobs[0].label == "sn2_R_opt"
        assert job.opt_jobs[0].molecule.positions[1][2] == pytest.approx(0.74)

    def test_no_path_search_skips_guess_when_r_p_ts_set(
        self, gaussian_jobrunner_no_scratch
    ):
        job = GaussianReactionJob(
            molecule=_h2(0.82),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
            reactants=[_h2(0.74)],
            products=[_h2(0.90)],
            no_path_search=True,
        )
        assert job.phase_by_name("Guess").should_skip()
        assert job.guess_job is None
        assert isinstance(job.ts_opt_job, GaussianTSJob)

    def test_multiple_reactants_use_numbered_labels(
        self, gaussian_jobrunner_no_scratch
    ):
        job = GaussianReactionJob(
            molecule=_h2(0.82),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
            reactants=[_h2(0.74, charge=-1), _h2(0.80)],
            no_path_search=True,
        )
        assert [child.label for child in job.reactant_opt_jobs] == [
            "sn2_R1_opt",
            "sn2_R2_opt",
        ]

    def test_sp_uses_solv_theory_and_opt_labels(
        self, gaussian_jobrunner_no_scratch
    ):
        job = GaussianReactionJob(
            molecule=_h2(0.82),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
            reactants=[_h2(0.74)],
            products=[_h2(0.90)],
            no_path_search=True,
            opt_settings=_opt_settings(),
            ts_settings=_ts_settings(),
            sp_settings=_sp_settings(),
        )
        sp_jobs = job._sp_jobs_for_phase()
        assert [child.label for child in sp_jobs] == [
            "sn2_R_sp",
            "sn2_TS_sp",
            "sn2_P_sp",
        ]
        assert all(
            isinstance(child, GaussianSinglePointJob) for child in sp_jobs
        )
        assert job.ts_sp_job.settings.functional == "wb97xd"
        assert job.ts_sp_job.settings.basis == "def2tzvp"
        assert job.ts_sp_job.settings.solvent_model == "smd"
        assert job.ts_sp_job.settings.freq is False
        assert job.ts_opt_job.settings.functional == "b3lyp"
        assert job.ts_opt_job.settings.basis == "def2svp"

    def test_no_sp_skips_sp_phase(self, gaussian_jobrunner_no_scratch):
        job = GaussianReactionJob(
            molecule=_h2(0.82),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
            no_sp=True,
        )
        assert job.phase_by_name("SP").should_skip()
        assert job.sp_jobs is None

    def test_missing_reactant_and_product_omitted(
        self, gaussian_jobrunner_no_scratch
    ):
        job = GaussianReactionJob(
            molecule=_h2(0.82),
            settings=_ts_settings(),
            label="sn2",
            jobrunner=gaussian_jobrunner_no_scratch,
        )
        assert job.reactant_opt_jobs == []
        assert job.product_opt_jobs == []
        assert job.opt_jobs == [job.ts_opt_job]


class TestORCAReactionJob:
    def test_case1_skips_guess_and_builds_ts_opt(
        self, orca_jobrunner_no_scratch
    ):
        job = ORCAReactionJob(
            molecule=_h2(0.82),
            settings=_orca_settings("ts"),
            label="sn2",
            jobrunner=orca_jobrunner_no_scratch,
            opt_settings=_orca_settings("opt"),
            ts_settings=_orca_settings("ts"),
            sp_settings=_orca_settings(
                "sp", functional="wb97x", basis="def2-tzvp", freq=False
            ),
        )
        assert isinstance(job, ORCAJob)
        assert job.TYPE == "orcareaction"
        assert job.phase_by_name("Guess").should_skip()
        assert job.guess_job is None
        assert isinstance(job.ts_opt_job, ORCATSJob)
        assert job.ts_opt_job.label == "sn2_TS_opt"
        assert job.ts_opt_job.settings.jobtype == "ts"
        assert job.ts_opt_job.settings.freq is True
        assert job.opt_jobs == [job.ts_opt_job]

    def test_case2_uses_neb_job_then_ts_job(self, orca_jobrunner_no_scratch):
        job = ORCAReactionJob(
            molecule=_h2(0.74),
            settings=_orca_settings("ts"),
            label="sn2",
            jobrunner=orca_jobrunner_no_scratch,
            products=[_h2(0.90)],
            opt_settings=_orca_settings("opt"),
            ts_settings=_orca_settings("ts"),
            neb_settings=_orca_settings("neb"),
        )
        assert not job.phase_by_name("Guess").should_skip()
        assert isinstance(job.guess_job, ORCANEBJob)
        assert job.guess_job.label == "sn2_neb"
        assert job.guess_job.settings.joboption == "NEB-TS"
        assert job.guess_job.settings.ending_xyzfile.endswith("sn2_neb_P.xyz")
        assert job.guess_job.settings.intermediate_xyzfile is None
        assert isinstance(job.ts_opt_job, ORCATSJob)
        assert isinstance(job.reactant_opt_jobs[0], ORCAOptJob)
        assert isinstance(job.product_opt_jobs[0], ORCAOptJob)

    def test_case2_neb_intermediate_when_ts_guess_given(
        self, orca_jobrunner_no_scratch
    ):
        job = ORCAReactionJob(
            molecule=_h2(0.74),
            settings=_orca_settings("opt"),
            label="sn2",
            jobrunner=orca_jobrunner_no_scratch,
            products=[_h2(0.90)],
            ts_guess=_h2(0.82),
        )
        assert job.guess_job.settings.intermediate_xyzfile.endswith(
            "sn2_neb_TS.xyz"
        )

    def test_case2_reactant_and_product_use_parent_as_neb_guess(
        self, orca_jobrunner_no_scratch
    ):
        job = ORCAReactionJob(
            molecule=_h2(0.82),
            settings=_orca_settings("opt"),
            label="sn2",
            jobrunner=orca_jobrunner_no_scratch,
            reactants=[_h2(0.74)],
            products=[_h2(0.90)],
        )
        assert job.guess_job.settings.intermediate_xyzfile.endswith(
            "sn2_neb_TS.xyz"
        )
        assert job.reactant_opt_jobs[0].label == "sn2_R_opt"

    def test_no_path_search_skips_guess_when_r_p_ts_set(
        self, orca_jobrunner_no_scratch
    ):
        job = ORCAReactionJob(
            molecule=_h2(0.82),
            settings=_orca_settings("ts"),
            label="sn2",
            jobrunner=orca_jobrunner_no_scratch,
            reactants=[_h2(0.74)],
            products=[_h2(0.90)],
            no_path_search=True,
        )
        assert job.phase_by_name("Guess").should_skip()
        assert job.guess_job is None

    def test_sp_uses_solv_theory(self, orca_jobrunner_no_scratch):
        job = ORCAReactionJob(
            molecule=_h2(0.82),
            settings=_orca_settings("ts"),
            label="sn2",
            jobrunner=orca_jobrunner_no_scratch,
            reactants=[_h2(0.74)],
            no_path_search=True,
            opt_settings=_orca_settings("opt"),
            ts_settings=_orca_settings("ts"),
            sp_settings=_orca_settings(
                "sp",
                functional="wb97x",
                basis="def2-tzvp",
                freq=False,
                solvent_model="smd",
                solvent_id="water",
            ),
        )
        job._sp_jobs_for_phase()
        assert job.ts_sp_job.settings.functional == "wb97x"
        assert job.ts_sp_job.settings.basis == "def2-tzvp"
        assert job.ts_sp_job.settings.solvent_id == "water"
        assert isinstance(job.ts_sp_job, ORCASinglePointJob)
        assert job.reactant_sp_jobs[0].label == "sn2_R_sp"

    def test_no_sp_skips_sp_phase(self, orca_jobrunner_no_scratch):
        job = ORCAReactionJob(
            molecule=_h2(0.82),
            settings=_orca_settings("ts"),
            label="sn2",
            jobrunner=orca_jobrunner_no_scratch,
            no_sp=True,
        )
        assert job.phase_by_name("SP").should_skip()
