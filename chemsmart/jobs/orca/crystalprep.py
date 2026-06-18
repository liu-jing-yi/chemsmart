"""ORCA ionic-crystal CrystalPrep prejob execution on the compute node."""

from __future__ import annotations

import logging
import os
import shlex
import subprocess
from collections.abc import Sequence
from typing import Optional

import click

from chemsmart.cli.orca.qmmm import (
    build_crystalprep_options_from_cli,
    write_crystalprep_template,
)
from chemsmart.jobs.orca.qmmm import (
    IonicCrystalPrepFilenames,
    build_ionic_crystal_prep_commands,
)
from chemsmart.jobs.orca.writer import CrystalPrepOptions

logger = logging.getLogger(__name__)

_IONIC_CRYSTAL_QMMM_JOBTYPE = "IONIC-CRYSTAL-QMMM"


class CrystalPrepPrejobError(RuntimeError):
    """Raised when ionic-crystal CrystalPrep prejob preparation fails."""


def potential_orca_crystalprep_prejob(cli_args: Sequence[str]) -> None:
    """
    Run ORCA CrystalPrep / ORCA_mm preparation when CLI args request it.

    No-op unless ``cli_args`` describe an ORCA QMMM job with
    ``--jobtype IONIC-CRYSTAL-QMMM`` and ``--geninput``.
    """
    options = _crystalprep_options_from_cli_args(cli_args)
    if options is None:
        return
    run_orca_crystalprep_prejob(options)


def _crystalprep_options_from_cli_args(
    cli_args: Sequence[str],
) -> Optional[CrystalPrepOptions]:
    args = list(cli_args)
    if not _should_run_crystalprep_prejob(args):
        return None

    structure_filename = _get_option_value(args, "-f", "--filename")
    try:
        return build_crystalprep_options_from_cli(
            cp_input_cif=_get_option_value(args, "--cp-input-cif"),
            structure_filename=structure_filename,
            cp_scdimension=_get_option_value(args, "--cp-scdimension"),
            cp_atomtype=_parse_cp_atomtypes(args),
            cp_docif=_get_optional_bool(args, "--cp-docif"),
            cp_dosupercell=_get_optional_bool(args, "--cp-dosupercell"),
            cp_doembedding=_get_optional_bool(args, "--cp-doembedding"),
            cp_dolayers=_get_optional_bool(args, "--cp-dolayers"),
            cp_doicqmmminput=_get_optional_bool(args, "--cp-doicqmmminput"),
            cp_qccharge=_get_optional_int(args, "--cp-qccharge"),
            cp_qcmult=_get_optional_int(args, "--cp-qcmult"),
            cp_neutralize=_get_optional_bool(args, "--cp-neutralize"),
            cp_template_out=_get_option_value(args, "--cp-template-out"),
        )
    except click.UsageError as exc:
        raise CrystalPrepPrejobError(str(exc)) from exc


def _should_run_crystalprep_prejob(args: Sequence[str]) -> bool:
    if "orca" not in args or "qmmm" not in args or "--geninput" not in args:
        return False
    jobtype = _get_option_value(args, "-j", "--jobtype")
    return (
        jobtype is not None and jobtype.upper() == _IONIC_CRYSTAL_QMMM_JOBTYPE
    )


def run_orca_crystalprep_prejob(options: CrystalPrepOptions) -> None:
    """Write the CrystalPrep template and run ORCA utility commands."""
    cp_template = write_crystalprep_template(options)
    _require_file(cp_template, "CrystalPrep template")

    filenames = IonicCrystalPrepFilenames.infer(
        options.input_cif,
        options.sc_dimension,
    )
    geninput_cmd, run_cmd, makeff_cmd = build_ionic_crystal_prep_commands(
        cp_template,
        filenames.supercell_xyz,
        options.atom_types,
    )

    _run_utility_command(geninput_cmd, "ORCA_crystalprep -geninput")
    _require_file(cp_template, "CrystalPrep template after -geninput")

    _run_utility_command(run_cmd, "ORCA_crystalprep")
    _require_file(filenames.supercell_xyz, "Supercell XYZ file")
    _require_file(filenames.supercell_pdb, "Supercell PDB file")

    _run_utility_command(makeff_cmd, "ORCA_mm -makeff")
    _require_file(filenames.orcaff_prms, "ORCAFF parameter file")

    if options.do_icqmmm_input:
        _require_file(filenames.icqmmm_inp, "IC-QMMM input file")


def _get_option_value(
    args: Sequence[str], *option_names: str
) -> Optional[str]:
    for name in option_names:
        if name not in args:
            continue
        index = args.index(name)
        if index + 1 >= len(args):
            raise CrystalPrepPrejobError(
                f"Missing value for option {name} in CLI_ARGS"
            )
        return args[index + 1]
    return None


def _get_optional_int(args: Sequence[str], option_name: str) -> Optional[int]:
    value = _get_option_value(args, option_name)
    if value is None:
        return None
    return int(value)


def _get_optional_bool(
    args: Sequence[str], option_name: str
) -> Optional[bool]:
    if option_name not in args:
        return None
    index = args.index(option_name)
    if index + 1 < len(args) and args[index + 1].lower() in {"true", "false"}:
        return args[index + 1].lower() == "true"
    return True


def _parse_cp_atomtypes(
    args: Sequence[str],
) -> list[tuple[str, int, float, float]]:
    atomtypes: list[tuple[str, int, float, float]] = []
    index = 0
    while index < len(args):
        if args[index] != "--cp-atomtype":
            index += 1
            continue
        if index + 4 >= len(args):
            raise CrystalPrepPrejobError(
                "Incomplete --cp-atomtype specification in CLI_ARGS"
            )
        atomtypes.append(
            (
                args[index + 1],
                int(args[index + 2]),
                float(args[index + 3]),
                float(args[index + 4]),
            )
        )
        index += 5
    return atomtypes


def _require_file(path: str, description: str) -> None:
    if not os.path.isfile(path):
        raise CrystalPrepPrejobError(f"{description} not found: {path}")


def _run_utility_command(command: str, description: str) -> None:
    argv = shlex.split(command)
    logger.info("Running %s: %s", description, command)
    try:
        subprocess.run(argv, check=True)
    except subprocess.CalledProcessError as exc:
        raise CrystalPrepPrejobError(
            f"{description} failed with exit code {exc.returncode}: {command}"
        ) from exc
