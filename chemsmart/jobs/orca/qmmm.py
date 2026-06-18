"""ORCA QM/MM jobs and ionic-crystal preparation helpers."""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Optional, Sequence

from chemsmart.jobs.orca.job import ORCAJob
from chemsmart.jobs.orca.writer import CrystalPrepAtomType, CrystalPrepOptions


def _default_prep_stem(input_cif: str, sc_dimension: str) -> str:
    """Build the default supercell file stem, e.g. ``NaCl.cif_15x15x15``."""
    cif_basename = os.path.basename(input_cif)
    return f"{cif_basename}_{sc_dimension}"


@dataclass(frozen=True)
class IonicCrystalPrepFilenames:
    """
    Expected filenames produced by ORCA_crystalprep / ORCA_mm for ionic crystals.

    Individual fields may be overridden when ORCA utility naming differs by
    version or workflow.
    """

    stem: str
    supercell_xyz: str
    supercell_pdb: str
    orcaff_prms: str
    icqmmm_inp: str

    @classmethod
    def infer(
        cls,
        input_cif: str,
        sc_dimension: str,
        *,
        stem: Optional[str] = None,
        supercell_xyz: Optional[str] = None,
        supercell_pdb: Optional[str] = None,
        orcaff_prms: Optional[str] = None,
        icqmmm_inp: Optional[str] = None,
    ) -> IonicCrystalPrepFilenames:
        """
        Infer expected prep filenames from a CIF and supercell dimension.

        Args:
            input_cif: Path or basename of the input CIF (e.g. ``NaCl.cif``).
            sc_dimension: Supercell dimensions (e.g. ``15x15x15``).
            stem: Optional override for the shared basename stem.
            supercell_xyz: Optional override for the supercell XYZ file.
            supercell_pdb: Optional override for the supercell PDB file.
            orcaff_prms: Optional override for the ORCAFF parameter file.
            icqmmm_inp: Optional override for the IC-QMMM input file.

        Returns:
            IonicCrystalPrepFilenames: Resolved filenames.
        """
        resolved_stem = stem or _default_prep_stem(input_cif, sc_dimension)
        return cls(
            stem=resolved_stem,
            supercell_xyz=supercell_xyz or f"{resolved_stem}.xyz",
            supercell_pdb=supercell_pdb or f"{resolved_stem}.pdb",
            orcaff_prms=orcaff_prms or f"{resolved_stem}.ORCAFF.prms",
            icqmmm_inp=icqmmm_inp or f"{resolved_stem}.xyz.ICQMMM.inp",
        )


def _shell_quote(path: str) -> str:
    """Quote a path for shell command generation."""
    return f'"{path}"'


def build_orca_crystalprep_geninput_command(cp_template: str) -> str:
    """Return the ``ORCA_crystalprep ... -geninput`` utility command."""
    return f"ORCA_crystalprep {_shell_quote(cp_template)} -geninput"


def build_orca_crystalprep_run_command(cp_template: str) -> str:
    """Return the ``ORCA_crystalprep`` utility command."""
    return f"ORCA_crystalprep {_shell_quote(cp_template)}"


def build_orca_crystalprep_commands(cp_template: str) -> tuple[str, str]:
    """
    Return CrystalPrep utility commands in execution order.

    Returns:
        tuple[str, str]: ``(-geninput, run)`` command pair.
    """
    return (
        build_orca_crystalprep_geninput_command(cp_template),
        build_orca_crystalprep_run_command(cp_template),
    )


def build_orca_mm_makeff_command(
    supercell_xyz: str,
    atom_types: Sequence[CrystalPrepAtomType],
) -> str:
    """
    Return the ``ORCA_mm -makeff`` command for the supercell XYZ file.

    ``-CEL`` flags are generated from :class:`CrystalPrepAtomType` entries.
    """
    if not atom_types:
        raise ValueError(
            "At least one CrystalPrepAtomType is required for ORCA_mm -makeff."
        )

    cel_flags = " ".join(
        f"-CEL {atom_type.symbol} {atom_type.formal_charge}"
        for atom_type in atom_types
    )
    return f"ORCA_mm -makeff {_shell_quote(supercell_xyz)} {cel_flags}"


def build_ionic_crystal_prep_commands(
    cp_template: str,
    supercell_xyz: str,
    atom_types: Sequence[CrystalPrepAtomType],
) -> tuple[str, str, str]:
    """
    Return all ionic-crystal prep utility commands in execution order.

    Returns:
        tuple[str, str, str]: ``(crystalprep_geninput, crystalprep_run, orca_mm_makeff)``
    """
    geninput_cmd, run_cmd = build_orca_crystalprep_commands(cp_template)
    makeff_cmd = build_orca_mm_makeff_command(supercell_xyz, atom_types)
    return geninput_cmd, run_cmd, makeff_cmd


class ORCAQMMMJob(ORCAJob):
    TYPE = "orcaqmmm"

    def __init__(
        self, molecule, settings, label, structure_filename=None, **kwargs
    ):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
        self.structure_filename = structure_filename

    @staticmethod
    def infer_prep_filenames(
        input_cif: str,
        sc_dimension: str,
        **overrides,
    ) -> IonicCrystalPrepFilenames:
        """Infer expected ionic-crystal prep filenames from CIF and dimensions."""
        return IonicCrystalPrepFilenames.infer(
            input_cif, sc_dimension, **overrides
        )

    @classmethod
    def infer_prep_filenames_from_options(
        cls,
        options: CrystalPrepOptions,
        **overrides,
    ) -> IonicCrystalPrepFilenames:
        """Infer expected prep filenames from :class:`CrystalPrepOptions`."""
        return cls.infer_prep_filenames(
            options.input_cif,
            options.sc_dimension,
            **overrides,
        )

    @staticmethod
    def build_prep_commands(
        cp_template: str,
        supercell_xyz: str,
        atom_types: Sequence[CrystalPrepAtomType],
    ) -> tuple[str, str, str]:
        """Build ORCA_crystalprep and ORCA_mm prep utility commands."""
        return build_ionic_crystal_prep_commands(
            cp_template, supercell_xyz, atom_types
        )

    @classmethod
    def build_prep_commands_from_options(
        cls,
        cp_template: str,
        options: CrystalPrepOptions,
        **filename_overrides,
    ) -> tuple[str, str, str]:
        """Build prep commands using filenames inferred from CrystalPrep options."""
        filenames = cls.infer_prep_filenames_from_options(
            options, **filename_overrides
        )
        return cls.build_prep_commands(
            cp_template,
            filenames.supercell_xyz,
            options.atom_types,
        )
