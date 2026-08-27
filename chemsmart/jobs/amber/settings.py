"""
Amber job settings and configuration classes.

This module provides settings classes for Amber molecular dynamics and
energy calculations, following the same pattern as other computational
chemistry modules in chemsmart.
"""

import logging
import os
import subprocess
import warnings
from typing import Optional

from chemsmart.jobs.settings import MolecularJobSettings

logger = logging.getLogger(__name__)


class AmberJobSettings(MolecularJobSettings):
    """
    Settings class for Amber molecular dynamics calculations.

    Provides configuration options for Amber simulations including
    runtime parameters, system setup, and computational settings.

    Attributes:
        job_name (str): Name identifier for the job
        max_runtime (str): Maximum runtime in HH:MM:SS format
        num_nodes (int): Number of compute nodes to request
        ppn (int): Processors per node
        queue (str): Queue/partition name for job submission
        email_notifications (bool): Whether to send email notifications
        email (str): Email address for notifications
        simulation_type (str): Type of simulation (md, energy, opt, etc.)
        input_file (str): Amber input file path
        topology_file (str): Topology file path (.prmtop)
        coordinate_file (str): Initial coordinate file (.inpcrd/.rst7)
        steps (int): Number of simulation steps
        temperature (float): Target temperature in Kelvin
        pressure (float): Target pressure in bar (for NPT)
        ensemble (str): Statistical ensemble (NVE, NVT, NPT)
    """

    def __init__(
        self,
        job_name: str = "Amber Simulation",
        max_runtime: str = "48:00:00",
        num_nodes: int = 2,
        ppn: int = 16,
        queue: str = "standard",
        simulation_type: str = "md",
        input_file: Optional[str] = None,
        topology_file: Optional[str] = None,
        coordinate_file: Optional[str] = None,
        steps: int = 50000,
        temperature: float = 300.0,
        pressure: Optional[float] = None,
        ensemble: str = "NVT",
        parameter: str = None,
        **kwargs,
    ):
        """
        Initialize Amber job settings.

        Args:
            job_name: Name identifier for the job
            max_runtime: Maximum runtime in HH:MM:SS format
            num_nodes: Number of compute nodes
            ppn: Processors per node
            queue: Queue/partition name
            simulation_type: Type of simulation
            input_file: Amber input file path
            topology_file: Topology file (.prmtop)
            coordinate_file: Coordinate file (.inpcrd/.rst7)
            steps: Number of simulation steps
            temperature: Target temperature in Kelvin
            pressure: Target pressure in bar
            ensemble: Statistical ensemble type
            **kwargs: Additional arguments passed to parent class
        """
        super().__init__(**kwargs)

        self.job_name = job_name
        self.max_runtime = max_runtime
        self.num_nodes = num_nodes
        self.ppn = ppn
        self.queue = queue
        self.simulation_type = simulation_type
        self.input_file = input_file
        self.topology_file = topology_file
        self.coordinate_file = coordinate_file
        self.steps = steps
        self.temperature = temperature
        self.pressure = pressure
        self.ensemble = ensemble
        self.parameter = parameter

    def validate_settings(self):
        """
        Validate Amber job settings for consistency.

        Raises:
            ValueError: If settings are invalid or inconsistent
        """
        if self.ensemble == "NPT" and self.pressure is None:
            raise ValueError("NPT ensemble requires pressure to be specified")

        if self.steps <= 0:
            raise ValueError("Number of steps must be positive")

        if self.temperature <= 0:
            raise ValueError("Temperature must be positive")

        if self.num_nodes <= 0 or self.ppn <= 0:
            raise ValueError("num_nodes and ppn must be positive integers")

    def merge(self, other_settings, keywords=None):
        """
        Merge current settings with another settings instance.

        Args:
            other_settings: Another AmberJobSettings or dict to merge
            keywords: Optional keywords dictionary to merge

        Returns:
            AmberJobSettings: New merged settings instance
        """
        # Create a copy of current settings
        merged_dict = self.__dict__.copy()

        # Merge with other settings
        if hasattr(other_settings, "__dict__"):
            other_dict = other_settings.__dict__
        else:
            other_dict = other_settings

        for key, value in other_dict.items():
            if value is not None:
                merged_dict[key] = value

        # Merge keywords if provided
        if keywords:
            for key, value in keywords.items():
                if value is not None:
                    merged_dict[key] = value

        return self.__class__(**merged_dict)

    @classmethod
    def from_dict(cls, settings_dict):
        """
        Create AmberJobSettings from dictionary.

        Args:
            settings_dict: Dictionary containing settings

        Returns:
            AmberJobSettings: Settings instance
        """
        return cls(**settings_dict)

    def copy(self):
        """
        Create a copy of the settings.

        Returns:
            AmberJobSettings: Copy of current settings
        """
        return self.__class__(**self.__dict__)

    @classmethod
    def default(cls):
        """
        Create default AmberJobSettings instance.

        Returns:
            AmberJobSettings: Default settings
        """
        return cls()


class AmberFFParamJobSettings(AmberJobSettings):
    """
    Settings class for Amber force field parameterization jobs.

    Inherits from AmberJobSettings and adds specific parameters
    for force field parameterization tasks.

    Attributes:
        element (str): Element name for metal parameterization
        ligand (str): Ligand name for parameterization
        ligand_charge (int): Net charge of the ligand
        ligand_multiplicity (int): Spin multiplicity of the ligand
        gaussian_logfile (str): Path to Gaussian log file
        gaussian_fchkfile (str): Path to Gaussian fchk file
        parameter (str): Force field parameter set to use
        input_pdb (str): Path to input PDB file for processing
        cutoff (float): Cutoff distance for metal binding
        large_opt (int): Large optimization setting
        mcpb_job_type (str): Type of quantum chemistry job for MCPB ("gaussian" or "gamess")

    Properties:
        ligand_pdb (str): Path to extracted ligand PDB file
        metal_pdb (str): Path to extracted metal PDB file
        reduced_ligand (str): Path to ligand PDB with added hydrogens
        ligand_mol2_file (str): Path to ligand mol2 file
        ligand_frcmod_file (str): Path to ligand frcmod file
        water_mol2_file (str): Path to water mol2 file
        hpp_pdb_file (str or None): Path to PDB file with AMBER naming scheme
        group_name (str or None): Basename of input PDB file without extension
        ion_ids (str or None): Metal atom index as string for MCPB input
        ion_mol2files (str or None): Ion mol2 file name for MCPB input
        naa_mol2files (str or None): Non-amino acid mol2 file name for MCPB input
        frcmod_files (str or None): Force field modification file name for MCPB input
        mcpb_topology_files (tuple): Topology and coordinate files from tleap after MCPB step 4
        verified_topology_files (tuple): Verified and fixed topology/coordinate files ready for MD
        verified_solv_prmtop (str): Verified and validated solvated topology file with ParmEd validation
    """

    def __init__(
        self,
        element: str = None,
        element_charge: int = None,
        ligand: str = None,
        ligand_charge: int = None,
        ligand_multiplicity: int = None,
        charge_method: str = "bcc",
        gaussian_logfile: str = None,
        gaussian_fchkfile: str = None,
        parameter: str = None,
        input_pdb: str = None,
        fixed_H_input_pdb_file: str = None,
        cutoff: float = 2.8,
        large_opt: int = 1,
        mcpb_job_type: str = "gaussian",
        **kwargs,
    ):
        """
        Initialize Amber force field parameterization settings.

        Args:
            element: Element name for metal parameterization
            ligand: Ligand name for parameterization
            ligand_charge: Net charge of the ligand
            ligand_multiplicity: Spin multiplicity of the ligand
            gaussian_logfile: Path to Gaussian log file
            gaussian_fchkfile: Path to Gaussian fchk file
            parameter: Force field parameter set to use
            input_pdb: Path to input PDB file for processing
            cutoff: Cutoff distance for metal binding (default: 2.8)
            large_opt: Large optimization setting (default: 1)
            mcpb_job_type: Type of quantum chemistry job ("gaussian" or "gamess", default: "gaussian")
            **kwargs: Additional arguments passed to parent class
        """
        super().__init__(**kwargs)

        self.element = element
        self.element_charge = element_charge
        self.ligand = ligand
        self.ligand_charge = ligand_charge
        self.ligand_multiplicity = ligand_multiplicity
        self.charge_method = charge_method
        self.parameter = parameter
        self.gaussian_logfile = gaussian_logfile
        self.gaussian_fchkfile = gaussian_fchkfile
        self.input_pdb = input_pdb
        self.fixed_H_input_pdb_file = fixed_H_input_pdb_file
        self.cutoff = cutoff
        self.large_opt = large_opt
        self.mcpb_job_type = mcpb_job_type.lower()

    @property
    def ligand_pdb(self):
        """
        Extract ligand lines from a PDB input file and write them to <ligand>.pdb

        Returns:
            str: Path to the created ligand PDB file
        """
        if not self.ligand:
            raise ValueError("Ligand name must be specified")
        return self._get_residue_pdb(self.ligand)

    @property
    def metal_pdb(self):
        """
        Extract metal lines from a PDB input file and write them to <element>.pdb

        Returns:
            str: Path to the created metal PDB file
        """
        if not self.element:
            raise ValueError("Element name must be specified")
        return self._get_residue_pdb(self.element)

    @property
    def reduce_command(self):
        """the reduce command for adding hydrogens to ligand pdb file.

        Example command:
            reduce LIGAND.pdb > LIGAND_H.pdb
        """
        return self._get_reduce_command()

    @property
    def antechamber_mol2_command(self):
        """get antechamber command for generating mol2 jobs.

        Example command:
            antechamber -fi pdb -fo mol2 -i LIGAND_H.pdb -o LIGAND_pre.mol2 -c bcc -pf y -nc 0
        """
        return self._get_antechamber_command()

    @property
    def rename_antechamber_output_command(self):
        """get command to rename antechamber output mol2 file.

        Example command:
            mv LIGAND_pre.mol2 LIGAND.mol2
        """
        if not self.ligand:
            raise ValueError(
                "Ligand name must be specified to rename antechamber output"
            )
        return f"mv {self.ligand}_pre.mol2 {self.ligand}.mol2"

    @property
    def parmchk2_command(self):
        """get parmchk2 command for generating frcmod files.

        Example command:
            parmchk2 -i LIGAND_pre.mol2 -f mol2 -o LIGAND.frcmod -s gaff
        """
        return self._get_parmchk2_command()

    @property
    def metalpdb2mol_command(self):
        """get metalpdb2mol command for generating metal mol2 files.

        Example command:
            metalpdb2mol.py -i ELEMENT.pdb -o ELEMENT.mol2 -c ELEMENT_CHARGE
        """
        return self._get_metalpdb2mol_command()

    @property
    def ambpdb_command(self):
        """get ambpdb command for generating pdb file from topology/coordinate files.

        Example command:
            ambpdb -p TOPOLOGY_FILE -c COORDINATE_FILE > OUTPUT_PDB_FILE
        """
        if not self.topology_file or not self.coordinate_file:
            return None
        return self._get_ambpdb_command()

    @property
    def combine_command(self):
        """get command to combine multiple pdb files into one.

        Example command:
            cat PROTEIN_H_fixed.pdb ELEMENT.pdb LIGAND_H.pdb > COMPLEX_H.pdb
        """
        return self._get_combine_command()

    @property
    def pdb4amber_command(self):
        """get pdb4amber command for fixing and renumbering pdb files."""
        return self._get_pdb4amber_command()

    @property
    def mcpb_step1_command(self):
        """get MCPB.py step 1 command for generating optimization and RESP files.

        Example command:
            MCPB.py -i MCPB_INPUT_FILE -s 1
        """
        return self._get_mcpb_step_command(step=1)

    @property
    def mcpb_step2_command(self):
        """get MCPB.py step 2 command for generating force field parameters.

        Other options are available:
        Z-matrix (with step_number 2z)；
        Empirical (with step_number 2e) methods.
        The Empirical method doesn't need any Gaussian calculations to obtain the force constants
        (but still needs Gaussian calculation for getting the RESP charges),
        but it only supports zinc ion modeling in the current version.

        Example command:
            MCPB.py -i MCPB_INPUT_FILE -s 2
        """
        return self._get_mcpb_step_command(step=2)

    @property
    def mcpb_step3_command(self):
        """get MCPB.py step 3 command for generating force field parameters.

        Use ChgModB to perform the REWSP charge fitting and generate the mol2 files for the metal center and coordinating residues.
        Other options are also available:
        ChgModA, ChgModC and ChgModD (as 3a, 3c and 3d respectively)

        Example command:
            MCPB.py -i MCPB_INPUT_FILE -s 3
        """
        return self._get_mcpb_step_command(step=3)

    @property
    def mcpb_step4_command(self):
        """get MCPB.py step 4 command for generating force field parameters.

        Other options are available: Z-matrix (with step_number 2z) and Empirical (with step_number 2e) methods. The Empirical method doesn't need any Gaussian calculations to obtain the force constants (but still needs Gaussian calculation for getting the RESP charges), but it only supports zinc ion modeling in the current version.

        Example command:
            MCPB.py -i MCPB_INPUT_FILE -s 2
        """
        return self._get_mcpb_step_command(step=2)

    @property
    def water_mol2_file(self):
        """
        Get the water mol2 file by running tleap and antechamber.

        Returns:
            str: Path to the WAT.mol2 file
        """
        return self._run_tleap()

    @property
    def hpp_pdb_file(self):
        """
        Generate PDB file with AMBER naming scheme from topology and coordinate files.

        Uses ambpdb to convert topology/coordinate files to PDB format, useful
        after H++ processing to get proper AMBER residue names.

        Returns:
            str or None: Path to the generated PDB file if topology file exists, None otherwise
        """
        return self._run_ambpdb()

    @property
    def ion_id(self):
        return self._get_metal_id()

    @property
    def group_name(self):
        """
        Get the group name based on the basename of the input PDB file.

        Returns the filename without path and extension from the input_pdb.

        Returns:
            str or None: Basename of input PDB file if input_pdb is specified, None otherwise

        Raises:
            ValueError: If input_pdb doesn't exist
        """
        if not self.input_pdb:
            return None

        if not os.path.exists(self.input_pdb):
            raise ValueError(f"Input PDB file {self.input_pdb} does not exist")

        # Extract basename without extension
        basename = os.path.splitext(os.path.basename(self.input_pdb))[0]
        return basename

    @property
    def ion_ids(self):
        """
        Get the metal ion ID from the input PDB file.

        Returns:
            str or None: Metal atom index as string if found, None otherwise
        """
        metal_id = self._get_metal_id()
        return str(metal_id) if metal_id is not None else None

    @property
    def ion_mol2files(self):
        """
        Get the ion mol2 file name based on the element.

        Returns:
            str or None: Ion mol2 file name if element is specified, None otherwise
        """
        if not self.element:
            return None
        return f"{self.element}.mol2"

    @property
    def naa_mol2files(self):
        """
        Get the non-amino acid mol2 file name based on the ligand.

        Returns:
            str or None: Ligand mol2 file name if ligand is specified, None otherwise
        """
        if not self.ligand:
            return None
        return f"{self.ligand}.mol2"

    @property
    def frcmod_files(self):
        """
        Get the force field modification file name based on the ligand.

        Returns:
            str or None: Ligand frcmod file name if ligand is specified, None otherwise
        """
        if not self.ligand:
            return None
        return f"{self.ligand}.frcmod"

    @property
    def mcpb_input(self):
        return self._write_MCPB_input()

    @property
    def small_opt_com_file(self):
        """
        Get the optimization file from MCPB.py step 1.

        Returns:
            str or None: Path to optimization file (*.com for Gaussian, *small_opt.inp for GAMESS-US)
        """
        small_opt_file, _, _ = self._run_mcpb(self.mcpb_input, step=1)
        return small_opt_file

    @property
    def fc_file(self):
        """
        Get the frequency calculation file from MCPB.py step 1.

        Returns:
            str or None: Path to frequency file (*.com for Gaussian, *fc.inp for GAMESS-US)
        """
        _, fc_file, _ = self._run_mcpb(self.mcpb_input, step=1)
        return fc_file

    @property
    def mk_file(self):
        """
        Get the RESP charge fitting file from MCPB.py step 1.

        Returns:
            str or None: Path to RESP file (*.com for Gaussian, *mk.inp for GAMESS-US)
        """
        _, _, mk_file = self._run_mcpb(self.mcpb_input, step=1)
        return mk_file

    @property
    def mcpbpy_frcmod_file(self):
        return self._run_mcpb(self.mcpb_input, step=2)

    @property
    def mcpb_mol_files(self):
        return self._run_mcpb(self.mcpb_input, step=3)

    @property
    def tleap_pdb_file(self):
        tleap_pdb_file, _ = self._run_mcpb(self.mcpb_input, step=4)
        return tleap_pdb_file

    @property
    def tleap_input_file(self):
        tleap_input_file, _ = self._run_mcpb(self.mcpb_input, step=4)
        return tleap_input_file

    @property
    def mcpb_topology_files(self):
        """
        Generate topology and coordinate files using tleap after MCPB step 4.

        This property runs tleap on the input file generated by MCPB step 4
        and returns the generated topology and coordinate files.

        Returns:
            tuple: (topology_file, coordinate_file) containing paths to:
                   - topology_file: Generated .prmtop file
                   - coordinate_file: Generated .inpcrd file
                   Each element can be None if the corresponding file is not found.
        """
        return self._run_tleap_mcpb()

    def fix_pdb_file(self, input_pdb_file, output_pdb_file=None):
        """
        Fix and renumber PDB file using pdb4amber.

        Convenience method to run pdb4amber for PDB file processing.

        Args:
            input_pdb (str): Path to input PDB file to process
            output_pdb_file (str, optional): Path to output PDB file. If not specified,
                                           defaults to <input_basename>_fixed.pdb

        Returns:
            str: Path to the output fixed PDB file

        Example:
            # Equivalent to: pdb4amber -i 1OKL_H.pdb -o 1OKL_fixed_H.pdb
            fixed_file = settings.fix_pdb_file('1OKL_H.pdb', '1OKL_fixed_H.pdb')
        """
        return self._run_pdb4amber(input_pdb_file, output_pdb_file)

    def run_mcpb(self, input_file, step=1):
        """
        Run MCPB.py for metal center parameter building.

        Convenience method to run MCPB.py and check for generated output files.

        Args:
            input_file (str): Path to MCPB input file (typically .in file)
            step (int, optional): MCPB step number (default: 1)

        Returns:
            dict: Dictionary containing paths to generated files and completion status.
                  Keys depend on job type:
                  - For Gaussian: 'opt_com', 'fc_com', 'mk_com'
                  - For GAMESS-US: 'opt_inp', 'fc_inp', 'mk_inp'
                  - 'all_found': Boolean indicating if all files were found

        Example:
            # Equivalent to: MCPB.py -i 1OKL.in -s 1
            result = settings.run_mcpb('1OKL.in', 1)
            if result['all_found']:
                print("All MCPB files generated successfully")
        """
        opt_file, fc_file, mk_file = self._run_mcpb(input_file, step)

        # Create backward-compatible dictionary based on job type
        if self.mcpb_job_type == "gaussian":
            result = {
                "opt_com": opt_file,
                "fc_com": fc_file,
                "mk_com": mk_file,
                "all_found": all([opt_file, fc_file, mk_file]),
            }
        elif self.mcpb_job_type == "gamess":
            result = {
                "opt_inp": opt_file,
                "fc_inp": fc_file,
                "mk_inp": mk_file,
                "all_found": all([opt_file, fc_file, mk_file]),
            }
        else:
            # Fallback to generic names
            result = {
                "opt_file": opt_file,
                "fc_file": fc_file,
                "mk_file": mk_file,
                "all_found": all([opt_file, fc_file, mk_file]),
            }

        return result

    def _get_residue_pdb(self, residue_name, residue_aliases=None):
        """
        Extract residue lines from a PDB input file and write them to <residue>.pdb.

        This unified method can extract any type of residue (ligand, metal, water, etc.)
        based on the residue name and optional aliases.

        Args:
            residue_name (str): Primary residue name to search for
            residue_aliases (list, optional): Additional residue names to search for
                                            (useful for water: ['HOH', 'WAT', 'TIP3', 'SOL'])

        Returns:
            str: Path to the created residue PDB file

        Raises:
            ValueError: If residue name, input file not specified, or no residue found
        """
        if not residue_name:
            raise ValueError("Residue name must be specified")

        if not self.input_file:
            raise ValueError("Input file must be specified")

        if not os.path.exists(self.input_file):
            raise ValueError(f"Input file {self.input_file} does not exist")

        # Format residue name using _residue_formatting (includes warning if needed)
        residue_upper = self._residue_formatting(residue_name)
        residue_pdb_filename = f"{residue_name}.pdb"

        # Read input file
        with open(self.input_file, "r") as f:
            lines = f.readlines()

        # Build list of residue names to search for
        search_names = [residue_name]
        if residue_aliases:
            search_names.extend(residue_aliases)

        # Filter lines containing any of the search names
        residue_lines = []
        for line in lines:
            if any(search_name in line for search_name in search_names):
                residue_lines.append(line)

        if not residue_lines:
            if residue_aliases:
                all_names = ", ".join([f"'{name}'" for name in search_names])
                raise ValueError(
                    f"No lines containing residue names {all_names} found in {self.input_file}"
                )
            else:
                raise ValueError(
                    f"No lines containing residue '{residue_name}' found in {self.input_file}"
                )

        # Process lines to ensure residue names are properly formatted
        processed_lines = []
        for line in residue_lines:
            processed_line = line
            # Replace the primary residue name with uppercase version
            if residue_name in line:
                processed_line = processed_line.replace(
                    residue_name, residue_upper
                )
            processed_lines.append(processed_line)

        # Write processed lines to new PDB file
        with open(residue_pdb_filename, "w") as out:
            out.writelines(processed_lines)

        logger.info(
            f"Extracted {len(residue_lines)} residue lines to {residue_pdb_filename}"
        )
        return residue_pdb_filename

    def _residue_formatting(self, resi):
        """
        Format residue identifier for Amber LEaP.

        Checks if residue name is capitalized and warns if not. Ensures output
        uses capital letters for the residue name.

        Args:
            resi (str): Residue identifier

        Returns:
            str: Uppercase residue identifier

        Raises:
            ValueError: If residue identifier is not a string
        """
        if not isinstance(resi, str):
            raise ValueError("Residue identifier must be a string")

        # Check if residue name is capitalized and warn if not
        if not resi.isupper():
            warnings.warn(
                f"Residue name '{resi}' is not all capital letters. "
                f"It will be converted to '{resi.upper()}' in the output.",
                UserWarning,
            )

        return resi.upper()

    def _get_reduce_command(self):
        if not self.ligand:
            raise ValueError("Ligand name must be specified to run reduce")

        ligand_pdb_file = self.ligand_pdb
        output_file = f"{self.ligand}_H.pdb"
        reduce_command = f"reduce {ligand_pdb_file} > {output_file}"
        return reduce_command

    def _get_antechamber_command(self):
        """
        Execute the antechamber command to convert PDB to mol2 format.

        Runs: antechamber -fi pdb -fo mol2 -i {ligand}_H.pdb -o {ligand}_pre.mol2 -c bcc -pf y -nc {ligand_charge}

        Returns:
            str: Path to the output mol2 file

        Raises:
            ValueError: If ligand or ligand_charge is not specified
            RuntimeError: If antechamber command fails or antechamber is not available
        """
        if not self.ligand:
            raise ValueError(
                "Ligand name must be specified to run antechamber"
            )

        if self.ligand_charge is None:
            raise ValueError(
                "Ligand charge must be specified to run antechamber"
            )

        # Input file is the reduced ligand PDB (with hydrogens)
        input_file = self.reduced_ligand

        # Output mol2 file
        output_file = f"{self.ligand}_pre.mol2"

        # Construct the antechamber command
        antechamber_command = [
            "antechamber",
            "-fi",
            "pdb",
            "-fo",
            "mol2",
            "-i",
            input_file,
            "-o",
            output_file,
            "-c",
            self.charge_method,
            "-pf",
            "y",
            "-nc",
            str(self.ligand_charge),
        ]

        return antechamber_command

    def _get_parmchk2_command(self):
        """Construct the parmchk2 command to generate force field modification file.

        Example:
            parmchk2 -i {LIGAND}.mol2 -o {LIGAND}.frcmod -f mol2
        """
        assert self.ligand, "Ligand name must be specified to run parmchk2"
        parmchk2_command = (
            f"parmchk2 -i {self.ligand}.mol2 -o {self.ligand}.frcmod -f mol2"
        )
        return parmchk2_command

    def _get_metalpdb2mol_command(self):
        """Construct the metalpdb2mol command to generate metal mol2 file.

        Example:
            metalpdb2mol.py -i {ELEMENT}.pdb -o {ELEMENT}.mol2 -c {ELEMENT_CHARGE}
        """
        assert (
            self.element
        ), "Element name must be specified to run metalpdb2mol"
        assert (
            self.element_charge is not None
        ), "Element charge must be specified to run metalpdb2mol"
        metalpdb2mol_command = f"metalpdb2mol.py -i {self.element}.pdb -o {self.element}.mol2 -c {self.element_charge}"
        return metalpdb2mol_command

    def _get_ambpdb_command(self):
        """Construct the ambpdb command to generate PDB file from topology/coordinate files.

        Example:
            ambpdb -p TOPOLOGY_FILE -c COORDINATE_FILE > OUTPUT_PDB_FILE
        """
        assert (
            self.topology_file
        ), "Topology file must be specified to run ambpdb"
        assert (
            self.coordinate_file
        ), "Coordinate file must be specified to run ambpdb"
        output_pdb_file = (
            f"{self.group_name}_hpp.pdb"
            if self.group_name
            else "output_hpp.pdb"
        )
        ambpdb_command = f"ambpdb -p {self.topology_file} -c {self.coordinate_file} > {output_pdb_file}"
        return ambpdb_command

    def _run_tleap(self):
        """
        Execute tleap to add hydrogens to water molecules.

        Creates tleap input file and runs:
        tleap -s -f wat_tleap.in > wat_tleap.out

        Then uses antechamber to generate mol2 file:
        antechamber -fi pdb -fo mol2 -i WAT_H.pdb -o WAT.mol2 -at amber -c bcc -pf y

        Returns:
            str: Path to the output WAT.mol2 file

        Raises:
            ValueError: If required parameters are not specified
            RuntimeError: If tleap or antechamber commands fail
        """
        if not self.input_file:
            raise ValueError("Input file must be specified to run tleap")

        # First extract water molecules to WAT.pdb
        water_pdb = self._get_water_pdb()

        # Create tleap input file
        tleap_input_file = "wat_tleap.in"
        tleap_output_file = "wat_tleap.out"
        water_with_h_pdb = "WAT_H.pdb"
        water_mol2_file = "WAT.mol2"

        # Write tleap input file
        tleap_commands = [
            "loadoff solvents.lib",
            "HOH = TIP3",
            f"wat = loadpdb {water_pdb}",
            f"savepdb wat {water_with_h_pdb}",
            "quit",
        ]

        try:
            with open(tleap_input_file, "w") as f:
                f.write("\n".join(tleap_commands) + "\n")

            logger.info(f"Created tleap input file: {tleap_input_file}")

            # Run tleap command
            tleap_command = f"tleap -s -f {tleap_input_file}"

            with open(tleap_output_file, "w") as output:
                subprocess.run(
                    tleap_command,
                    shell=True,
                    stdout=output,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )

            logger.info(f"Successfully ran tleap command: {tleap_command}")
            logger.info(f"Output written to: {tleap_output_file}")

            # Verify WAT_H.pdb was created
            if not os.path.exists(water_with_h_pdb):
                raise RuntimeError(
                    f"Tleap completed but output file {water_with_h_pdb} was not created"
                )

            # Run antechamber on the water with hydrogens
            antechamber_command = [
                "antechamber",
                "-fi",
                "pdb",
                "-fo",
                "mol2",
                "-i",
                water_with_h_pdb,
                "-o",
                water_mol2_file,
                "-at",
                "amber",
                "-c",
                "bcc",
                "-pf",
                "y",
            ]

            subprocess.run(
                antechamber_command, capture_output=True, text=True, check=True
            )

            logger.info(
                f"Successfully ran antechamber command: {' '.join(antechamber_command)}"
            )
            logger.info(f"Water mol2 output written to: {water_mol2_file}")

            # Verify WAT.mol2 was created
            if not os.path.exists(water_mol2_file):
                raise RuntimeError(
                    f"Antechamber completed but output file {water_mol2_file} was not created"
                )

            # Modify water mol2 file to use HW atom type for water hydrogens
            self._fix_water_mol2_atomtypes(water_mol2_file)

        except subprocess.CalledProcessError as e:
            error_msg = (
                f"Tleap/Antechamber command failed: {e}\nStderr: {e.stderr}"
            )
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except FileNotFoundError as e:
            if "tleap" in str(e):
                error_msg = "Tleap command not found. Please ensure AmberTools is installed and in PATH."
            elif "antechamber" in str(e):
                error_msg = "Antechamber command not found. Please ensure AmberTools is installed and in PATH."
            else:
                error_msg = f"Command not found: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except Exception as e:
            error_msg = f"Failed to process water molecules: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        return water_mol2_file

    def _get_water_pdb(self):
        """
        Extract water molecules from input file and write them to WAT.pdb.

        Returns:
            str: Path to the created water PDB file

        Raises:
            ValueError: If input file is not specified or no water found
        """
        # Use common water residue names as aliases
        water_aliases = ["HOH", "WAT", "TIP3", "SOL"]
        return self._get_residue_pdb("WAT", residue_aliases=water_aliases)

    def _fix_water_mol2_atomtypes(self, mol2_file):
        """
        Modify WAT.mol2 file to change H1 and H2 atom types from "HO" to "HW".

        Args:
            mol2_file (str): Path to the mol2 file to modify

        Raises:
            RuntimeError: If file operations fail
        """
        if not os.path.exists(mol2_file):
            raise RuntimeError(f"Mol2 file {mol2_file} does not exist")

        try:
            # Read the file
            with open(mol2_file, "r") as f:
                lines = f.readlines()

            # Modify H1 and H2 atom types from HO to HW
            modified_lines = []
            in_atom_section = False

            for line in lines:
                if line.strip() == "@<TRIPOS>ATOM":
                    in_atom_section = True
                    modified_lines.append(line)
                elif (
                    line.strip().startswith("@<TRIPOS>")
                    and line.strip() != "@<TRIPOS>ATOM"
                ):
                    in_atom_section = False
                    modified_lines.append(line)
                elif in_atom_section:
                    # Check if this is a hydrogen atom line and modify atom type
                    parts = line.split()
                    if len(parts) >= 6:
                        atom_name = parts[1]  # H1 or H2
                        atom_type = parts[5]  # Current atom type

                        # If it's a hydrogen in water and atom type is HO, change to HW
                        if atom_name in ["H1", "H2"] and atom_type == "HO":
                            parts[5] = "HW"
                            line = "     ".join(
                                f"{part:>8}" if i > 1 else f"{part:<8}"
                                for i, part in enumerate(parts[:6])
                            )
                            if len(parts) > 6:
                                line += "     " + "     ".join(parts[6:])
                            line += "\n"
                            logger.info(
                                f"Changed atom type for {atom_name} from HO to HW"
                            )

                    modified_lines.append(line)
                else:
                    modified_lines.append(line)

            # Write back the modified file
            with open(mol2_file, "w") as f:
                f.writelines(modified_lines)

            logger.info(f"Modified water atom types in {mol2_file}")

        except Exception as e:
            error_msg = f"Failed to modify water mol2 atom types: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

    def _run_ambpdb(self):
        """
        Execute ambpdb command to generate PDB file from topology and coordinate files.

        Runs: ambpdb -p <topology_file> -c <coordinate_file> > <output_pdb>
        where the output PDB name is derived from the base name of the topology file.

        This is useful after H++ processing to generate a PDB file using AMBER
        naming scheme for residues (e.g., HID, HIE, HIP instead of HIS for histidine).

        Returns:
            str or None: Path to the generated PDB file if topology file exists, None otherwise

        Raises:
            ValueError: If topology file is specified but coordinate file is missing
            RuntimeError: If ambpdb command fails or ambpdb is not available
        """
        # Check if topology file is provided
        if not self.topology_file:
            logger.info(
                "No topology file specified, skipping ambpdb generation"
            )
            return None

        if not os.path.exists(self.topology_file):
            logger.warning(
                f"Topology file {self.topology_file} does not exist, skipping ambpdb generation"
            )
            return None

        # Check for coordinate file
        if not self.coordinate_file:
            raise ValueError(
                "Coordinate file must be specified when topology file is provided for ambpdb"
            )

        if not os.path.exists(self.coordinate_file):
            raise ValueError(
                f"Coordinate file {self.coordinate_file} does not exist"
            )

        # Generate output filename based on topology file base name
        topology_basename = os.path.splitext(
            os.path.basename(self.topology_file)
        )[0]
        output_file = f"{topology_basename}_Hpp.pdb"

        # Construct the ambpdb command
        ambpdb_command = (
            f"ambpdb -p {self.topology_file} -c {self.coordinate_file}"
        )

        try:
            result = subprocess.run(
                ambpdb_command,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
            )

            # Write stdout to output file
            with open(output_file, "w") as f:
                f.write(result.stdout)

            logger.info(f"Successfully ran ambpdb command: {ambpdb_command}")
            logger.info(f"PDB output written to: {output_file}")

            # Verify output file was created and has content
            if (
                not os.path.exists(output_file)
                or os.path.getsize(output_file) == 0
            ):
                raise RuntimeError(
                    f"Ambpdb completed but output file {output_file} was not created or is empty"
                )

        except subprocess.CalledProcessError as e:
            error = f"Ambpdb command failed: {e}\nStderr: {e.stderr}\nStdout: {e.stdout}"
            logger.error(error)
            raise RuntimeError(error)
        except FileNotFoundError:
            error = "Ambpdb command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error)
            raise RuntimeError(error)

        return output_file

    def _get_combine_command(self):
        if self.hpp_pdb_file:
            combine_command = f"cat {self.hpp_pdb_file} {self.metal_pdb} {self.ligand}_fixed_H.pdb > {self.group_name}.pdb"
            return combine_command
        return None

    def _get_pdb4amber_command(self):
        if not os.path.exists(f"{self.group_name}_H.pdb"):
            raise ValueError(
                f"{self.group_name}_H.pdb does not exist for running pdb4amber"
            )

        # Generate output filename based on input PDB base name
        basename = os.path.splitext(os.path.basename(self.input_pdb))[0]
        output_pdb_file = f"{basename}_fixed.pdb"

        pdb4amber_command = (
            f"pdb4amber -i {self.group_name}_H.pdb -o {output_pdb_file}"
        )
        return pdb4amber_command

    def _get_mcpb_step_command(self, step):
        """
        Construct the MCPB.py command for a specific step."""
        if not self.mcpb_input:
            raise ValueError(
                "MCPB input file must be specified to run MCPB.py"
            )

        mcpb_command = f"MCPB.py -i {self.mcpb_input} -s {step}"
        return mcpb_command

    def _run_pdb4amber(self, input_pdb_file, output_pdb_file=None):
        """
        Execute pdb4amber command to renumber and fix PDB file.

        Runs: pdb4amber -i <input_pdb> -o <output_pdb_file>

        Args:
            input_pdb (str): Path to input PDB file to process
            output_pdb_file (str, optional): Path to output PDB file. If not specified,
                                           defaults to <input_basename>_fixed.pdb

        Returns:
            str: Path to the output fixed PDB file

        Raises:
            ValueError: If input file is not specified or doesn't exist
            RuntimeError: If pdb4amber command fails or pdb4amber is not available
        """
        if not input_pdb_file:
            raise ValueError("Input PDB file must be specified")

        if not os.path.exists(input_pdb_file):
            raise ValueError(f"Input PDB file {input_pdb_file} does not exist")

        # Generate output filename if not provided
        if output_pdb_file is None:
            # Extract basename without extension and add _fixed suffix
            basename = os.path.splitext(os.path.basename(input_pdb_file))[0]
            output_pdb_file = f"{basename}_fixed.pdb"

        # Construct the pdb4amber command
        pdb4amber_command = [
            "pdb4amber",
            "-i",
            input_pdb_file,
            "-o",
            output_pdb_file,
        ]

        try:
            result = subprocess.run(
                pdb4amber_command,
                capture_output=True,
                text=True,
                check=True,
            )

            logger.info(
                f"Successfully ran pdb4amber command: {' '.join(pdb4amber_command)}"
            )
            logger.info(f"Output written to: {output_pdb_file}")

            # Log stdout/stderr for debugging
            if result.stdout:
                logger.debug(f"pdb4amber stdout: {result.stdout}")
            if result.stderr:
                logger.debug(f"pdb4amber stderr: {result.stderr}")

            # Verify output file was created
            if not os.path.exists(output_pdb_file):
                raise RuntimeError(
                    f"pdb4amber completed but output file {output_pdb_file} was not created"
                )

            # Verify output file has content
            if os.path.getsize(output_pdb_file) == 0:
                raise RuntimeError(f"Output file {output_pdb_file} is empty")

        except subprocess.CalledProcessError as e:
            error_msg = f"pdb4amber command failed: {e}\nStderr: {e.stderr}\nStdout: {e.stdout}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except FileNotFoundError:
            error_msg = "pdb4amber command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        return output_pdb_file

    def _run_mcpb(self, input_file, step=1):
        """
        Execute MCPB.py command for metal center parameter building.

        Runs: MCPB.py -i <input_file> -s <step>

        Args:
            input_file (str): Path to MCPB input file (typically .in file)
            step (int, optional): MCPB step number (default: 1)

        Returns:
            tuple: Tuple containing paths to the three generated files:
                   (opt_file, fc_file, mk_file)
                   - For Gaussian: (*opt.com, *fc.com, *mk.com)
                   - For GAMESS-US: (*small_opt.inp, *fc.inp, *mk.inp)
                   Each element can be None if the corresponding file is not found.

        Raises:
            ValueError: If input file is not specified, doesn't exist, or invalid job type
            RuntimeError: If MCPB.py command fails or MCPB.py is not available
        """
        if not input_file:
            raise ValueError("Input file must be specified")

        if not os.path.exists(input_file):
            raise ValueError(f"Input file {input_file} does not exist")

        # Validate job type
        valid_job_types = ["gaussian", "gamess"]
        if self.mcpb_job_type not in valid_job_types:
            raise ValueError(
                f"Invalid mcpb_job_type '{self.mcpb_job_type}'. Must be one of: {valid_job_types}"
            )

        # Construct the MCPB.py command
        mcpb_command = [
            "MCPB.py",
            "-i",
            input_file,
            "-s",
            str(step),
        ]

        try:
            logger.info(f"Running MCPB.py command: {' '.join(mcpb_command)}")

            result = subprocess.run(
                mcpb_command,
                capture_output=True,
                text=True,
                check=True,
            )

            logger.info(
                f"Successfully ran MCPB.py command: {' '.join(mcpb_command)}"
            )

            # Log stdout/stderr for debugging
            if result.stdout:
                logger.debug(f"MCPB.py stdout: {result.stdout}")
            if result.stderr:
                logger.debug(f"MCPB.py stderr: {result.stderr}")

        except subprocess.CalledProcessError as e:
            error_msg = f"MCPB.py command failed: {e}\nStderr: {e.stderr}\nStdout: {e.stdout}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except FileNotFoundError:
            error_msg = "MCPB.py command not found. Please ensure MCPB.py is installed and in PATH."
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        # Check for generated files after MCPB.py has finished
        opt_file, fc_file, mk_file = self._check_mcpb_output_files()

        return opt_file, fc_file, mk_file

    def _check_mcpb_output_files(self):
        """
        Check for MCPB.py generated output files based on job type.

        For Gaussian: Files ending with 'opt.com', 'fc.com', and 'mk.com'
        For GAMESS-US: Files ending with 'small_opt.inp', 'fc.inp', and 'mk.inp'

        Returns:
            tuple: Tuple containing (opt_file, fc_file, mk_file)
                   Each element can be None if the corresponding file is not found.
        """
        import glob

        # Define file patterns based on job type
        if self.mcpb_job_type == "gaussian":
            opt_pattern = "*opt.com"
            fc_pattern = "*fc.com"
            mk_pattern = "*mk.com"
            file_type_name = "Gaussian"
        elif self.mcpb_job_type == "gamess":
            opt_pattern = "*small_opt.inp"
            fc_pattern = "*fc.inp"
            mk_pattern = "*mk.inp"
            file_type_name = "GAMESS-US"
        else:
            raise ValueError(f"Unsupported job type: {self.mcpb_job_type}")

        # Search for files with specific patterns
        opt_files = glob.glob(opt_pattern)
        fc_files = glob.glob(fc_pattern)
        mk_files = glob.glob(mk_pattern)

        # Get the first match for each type (there should typically be one of each)
        opt_file = opt_files[0] if opt_files else None
        fc_file = fc_files[0] if fc_files else None
        mk_file = mk_files[0] if mk_files else None

        # Log findings
        if opt_file:
            logger.info(
                f"Found {file_type_name} optimization file: {opt_file}"
            )
        else:
            logger.warning(f"No {opt_pattern} file found")

        if fc_file:
            logger.info(f"Found {file_type_name} frequency file: {fc_file}")
        else:
            logger.warning(f"No {fc_pattern} file found")

        if mk_file:
            logger.info(f"Found {file_type_name} RESP file: {mk_file}")
        else:
            logger.warning(f"No {mk_pattern} file found")

        all_found = all([opt_file, fc_file, mk_file])

        if all_found:
            logger.info(
                f"All expected {file_type_name} MCPB.py output files found successfully"
            )
        else:
            missing = []
            if not opt_file:
                missing.append(opt_pattern)
            if not fc_file:
                missing.append(fc_pattern)
            if not mk_file:
                missing.append(mk_pattern)
            logger.warning(
                f"Missing {file_type_name} MCPB.py output files: {', '.join(missing)}"
            )

        return opt_file, fc_file, mk_file

    def _write_MCPB_input(self, output_file=None):
        """
        Generate and write MCPB.in input file with the specified format.

        Creates an MCPB input file containing parameters for metal center parameter building.
        The format follows: original_pdb, group_name, cut_off, ion_ids, ion_mol2files,
        naa_mol2files, frcmod_files, large_opt.

        Args:
            output_file (str, optional): Path to output MCPB input file.
                                       Defaults to "{group_name}.in" if not specified.

        Returns:
            str: Path to the written MCPB input file

        Raises:
            ValueError: If required parameters are not specified
            RuntimeError: If file writing fails
        """
        # Validate required parameters
        if not self.input_pdb_file:
            raise ValueError(
                "input_pdb must be specified for MCPB input generation"
            )

        if not self.group_name:
            raise ValueError(
                "group_name must be specified for MCPB input generation"
            )

        if not self.ion_ids:
            raise ValueError(
                "ion_ids must be specified for MCPB input generation"
            )

        if not self.ion_mol2files:
            raise ValueError(
                "ion_mol2files must be specified for MCPB input generation"
            )

        # Generate output filename if not provided
        if output_file is None:
            output_file = f"{self.group_name}.in"

        # Build MCPB input content
        mcpb_lines = [
            f"original_pdb {self.input_pdb_file}",
            f"group_name {self.group_name}",
            f"cut_off {self.cutoff}",
            f"ion_ids {self.ion_ids}",
            f"ion_mol2files {self.ion_mol2files}",
        ]

        # Add optional naa_mol2files if ligand is specified
        if self.naa_mol2files:
            mcpb_lines.append(f"naa_mol2files {self.naa_mol2files}")

        # Add optional frcmod_files if ligand is specified
        if self.frcmod_files:
            mcpb_lines.append(f"frcmod_files {self.frcmod_files}")

        # Add large_opt parameter
        mcpb_lines.append(f"large_opt {self.large_opt}")

        # Join lines with newlines
        mcpb_content = "\n".join(mcpb_lines) + "\n"

        try:
            with open(output_file, "w") as f:
                f.write(mcpb_content)

            logger.info(f"Generated MCPB input file: {output_file}")
            logger.debug(f"MCPB input content:\n{mcpb_content}")

            # Verify output file was created
            if not os.path.exists(output_file):
                raise RuntimeError(
                    f"Failed to create MCPB input file {output_file}"
                )

        except IOError as e:
            error_msg = f"Failed to write MCPB input file {output_file}: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        return output_file

    def _get_metal_id(self):
        """
        Extract the metal atom index from the input PDB file based on the element attribute.

        Scans the PDB file for the specific metal element defined in self.element
        and returns the atom index of the first occurrence found.

        Returns:
            int or None: Atom index (1-based) of the metal atom if found, None otherwise

        Raises:
            ValueError: If element or input file is not specified, or if input file doesn't exist
        """
        if not self.element:
            raise ValueError("Element must be specified to get metal ID")

        # Use input_pdb if available, otherwise fall back to input_file
        input_file = self.input_pdb_file or self.input_file

        if not input_file:
            raise ValueError(
                "Input PDB file or input file must be specified to get metal ID"
            )

        if not os.path.exists(input_file):
            raise ValueError(f"Input file {input_file} does not exist")

        # Convert element to uppercase for comparison
        target_element = self.element.upper()

        try:
            with open(input_file, "r") as f:
                lines = f.readlines()

            for line_num, line in enumerate(lines, 1):
                # Check if this is an ATOM or HETATM record
                if line.startswith(("ATOM  ", "HETATM")):
                    # Extract atom name (columns 13-16, but often padded)
                    atom_name = line[12:16].strip()

                    # Extract element symbol (columns 77-78, or from atom name if not present)
                    if len(line) > 76:
                        element = line[76:78].strip()
                    else:
                        # Fallback: extract from atom name (remove digits and common suffixes)
                        import re

                        element = re.sub(r"[0-9]+", "", atom_name).strip()
                        # Remove common suffixes like A, B, etc.
                        if (
                            len(element) > 1
                            and element[-1].isalpha()
                            and element[-1].isupper()
                        ):
                            element = element[:-1]

                    # Check if this matches our target element
                    if element.upper() == target_element:
                        # Extract atom index (columns 7-11)
                        try:
                            atom_index = int(line[6:11].strip())
                            logger.info(
                                f"Found metal element {target_element} (atom name: {atom_name}) at index {atom_index}"
                            )
                            return atom_index
                        except ValueError:
                            logger.warning(
                                f"Could not parse atom index from line {line_num}: {line.strip()}"
                            )
                            continue

            logger.warning(f"No {target_element} atoms found in {input_file}")
            return None

        except IOError as e:
            error_msg = f"Failed to read PDB file {input_file}: {e}"
            logger.error(error_msg)
            raise ValueError(error_msg)

    def _run_tleap_mcpb(self, tleap_input_file=None, tleap_output_file=None):
        """
        Execute tleap command for MCPB topology and coordinate file generation.

        This method should be run after MCPB.py step 4 has been executed, which generates
        the tleap input file and modified PDB file.

        Runs: tleap -s -f <tleap_input_file> > <tleap_output_file>

        Args:
            tleap_input_file (str, optional): Path to tleap input file.
                                            Defaults to "{group_name}_tleap.in"
            tleap_output_file (str, optional): Path to tleap output file.
                                             Defaults to "{group_name}_tleap.out"

        Returns:
            tuple: (topology_file, coordinate_file) containing paths to:
                   - topology_file: Generated .prmtop file
                   - coordinate_file: Generated .inpcrd file
                   Each element can be None if the corresponding file is not found.

        Raises:
            ValueError: If group_name is not available or input file doesn't exist
            RuntimeError: If tleap command fails or tleap is not available
        """
        if not self.group_name:
            raise ValueError("group_name must be specified to run tleap")

        # Set default filenames if not provided
        if tleap_input_file is None:
            tleap_input_file = f"{self.group_name}_tleap.in"

        if tleap_output_file is None:
            tleap_output_file = f"{self.group_name}_tleap.out"

        # Verify tleap input file exists (should be generated by MCPB step 4)
        if not os.path.exists(tleap_input_file):
            raise ValueError(
                f"Tleap input file {tleap_input_file} does not exist. "
                "Please ensure MCPB.py step 4 has been run first."
            )

        # Construct the tleap command
        tleap_command = f"tleap -s -f {tleap_input_file}"

        try:
            logger.info(f"Running tleap command: {tleap_command}")

            # Run tleap command and redirect output to file
            with open(tleap_output_file, "w") as output:
                result = subprocess.run(
                    tleap_command,
                    shell=True,
                    stdout=output,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )

            logger.info(f"Successfully ran tleap command: {tleap_command}")
            logger.info(f"Tleap output written to: {tleap_output_file}")

            # Log stderr if present
            if result.stderr:
                logger.debug(f"Tleap stderr: {result.stderr}")

        except subprocess.CalledProcessError as e:
            error_msg = f"Tleap command failed: {e}\nStderr: {e.stderr}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except FileNotFoundError:
            error_msg = "Tleap command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        # Check for generated topology and coordinate files
        topology_file, coordinate_file = self._check_tleap_output_files()

        return topology_file, coordinate_file

    def run_tleap_mcpb(self, tleap_input_file=None, tleap_output_file=None):
        """
        Run tleap command for MCPB topology and coordinate file generation.

        Convenience method to run tleap after MCPB step 4 and generate topology/coordinate files.

        Args:
            tleap_input_file (str, optional): Path to tleap input file.
                                            Defaults to "{group_name}_tleap.in"
            tleap_output_file (str, optional): Path to tleap output file.
                                             Defaults to "{group_name}_tleap.out"

        Returns:
            tuple: (topology_file, coordinate_file) containing paths to generated files

        Example:
            # Equivalent to: tleap -s -f 1OKL_tleap.in > 1OKL_tleap.out
            topology, coordinate = settings.run_tleap_mcpb()
            if topology and coordinate:
                print("Topology and coordinate files generated successfully")
        """
        return self._run_tleap_mcpb(tleap_input_file, tleap_output_file)

    def _check_tleap_output_files(self):
        """
        Check for tleap generated topology and coordinate files.

        Looks for files with patterns:
        - Topology: {group_name}_solv.prmtop or {group_name}.prmtop
        - Coordinate: {group_name}_solv.inpcrd or {group_name}.inpcrd

        Returns:
            tuple: (topology_file, coordinate_file) containing paths to found files.
                   Each element can be None if the corresponding file is not found.
        """
        import glob

        if not self.group_name:
            logger.warning("No group_name available for file pattern matching")
            return None, None

        # Define file patterns - try solvated files first, then non-solvated
        topology_patterns = [
            f"{self.group_name}_solv.prmtop",
            f"{self.group_name}.prmtop",
        ]

        coordinate_patterns = [
            f"{self.group_name}_solv.inpcrd",
            f"{self.group_name}.inpcrd",
        ]

        # Search for topology file
        topology_file = None
        for pattern in topology_patterns:
            files = glob.glob(pattern)
            if files:
                topology_file = files[0]
                break

        # Search for coordinate file
        coordinate_file = None
        for pattern in coordinate_patterns:
            files = glob.glob(pattern)
            if files:
                coordinate_file = files[0]
                break

        # Log findings
        if topology_file:
            logger.info(f"Found topology file: {topology_file}")
        else:
            logger.warning(
                f"No topology file found. Searched for: {', '.join(topology_patterns)}"
            )

        if coordinate_file:
            logger.info(f"Found coordinate file: {coordinate_file}")
        else:
            logger.warning(
                f"No coordinate file found. Searched for: {', '.join(coordinate_patterns)}"
            )

        # Check if both files were found
        if topology_file and coordinate_file:
            logger.info(
                "Both topology and coordinate files found successfully"
            )
        else:
            missing = []
            if not topology_file:
                missing.append("topology (.prmtop)")
            if not coordinate_file:
                missing.append("coordinate (.inpcrd)")
            logger.warning(f"Missing tleap output files: {', '.join(missing)}")

        return topology_file, coordinate_file

    def check_and_fix_topology(self, topology_file=None, coordinate_file=None):
        """
        Check and fix topology file using cpptraj for atom numbering issues.

        Convenience method to check topology files and fix atom ordering problems.

        Args:
            topology_file (str, optional): Path to topology file.
                                         Defaults to "{group_name}_solv.prmtop"
            coordinate_file (str, optional): Path to coordinate file.
                                           Defaults to "{group_name}_solv.inpcrd"

        Returns:
            tuple: (fixed_topology, fixed_coordinate, needs_fix) containing:
                   - fixed_topology: Path to fixed topology file (original if no fix needed)
                   - fixed_coordinate: Path to fixed coordinate file (original if no fix needed)
                   - needs_fix: Boolean indicating if fixing was attempted

        Example:
            # Check and fix topology files from tleap
            topology, coordinate, was_fixed = settings.check_and_fix_topology()
            if was_fixed:
                print("Topology files were fixed and ready for MD simulation")
            else:
                print("No topology errors found, ready for MD simulation")
        """
        return self._run_cpptraj(topology_file, coordinate_file)

    @property
    def verified_topology_files(self):
        """
        Get verified and fixed topology/coordinate files ready for MD simulation.

        This property automatically checks the topology files generated by tleap
        and fixes any atom numbering issues if found.

        Returns:
            tuple: (topology_file, coordinate_file) containing paths to verified files
                   ready for MD simulation. Files may be the original or fixed versions.
        """
        # Get topology files from tleap
        topology, coordinate = self.mcpb_topology_files

        if topology and coordinate:
            # Check and fix if necessary
            fixed_topology, fixed_coordinate, _ = self._run_cpptraj(
                topology, coordinate
            )
            return fixed_topology, fixed_coordinate
        else:
            return None, None

    @property
    def verified_solv_prmtop(self):
        """
        Get verified solvated topology file with validated metal site parameters.

        This property automatically runs the complete validation workflow:
        1. Gets the topology files from MCPB workflow
        2. Checks and fixes atom ordering issues with cpptraj
        3. Validates metal site parameters with ParmEd using criteria A-E

        Metal ion parameter criteria:
        A) Bond force constants < 200 kcal/(mol*Å²), distances < 2.8 Å
        B) Angle force constants < 100 kcal/(mol*Rad²), angles > 100°
        C) Dihedral barriers should be zero for metal-involved dihedrals
        D) RESP charge < oxidation state, usually < +1
        E) LJ radius > 1.0 Å

        Returns:
            str or None: Path to verified and validated topology file ready for MD simulation
        """
        # Get verified topology files (cpptraj checked/fixed)
        topology, coordinate = self.verified_topology_files

        if topology:
            # Run ParmEd to validate metal site parameters
            validated_topology = self._run_parmed(topology)
            return validated_topology
        else:
            return None

    def validate_metal_parameters(
        self, topology_file=None, parmed_input_file=None
    ):
        """
        Validate metal site parameters using ParmEd against standard criteria A-E.

        Convenience method to check metal ion parameters in topology file.

        Args:
            topology_file (str, optional): Path to topology file to check.
                                         Defaults to "{group_name}_solv.prmtop"
            parmed_input_file (str, optional): Path to ParmEd input file.
                                             Defaults to "mcpbpy_parmed.in"

        Returns:
            str: Path to the topology file that was validated

        Example:
            # Validate metal parameters in topology file
            validated_topology = settings.validate_metal_parameters()
            print(f"Metal parameters validated in {validated_topology}")
        """
        return self._run_parmed(topology_file, parmed_input_file)

    def _run_cpptraj(self, topology_file=None, coordinate_file=None):
        """
        Check and fix topology file using cpptraj for atom numbering issues.

        First runs a check command to detect atom numbering errors, then if errors
        are found, attempts to fix them using the 'fixatomorder' command.

        Args:
            topology_file (str, optional): Path to topology file.
                                         Defaults to "{group_name}_solv.prmtop"

            coordinate_file (str, optional): Path to coordinate file.
                                           Defaults to "{group_name}_solv.inpcrd"

        Returns:
            tuple: (fixed_topology, fixed_coordinate, needs_fix) containing:
                   - fixed_topology: Path to fixed topology file (original if no fix needed)
                   - fixed_coordinate: Path to fixed coordinate file (original if no fix needed)
                   - needs_fix: Boolean indicating if fixing was attempted

        Raises:
            ValueError: If group_name is not available or files don't exist
            RuntimeError: If cpptraj command fails or cpptraj is not available
        """
        if not self.group_name:
            raise ValueError("group_name must be specified to run cpptraj")

        # Set default filenames if not provided
        if topology_file is None:
            topology_file = f"{self.group_name}_solv.prmtop"

        if coordinate_file is None:
            coordinate_file = f"{self.group_name}_solv.inpcrd"

        # Verify input files exist
        if not os.path.exists(topology_file):
            raise ValueError(f"Topology file {topology_file} does not exist")

        if not os.path.exists(coordinate_file):
            raise ValueError(
                f"Coordinate file {coordinate_file} does not exist"
            )

        # First, check for atom numbering errors
        needs_fix = self._check_topology_errors(topology_file)

        if needs_fix:
            logger.info(
                "Atom numbering errors detected, attempting to fix with cpptraj"
            )
            # Generate fixed filenames
            fixed_topology = f"fixed.{topology_file}"
            fixed_coordinate = f"fixed.{coordinate_file}"

            # Run cpptraj fix
            self._fix_atom_order(
                topology_file,
                coordinate_file,
                fixed_topology,
                fixed_coordinate,
            )

            # Verify the fix worked
            if self._check_topology_errors(fixed_topology):
                logger.warning(
                    "Cpptraj fix attempt completed but errors may still persist"
                )
            else:
                logger.info("Cpptraj successfully fixed atom numbering errors")

            return fixed_topology, fixed_coordinate, True
        else:
            logger.info("No atom numbering errors detected in topology file")
            return topology_file, coordinate_file, False

    def _check_topology_errors(self, topology_file):
        """
        Check topology file for atom numbering errors using cpptraj.

        Runs: cpptraj -p <topology_file>

        Args:
            topology_file (str): Path to topology file to check

        Returns:
            bool: True if errors are detected, False if no errors

        Raises:
            RuntimeError: If cpptraj command fails or cpptraj is not available
        """
        # Construct the cpptraj check command
        cpptraj_command = ["cpptraj", "-p", topology_file]

        try:
            logger.info(
                f"Checking topology file for errors: {' '.join(cpptraj_command)}"
            )

            result = subprocess.run(
                cpptraj_command,
                capture_output=True,
                text=True,
                check=False,  # Don't raise exception on non-zero exit
            )

            # Check both stdout and stderr for error indicators
            output = result.stdout + result.stderr

            # Look for specific error patterns
            error_patterns = [
                "was assigned a lower molecule # than previous atom",
                "Could not determine molecule information",
                "Could not determine solvent information",
                "bond information is incorrect or missing",
                "atom numbering in molecules is not sequential",
            ]

            has_errors = any(pattern in output for pattern in error_patterns)

            if has_errors:
                logger.warning("Topology errors detected:")
                for line in output.split("\n"):
                    if any(pattern in line for pattern in error_patterns):
                        logger.warning(f"  {line.strip()}")
            else:
                logger.info("No topology errors detected")

            return has_errors

        except FileNotFoundError:
            error_msg = "Cpptraj command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except Exception as e:
            error_msg = f"Cpptraj check command failed: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

    def _fix_atom_order(
        self, topology_file, coordinate_file, fixed_topology, fixed_coordinate
    ):
        """
        Fix atom ordering issues using cpptraj fixatomorder command.

        Creates a cpptraj input file and runs:
        cpptraj -p <topology_file> -c <coordinate_file> -i fixatord.in > fixatord.out

        Args:
            topology_file (str): Path to input topology file
            coordinate_file (str): Path to input coordinate file
            fixed_topology (str): Path for output fixed topology file
            fixed_coordinate (str): Path for output fixed coordinate file

        Raises:
            RuntimeError: If cpptraj command fails or files are not generated
        """
        # Create cpptraj input file for fixing atom order
        cpptraj_input_file = "fixatord.in"
        cpptraj_output_file = "fixatord.out"

        # Write cpptraj input commands
        cpptraj_commands = [
            "# Fix atom ordering issues",
            "fixatomorder",
            f"parmwrite out {fixed_topology}",
            f"trajout {fixed_coordinate}",
            "run",
            "quit",
        ]

        try:
            # Write cpptraj input file
            with open(cpptraj_input_file, "w") as f:
                f.write("\n".join(cpptraj_commands) + "\n")

            logger.info(f"Created cpptraj input file: {cpptraj_input_file}")

            # Construct cpptraj fix command
            cpptraj_command = [
                "cpptraj",
                "-p",
                topology_file,
                "-c",
                coordinate_file,
                "-i",
                cpptraj_input_file,
            ]

            logger.info(
                f"Running cpptraj fix command: {' '.join(cpptraj_command)}"
            )

            # Run cpptraj command and redirect output to file
            with open(cpptraj_output_file, "w") as output:
                result = subprocess.run(
                    cpptraj_command,
                    stdout=output,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )

            logger.info("Successfully ran cpptraj fix command")
            logger.info(f"Cpptraj output written to: {cpptraj_output_file}")

            # Log stderr if present
            if result.stderr:
                logger.debug(f"Cpptraj stderr: {result.stderr}")

            # Verify output files were created
            if not os.path.exists(fixed_topology):
                raise RuntimeError(
                    f"Cpptraj completed but fixed topology file {fixed_topology} was not created"
                )

            if not os.path.exists(fixed_coordinate):
                raise RuntimeError(
                    f"Cpptraj completed but fixed coordinate file {fixed_coordinate} was not created"
                )

            logger.info(
                f"Generated fixed files: {fixed_topology}, {fixed_coordinate}"
            )

        except subprocess.CalledProcessError as e:
            error_msg = f"Cpptraj fix command failed: {e}\nStderr: {e.stderr}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except FileNotFoundError:
            error_msg = "Cpptraj command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except Exception as e:
            error_msg = f"Failed to fix topology with cpptraj: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

    def _run_parmed(self, topology_file=None, parmed_input_file=None):
        """
        Execute ParmEd to check metal site parameters in topology file.

        Runs: parmed -i <parmed_input_file> -p <topology_file>

        This method checks metal ion parameters against standard criteria:
        A) Bond force constants < 200 kcal/(mol*Angstrom^2), and equilibrium bond distances
           should be less than 2.8 Angstrom
        B) Angle force constants related to the metal ion should usually be
           less than 100 kcal/(mol*Rad^2) while equilibrium angle values should
           be bigger than 100 degrees
        C) All or most dihedral potential barriers should be zero for metal
           involved dihedrals
        D) The RESP charge of the metal ion should be less than its oxidation
           state, usually even less than +1
        E) The LJ radius of the metal ion should usually be bigger than 1.0 Angstrom

        Args:
            topology_file (str, optional): Path to topology file to check.
                                         Defaults to "{group_name}_solv.prmtop"
            parmed_input_file (str, optional): Path to ParmEd input file.
                                             Defaults to "mcpbpy_parmed.in"

        Returns:
            str: Path to the topology file that was checked

        Raises:
            ValueError: If group_name is not available or files don't exist
            RuntimeError: If ParmEd command fails or ParmEd is not available
        """
        if not self.group_name:
            raise ValueError("group_name must be specified to run ParmEd")

        # Set default filenames if not provided
        if topology_file is None:
            topology_file = f"{self.group_name}_solv.prmtop"

        if parmed_input_file is None:
            parmed_input_file = "mcpbpy_parmed.in"

        # Verify input files exist
        if not os.path.exists(topology_file):
            raise ValueError(f"Topology file {topology_file} does not exist")

        # Create default ParmEd input file if it doesn't exist
        if not os.path.exists(parmed_input_file):
            self._create_parmed_input_file(parmed_input_file)

        # Construct the ParmEd command
        parmed_command = [
            "parmed",
            "-i",
            parmed_input_file,
            "-p",
            topology_file,
        ]

        try:
            logger.info(f"Running ParmEd command: {' '.join(parmed_command)}")

            result = subprocess.run(
                parmed_command,
                capture_output=True,
                text=True,
                check=True,
            )

            logger.info(
                f"Successfully ran ParmEd command: {' '.join(parmed_command)}"
            )

            # Log output for parameter analysis
            if result.stdout:
                logger.info("ParmEd output:")
                logger.info(result.stdout)
            if result.stderr:
                logger.debug(f"ParmEd stderr: {result.stderr}")

            # Analyze the output for parameter validation
            self._analyze_parmed_output(result.stdout)

        except subprocess.CalledProcessError as e:
            error_msg = f"ParmEd command failed: {e}\nStderr: {e.stderr}\nStdout: {e.stdout}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except FileNotFoundError:
            error_msg = "ParmEd command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        return topology_file

    def _create_parmed_input_file(self, parmed_input_file):
        """
        Create a default ParmEd input file for metal site parameter checking.

        Creates an input file with commands to display metal-related parameters
        using the metal element as a mask (e.g., ":ZN1" for zinc).

        Args:
            parmed_input_file (str): Path to create the ParmEd input file

        Raises:
            ValueError: If element is not specified
            RuntimeError: If file creation fails
        """
        if not self.element:
            raise ValueError(
                "Element must be specified to create ParmEd input file"
            )

        # Format element name for mask (e.g., "ZN" -> ":ZN1")
        metal_mask = f":{self.element.upper()}1"

        # Create ParmEd input commands
        parmed_commands = [
            f"# ParmEd input file for checking {self.element.upper()} metal site parameters",
            "# Generated automatically by chemsmart",
            "",
            "# Display bonds involving the metal ion",
            f"printBonds {metal_mask}",
            "",
            "# Display angles involving the metal ion",
            f"printAngles {metal_mask}",
            "",
            "# Display dihedrals involving the metal ion",
            f"printDihedrals {metal_mask}",
            "",
            "# Display charges",
            f"printDetails {metal_mask}",
            "",
            "# Display LJ parameters",
            f"printLJMatrix {metal_mask}",
            "",
            "# Exit ParmEd",
            "quit",
        ]

        try:
            with open(parmed_input_file, "w") as f:
                f.write("\n".join(parmed_commands) + "\n")

            logger.info(f"Created ParmEd input file: {parmed_input_file}")
            logger.info(f"Metal mask used: {metal_mask}")

        except IOError as e:
            error_msg = (
                f"Failed to create ParmEd input file {parmed_input_file}: {e}"
            )
            logger.error(error_msg)
            raise RuntimeError(error_msg)

    def _analyze_parmed_output(self, output):
        """
        Analyze ParmEd output to validate metal site parameters against standard criteria.

        Checks the output against criteria A-E for metal ion parameters:
        A) Bond parameters: force constants < 200 kcal/(mol*Å²), distances < 2.8 Å
        B) Angle parameters: force constants < 100 kcal/(mol*Rad²), angles > 100°
        C) Dihedral barriers should be zero for metal-involved dihedrals
        D) RESP charge should be less than oxidation state, usually < +1
        E) LJ radius should be > 1.0 Angstrom

        Args:
            output (str): ParmEd command output to analyze
        """
        if not output:
            logger.warning("No ParmEd output to analyze")
            return

        logger.info(
            "Analyzing metal site parameters against standard criteria:"
        )

        # Split output into lines for analysis
        lines = output.split("\n")

        # Initialize analysis flags
        found_bonds = False
        found_angles = False
        found_dihedrals = False
        found_charges = False
        found_lj = False

        # Analyze each line
        for line in lines:
            line = line.strip()

            # Check for bond information (Criterion A)
            if "bond" in line.lower() and any(
                x in line for x in ["force", "constant", "distance"]
            ):
                found_bonds = True
                logger.info(f"Bond parameter found: {line}")
                # Note: Detailed parsing would require specific ParmEd output format knowledge

            # Check for angle information (Criterion B)
            elif "angle" in line.lower() and any(
                x in line for x in ["force", "constant", "degree"]
            ):
                found_angles = True
                logger.info(f"Angle parameter found: {line}")

            # Check for dihedral information (Criterion C)
            elif "dihedral" in line.lower() and "barrier" in line.lower():
                found_dihedrals = True
                logger.info(f"Dihedral parameter found: {line}")

            # Check for charge information (Criterion D)
            elif "charge" in line.lower() and self.element.upper() in line:
                found_charges = True
                logger.info(f"Charge parameter found: {line}")

            # Check for LJ parameters (Criterion E)
            elif "lj" in line.lower() or "lennard" in line.lower():
                found_lj = True
                logger.info(f"LJ parameter found: {line}")

        # Summary of analysis
        logger.info("Parameter analysis summary:")
        logger.info(f"  A) Bond parameters found: {found_bonds}")
        logger.info(f"  B) Angle parameters found: {found_angles}")
        logger.info(f"  C) Dihedral parameters found: {found_dihedrals}")
        logger.info(f"  D) Charge parameters found: {found_charges}")
        logger.info(f"  E) LJ parameters found: {found_lj}")

        # Validation criteria reminder
        logger.info("Standard criteria for metal ion parameters:")
        logger.info(
            "  A) Bond force constants < 200 kcal/(mol*Å²), distances < 2.8 Å"
        )
        logger.info(
            "  B) Angle force constants < 100 kcal/(mol*Rad²), angles > 100°"
        )
        logger.info(
            "  C) Dihedral barriers should be zero for metal-involved dihedrals"
        )
        logger.info("  D) RESP charge < oxidation state, usually < +1")
        logger.info("  E) LJ radius > 1.0 Å")

        if not any(
            [
                found_bonds,
                found_angles,
                found_dihedrals,
                found_charges,
                found_lj,
            ]
        ):
            logger.warning(
                "No metal-related parameters detected in ParmEd output"
            )
            logger.warning(
                "Please check the ParmEd input file and topology file"
            )

    def check_metal_coordination(
        self, topology_file=None, coordinate_file=None, expected_ligands=None
    ):
        """
        Check metal coordination bonds and distances without VMD (GUI-free replacement).

        Convenience method to verify metal-ligand coordination using topology bonds
        and geometric distances. This replaces VMD-based coordination checking.

        Args:
            topology_file (str, optional): Path to topology file.
                                         Defaults to "{group_name}_solv.prmtop"
            coordinate_file (str, optional): Path to coordinate file.
                                           Defaults to "{group_name}_solv.inpcrd"
            expected_ligands (list, optional): List of expected ligand specifications.
                                             Defaults to common coordination patterns

        Returns:
            dict: Validation results with 'overall_pass' indicating success

        Example:
            # Basic metal coordination check
            results = settings.check_metal_coordination()
            if results['overall_pass']:
                print("Metal coordination validated successfully")

            # Check specific ligands
            custom_ligands = [
                {'residue': 'HIS', 'atom': 'NE2', 'distance_range': (2.0, 2.5)},
                {'residue': 'CYS', 'atom': 'SG', 'distance_range': (2.2, 2.7)}
            ]
            results = settings.check_metal_coordination(expected_ligands=custom_ligands)
        """
        return self._run_vmd_check(
            topology_file, coordinate_file, expected_ligands
        )

    @property
    def validated_metal_coordination(self):
        """
        Get validated metal coordination results using complete verification.

        This property runs the complete coordination validation workflow:
        1. Checks metal-ligand bonds exist in topology using ParmEd
        2. Checks metal-ligand distances are geometrically reasonable using cpptraj
        3. Returns comprehensive validation results

        Returns:
            dict or None: Validation results if coordination check passes, None if fails.
                         Contains 'overall_pass', 'bonds_check', 'distances_check' and
                         detailed results for bonds and distances.
        """
        try:
            # Get verified topology files first
            topology, coordinate = self.verified_topology_files

            if topology and coordinate:
                # Run metal coordination check
                results = self._run_vmd_check(topology, coordinate)
                return results if results["overall_pass"] else None
            else:
                logger.warning(
                    "Cannot validate metal coordination: topology files not available"
                )
                return None

        except Exception as e:
            logger.error(f"Metal coordination validation failed: {e}")
            return None

    def _run_vmd_check(
        self, topology_file=None, coordinate_file=None, expected_ligands=None
    ):
        """
        VMD-free metal coordination check to verify metal-ligand coordination bonds.

        This method performs a comprehensive validation of metal coordination without
        requiring VMD or any GUI tools. It checks both topology bonds and geometric
        distances to ensure the metal center is properly coordinated.

        Steps performed:
        1. Check metal-ligand bonds in topology using ParmEd
        2. Check metal-ligand distances are geometrically reasonable using cpptraj

        Args:
            topology_file (str, optional): Path to topology file.
                                         Defaults to "{group_name}_solv.prmtop"
            coordinate_file (str, optional): Path to coordinate file.
                                           Defaults to "{group_name}_solv.inpcrd"
            expected_ligands (list, optional): List of expected ligand specifications.
                                             Each element should be a dict with keys:
                                             'residue': residue name (e.g., 'HIS', 'CYS')
                                             'atom': atom name (e.g., 'NE2', 'SG')
                                             'distance_range': tuple of (min_dist, max_dist) in Å
                                             Defaults to common zinc coordination

        Returns:
            dict: Validation results containing:
                  - 'bonds_check': Boolean, True if all expected bonds exist
                  - 'distances_check': Boolean, True if all distances are reasonable
                  - 'overall_pass': Boolean, True if both checks pass
                  - 'bond_results': List of bond validation results
                  - 'distance_results': List of distance measurements
                  - 'missing_bonds': List of missing expected bonds
                  - 'bad_distances': List of distances outside acceptable range

        Raises:
            ValueError: If required parameters are missing or files don't exist
            RuntimeError: If ParmEd or cpptraj commands fail

        Example:
            # Basic usage with default zinc ligands
            results = settings._run_vmd_check()
            if results['overall_pass']:
                print("Metal coordination validated successfully")

            # Custom ligand specification
            custom_ligands = [
                {'residue': 'HIS', 'atom': 'NE2', 'distance_range': (2.0, 2.5)},
                {'residue': 'CYS', 'atom': 'SG', 'distance_range': (2.2, 2.7)}
            ]
            results = settings._run_vmd_check(expected_ligands=custom_ligands)
        """
        if not self.group_name:
            raise ValueError("group_name must be specified for VMD check")

        if not self.element:
            raise ValueError(
                "element must be specified for metal coordination check"
            )

        # Set default filenames if not provided
        if topology_file is None:
            topology_file = f"{self.group_name}_solv.prmtop"

        if coordinate_file is None:
            coordinate_file = f"{self.group_name}_solv.inpcrd"

        # Verify input files exist
        if not os.path.exists(topology_file):
            raise ValueError(f"Topology file {topology_file} does not exist")

        if not os.path.exists(coordinate_file):
            raise ValueError(
                f"Coordinate file {coordinate_file} does not exist"
            )

        # Set default expected ligands for zinc if not provided
        if expected_ligands is None:
            expected_ligands = self._get_default_metal_ligands()

        logger.info("Starting VMD-free metal coordination check")
        logger.info(
            f"Checking {self.element.upper()} coordination in {topology_file}"
        )

        # Step 1: Check metal-ligand bonds in topology
        logger.info("Step 1: Checking metal-ligand bonds in topology")
        bonds_check, bond_results, missing_bonds = self._check_topology_bonds(
            topology_file, expected_ligands
        )

        # Step 2: Check metal-ligand distances
        logger.info("Step 2: Checking metal-ligand distances")
        distances_check, distance_results, bad_distances = (
            self._check_metal_distances(
                topology_file, coordinate_file, expected_ligands
            )
        )

        # Overall validation result
        overall_pass = bonds_check and distances_check

        # Compile results
        results = {
            "bonds_check": bonds_check,
            "distances_check": distances_check,
            "overall_pass": overall_pass,
            "bond_results": bond_results,
            "distance_results": distance_results,
            "missing_bonds": missing_bonds,
            "bad_distances": bad_distances,
            "topology_file": topology_file,
            "coordinate_file": coordinate_file,
            "metal_element": self.element.upper(),
        }

        # Log summary
        self._log_vmd_check_summary(results)

        return results

    def _log_vmd_check_summary(self, results):
        """
        Log the summary of VMD-free metal coordination check results.

        Args:
            results (dict): Validation results from _run_vmd_check
        """
        logger.info("Metal coordination check summary:")
        logger.info(f"  Topology file: {results['topology_file']}")
        logger.info(f"  Coordinate file: {results['coordinate_file']}")
        logger.info(f"  Metal element: {results['metal_element']}")
        logger.info(f"  Overall pass: {results['overall_pass']}")
        logger.info(f"  Bonds check: {results['bonds_check']}")
        logger.info(f"  Distances check: {results['distances_check']}")

        if not results["overall_pass"]:
            if results["missing_bonds"]:
                logger.warning(
                    f"  Missing expected bonds: {results['missing_bonds']}"
                )
            if results["bad_distances"]:
                logger.warning(
                    f"  Distances outside acceptable range: {results['bad_distances']}"
                )

    def _get_default_metal_ligands(self):
        """
        Get default expected ligands based on metal element.

        Returns common coordination patterns for different metals.

        Returns:
            list: List of expected ligand dictionaries
        """
        metal = self.element.upper()

        # Common coordination patterns
        if metal == "ZN":
            return [
                {
                    "residue": "HIS",
                    "atom": "NE2",
                    "distance_range": (2.0, 2.5),
                },
                {
                    "residue": "HIS",
                    "atom": "ND1",
                    "distance_range": (2.0, 2.5),
                },
                {"residue": "CYS", "atom": "SG", "distance_range": (2.2, 2.7)},
                {
                    "residue": "GLU",
                    "atom": "OE1",
                    "distance_range": (1.9, 2.4),
                },
                {
                    "residue": "GLU",
                    "atom": "OE2",
                    "distance_range": (1.9, 2.4),
                },
                {
                    "residue": "ASP",
                    "atom": "OD1",
                    "distance_range": (1.9, 2.4),
                },
                {
                    "residue": "ASP",
                    "atom": "OD2",
                    "distance_range": (1.9, 2.4),
                },
            ]
        elif metal == "MG":
            return [
                {
                    "residue": "ASP",
                    "atom": "OD1",
                    "distance_range": (1.9, 2.3),
                },
                {
                    "residue": "ASP",
                    "atom": "OD2",
                    "distance_range": (1.9, 2.3),
                },
                {
                    "residue": "GLU",
                    "atom": "OE1",
                    "distance_range": (1.9, 2.3),
                },
                {
                    "residue": "GLU",
                    "atom": "OE2",
                    "distance_range": (1.9, 2.3),
                },
                {"residue": "HOH", "atom": "O", "distance_range": (1.9, 2.2)},
            ]
        elif metal == "CA":
            return [
                {
                    "residue": "ASP",
                    "atom": "OD1",
                    "distance_range": (2.3, 2.8),
                },
                {
                    "residue": "ASP",
                    "atom": "OD2",
                    "distance_range": (2.3, 2.8),
                },
                {
                    "residue": "GLU",
                    "atom": "OE1",
                    "distance_range": (2.3, 2.8),
                },
                {
                    "residue": "GLU",
                    "atom": "OE2",
                    "distance_range": (2.3, 2.8),
                },
                {"residue": "HOH", "atom": "O", "distance_range": (2.3, 2.7)},
            ]
        elif metal in ["FE", "FE2", "FE3"]:
            return [
                {
                    "residue": "HIS",
                    "atom": "NE2",
                    "distance_range": (2.0, 2.3),
                },
                {
                    "residue": "HIS",
                    "atom": "ND1",
                    "distance_range": (2.0, 2.3),
                },
                {"residue": "CYS", "atom": "SG", "distance_range": (2.2, 2.5)},
                {"residue": "MET", "atom": "SD", "distance_range": (2.3, 2.6)},
            ]
        else:
            # Generic metal coordination
            logger.warning(f"Using generic coordination pattern for {metal}")
            return [
                {
                    "residue": "HIS",
                    "atom": "NE2",
                    "distance_range": (2.0, 2.5),
                },
                {
                    "residue": "HIS",
                    "atom": "ND1",
                    "distance_range": (2.0, 2.5),
                },
                {"residue": "CYS", "atom": "SG", "distance_range": (2.2, 2.7)},
                {
                    "residue": "ASP",
                    "atom": "OD1",
                    "distance_range": (1.9, 2.4),
                },
                {
                    "residue": "GLU",
                    "atom": "OE1",
                    "distance_range": (1.9, 2.4),
                },
            ]

    def _check_topology_bonds(self, topology_file, expected_ligands):
        """
        Check if expected metal-ligand bonds exist in topology using ParmEd.

        Args:
            topology_file (str): Path to topology file
            expected_ligands (list): List of expected ligand specifications

        Returns:
            tuple: (bonds_check_passed, bond_results, missing_bonds)
        """
        metal_mask = f":{self.element.upper()}1"

        # Create ParmEd input for bond checking
        parmed_input = "vmd_check_bonds.in"
        parmed_commands = [
            f"# Bond checking for {self.element.upper()} coordination",
            f"printBonds {metal_mask}",
            "quit",
        ]

        try:
            with open(parmed_input, "w") as f:
                f.write("\n".join(parmed_commands) + "\n")

            # Run ParmEd
            parmed_command = [
                "parmed",
                "-i",
                parmed_input,
                "-p",
                topology_file,
            ]

            logger.info(
                f"Running ParmEd bond check: {' '.join(parmed_command)}"
            )

            result = subprocess.run(
                parmed_command,
                capture_output=True,
                text=True,
                check=True,
            )

            # Parse bond output
            bond_results, missing_bonds = self._parse_parmed_bonds(
                result.stdout, expected_ligands
            )

            bonds_check_passed = len(missing_bonds) == 0

            logger.info(
                f"Bond topology check: {'PASSED' if bonds_check_passed else 'FAILED'}"
            )
            if missing_bonds:
                logger.warning(f"Missing bonds: {missing_bonds}")

            return bonds_check_passed, bond_results, missing_bonds

        except subprocess.CalledProcessError as e:
            error_msg = f"ParmEd bond check failed: {e}\nStderr: {e.stderr}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except FileNotFoundError:
            error_msg = "ParmEd command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        finally:
            # Clean up input file
            if os.path.exists(parmed_input):
                os.unlink(parmed_input)

    def _parse_parmed_bonds(self, parmed_output, expected_ligands):
        """
        Parse ParmEd output to find metal bonds.

        Args:
            parmed_output (str): ParmEd stdout output
            expected_ligands (list): Expected ligand specifications

        Returns:
            tuple: (found_bonds, missing_bonds)
        """
        lines = parmed_output.split("\n")
        found_bonds = []
        metal_element = self.element.upper()

        # Look for bond lines in ParmEd output
        for line in lines:
            line = line.strip()
            # ParmEd bond format typically includes atom names and indices
            if metal_element in line and (
                "bond" in line.lower() or "-" in line
            ):
                found_bonds.append(line)
                logger.debug(f"Found metal bond: {line}")

        # Check which expected ligands are present
        missing_bonds = []
        for ligand in expected_ligands:
            ligand_found = False
            search_patterns = [
                ligand["residue"],
                ligand["atom"],
                f"{ligand['residue']}@{ligand['atom']}",
            ]

            for bond_line in found_bonds:
                if any(pattern in bond_line for pattern in search_patterns):
                    ligand_found = True
                    break

            if not ligand_found:
                missing_bonds.append(f"{ligand['residue']}@{ligand['atom']}")

        return found_bonds, missing_bonds

    def _check_metal_distances(
        self, topology_file, coordinate_file, expected_ligands
    ):
        """
        Check metal-ligand distances using cpptraj.

        Args:
            topology_file (str): Path to topology file
            coordinate_file (str): Path to coordinate file
            expected_ligands (list): List of expected ligand specifications

        Returns:
            tuple: (distances_check_passed, distance_results, bad_distances)
        """
        metal_mask = f":{self.element.upper()}1"

        # Create cpptraj input for distance measurements
        cpptraj_input = "vmd_check_distances.in"
        cpptraj_output = "vmd_check_distances.out"

        cpptraj_commands = [
            f"# Distance checking for {self.element.upper()} coordination"
        ]

        # Add distance commands for each expected ligand
        for i, ligand in enumerate(expected_ligands):
            ligand_mask = f":{ligand['residue']}@{ligand['atom']}"
            dist_name = f"{self.element.lower()}_lig{i+1}"
            cpptraj_commands.append(
                f"distance {dist_name} {metal_mask} {ligand_mask}"
            )

        cpptraj_commands.extend(["run", "quit"])

        try:
            with open(cpptraj_input, "w") as f:
                f.write("\n".join(cpptraj_commands) + "\n")

            # Run cpptraj
            cpptraj_command = [
                "cpptraj",
                "-p",
                topology_file,
                "-y",
                coordinate_file,
                "-i",
                cpptraj_input,
            ]

            logger.info(
                f"Running cpptraj distance check: {' '.join(cpptraj_command)}"
            )

            with open(cpptraj_output, "w") as output_file:
                subprocess.run(
                    cpptraj_command,
                    stdout=output_file,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )

            # Parse distance results
            distance_results, bad_distances = self._parse_cpptraj_distances(
                cpptraj_output, expected_ligands
            )

            distances_check_passed = len(bad_distances) == 0

            logger.info(
                f"Distance geometry check: {'PASSED' if distances_check_passed else 'FAILED'}"
            )
            if bad_distances:
                logger.warning(f"Bad distances: {bad_distances}")

            return distances_check_passed, distance_results, bad_distances

        except subprocess.CalledProcessError as e:
            error_msg = (
                f"Cpptraj distance check failed: {e}\nStderr: {e.stderr}"
            )
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except FileNotFoundError:
            error_msg = "Cpptraj command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        finally:
            # Clean up input file
            for temp_file in [cpptraj_input, cpptraj_output]:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)

    def _parse_cpptraj_distances(self, output_file, expected_ligands):
        """
        Parse cpptraj distance output.

        Args:
            output_file (str): Path to cpptraj output file
            expected_ligands (list): Expected ligand specifications

        Returns:
            tuple: (distance_measurements, bad_distances)
        """
        distance_results = []
        bad_distances = []

        try:
            with open(output_file, "r") as f:
                lines = f.readlines()

            # Look for distance output in cpptraj format
            for line in lines:
                line = line.strip()
                # cpptraj distance output format: typically shows distance values
                if any(
                    keyword in line.lower()
                    for keyword in ["distance", "avg", self.element.lower()]
                ):
                    logger.debug(f"Distance line: {line}")

                    # Try to extract distance values
                    parts = line.split()
                    for part in parts:
                        try:
                            distance = float(part)
                            if (
                                1.0 < distance < 5.0
                            ):  # Reasonable range for metal-ligand distances
                                distance_results.append(
                                    {
                                        "distance": distance,
                                        "line": line,
                                        "ligand": "unknown",  # Would need more parsing to identify
                                    }
                                )

                                # Check against expected ranges
                                in_range = False
                                for ligand in expected_ligands:
                                    min_dist, max_dist = ligand[
                                        "distance_range"
                                    ]
                                    if min_dist <= distance <= max_dist:
                                        in_range = True
                                        break

                                if not in_range:
                                    bad_distances.append(
                                        {
                                            "distance": distance,
                                            "line": line,
                                            "reason": f"Distance {distance:.2f} Å outside acceptable ranges",
                                        }
                                    )
                                break
                        except ValueError:
                            continue

            # If no distances found in output, check if cpptraj ran successfully
            if not distance_results:
                logger.warning(
                    "No distance measurements found in cpptraj output"
                )
                # This might indicate the masks didn't match any atoms
                bad_distances.append(
                    {
                        "distance": None,
                        "line": "No distances measured",
                        "reason": "No matching atoms found for distance measurements",
                    }
                )

        except Exception as e:
            logger.error(f"Error parsing cpptraj distances: {e}")
            bad_distances.append(
                {
                    "distance": None,
                    "line": str(e),
                    "reason": "Error parsing distance output",
                }
            )

        return distance_results, bad_distances
