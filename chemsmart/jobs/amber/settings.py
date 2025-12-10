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
    """

    def __init__(
        self,
        element: str = None,
        ligand: str = None,
        ligand_charge: int = None,
        ligand_multiplicity: int = None,
        charge_method: str = "bcc",
        gaussian_logfile: str = None,
        gaussian_fchkfile: str = None,
        parameter: str = None,
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
            **kwargs: Additional arguments passed to parent class
        """
        super().__init__(**kwargs)

        self.element = element
        self.ligand = ligand
        self.ligand_charge = ligand_charge
        self.ligand_multiplicity = ligand_multiplicity
        self.charge_method = charge_method
        self.parameter = parameter
        self.gaussian_logfile = gaussian_logfile
        self.gaussian_fchkfile = gaussian_fchkfile

    @property
    def ligand_pdb(self):
        return self._get_ligand_pdb()

    @property
    def metal_pdb(self):
        return self._get_metal_pdb()

    @property
    def reduced_ligand(self):
        """
        Run the reduce command to add hydrogens to the ligand PDB file.
        """
        return self._run_reduce()

    @property
    def ligand_mol2_file(self):
        """
        Get the ligand mol2 file name based on the ligand name.

        Returns:
            str: Ligand mol2 file name
        """

        return self._run_antechamber()

    @property
    def ligand_frcmod_file(self):
        """
        Get the ligand frcmod file name based on the ligand name.

        Returns:
            str: Ligand frcmod file name
        """

        return self._run_parmchk2()

    def _get_ligand_pdb(self):
        """
        Extract ligand lines from a PDB input file and write them to <ligand>.pdb

        Returns:
            str: Path to the created ligand PDB file
        """
        if not self.ligand:
            raise ValueError("Ligand name must be specified")

        if not self.input_file:
            raise ValueError("Input file must be specified")

        if not os.path.exists(self.input_file):
            raise ValueError(f"Input file {self.input_file} does not exist")

        # Format ligand name )
        ligand_upper = self._residue_formatting(self.ligand)
        ligand_pdb_filename = f"{self.ligand}.pdb"
        with open(self.input_file, "r") as f:
            lines = f.readlines()

        ligand_lines = [line for line in lines if self.ligand in line]

        if not ligand_lines:
            raise ValueError(
                f"No lines containing ligand '{self.ligand}' found in {self.input_file}"
            )

        # Ensure ligand names are capitalized
        processed_lines = []
        for line in ligand_lines:
            processed_line = line.replace(self.ligand, ligand_upper)
            processed_lines.append(processed_line)

        with open(ligand_pdb_filename, "w") as out:
            out.writelines(processed_lines)

        return ligand_pdb_filename

    def _get_metal_pdb(self):
        """
        Extract metal lines from a PDB input file and write them to <element>.pdb
        using pure Python (no grep/subprocess).

        Returns:
            str: Path to the created metal PDB file
        """
        if not self.element:
            raise ValueError("Element name must be specified")

        if not self.input_file:
            raise ValueError("Input file must be specified")

        if not os.path.exists(self.input_file):
            raise ValueError(f"Input file {self.input_file} does not exist")

        # Format element name using _residue_formatting (includes warning if needed)
        element_upper = self._residue_formatting(self.element)

        metal_pdb_filename = f"{self.element}.pdb"

        # Read and filter lines containing element name
        with open(self.input_file, "r") as f:
            lines = f.readlines()

        metal_lines = [line for line in lines if self.element in line]

        if not metal_lines:
            raise ValueError(
                f"No lines containing element '{self.element}' found in {self.input_file}"
            )

        # Process lines to ensure element names are capitalized
        processed_lines = []
        for line in metal_lines:
            # Replace element name with uppercase version in the line
            processed_line = line.replace(self.element, element_upper)
            processed_lines.append(processed_line)

        # Write processed lines to new PDB file
        with open(metal_pdb_filename, "w") as out:
            out.writelines(processed_lines)

        return metal_pdb_filename

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

    def _run_reduce(self):
        """
        Execute the reduce command to add hydrogens to ligand PDB file.

        Runs: reduce {ligand}.pdb > {ligand}_H.pdb

        Returns:
            str: Path to the output file with added hydrogens

        Raises:
            ValueError: If ligand is not specified
            RuntimeError: If reduce command fails or reduce is not available
        """

        if not self.ligand:
            raise ValueError("Ligand name must be specified to run reduce")

        ligand_pdb_file = self.ligand_pdb
        output_file = f"{self.ligand}_H.pdb"
        reduce_command = f"reduce {ligand_pdb_file}"

        try:
            result = subprocess.run(
                reduce_command,
                shell=True,
                capture_output=True,
                text=True,
                check=True,
            )
            with open(output_file, "w") as f:
                f.write(result.stdout)

            logger.info(f"Successfully ran reduce command: {reduce_command}")
            logger.info(f"Output written to: {output_file}")

        except subprocess.CalledProcessError as e:
            error = f"Reduce command failed: {e}\nStderr: {e.stderr}"
            logger.error(error)
            raise RuntimeError(error)
        except FileNotFoundError:
            error = "Reduce command not found. Please ensure reduce is installed and in PATH."
            logger.error(error)
            raise RuntimeError(error)

        return output_file

    def _run_antechamber(self):
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

        try:
            subprocess.run(
                antechamber_command, capture_output=True, text=True, check=True
            )

            logger.info(
                f"Successfully ran antechamber command: {' '.join(antechamber_command)}"
            )
            logger.info(f"Output written to: {output_file}")

            # Verify output file was created
            if not os.path.exists(output_file):
                raise RuntimeError(
                    f"Antechamber completed but output file {output_file} was not created"
                )

        except subprocess.CalledProcessError as e:
            error = f"Antechamber command failed: {e}\nStderr: {e.stderr}\nStdout: {e.stdout}"
            logger.error(error)
            raise RuntimeError(error)
        except FileNotFoundError:
            error = "Antechamber command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error)
            raise RuntimeError(error)

        return output_file

    def _run_parmchk2(self):
        """
        Execute the parmchk2 command to generate force field modification file.

        First renames {ligand}_pre.mol2 to {ligand}.mol2, then runs:
        parmchk2 -i {ligand}.mol2 -o {ligand}.frcmod -f mol2

        Returns:
            str: Path to the output frcmod file

        Raises:
            ValueError: If ligand is not specified
            RuntimeError: If parmchk2 command fails or parmchk2 is not available
        """
        if not self.ligand:
            raise ValueError("Ligand name must be specified to run parmchk2")

        pre_mol2_file = self.ligand_mol2_file

        # Rename from _pre.mol2 to .mol2
        mol2_file = self._rename_mol2_file(pre_mol2_file)

        output_file = f"{self.ligand}.frcmod"

        # Construct the parmchk2 command
        parmchk2_command = [
            "parmchk2",
            "-i",
            mol2_file,
            "-o",
            output_file,
            "-f",
            "mol2",
        ]

        try:
            subprocess.run(
                parmchk2_command, capture_output=True, text=True, check=True
            )

            logger.info(
                f"Successfully ran parmchk2 command: {' '.join(parmchk2_command)}"
            )
            logger.info(f"Output written to: {output_file}")

            # Verify output file was created
            if not os.path.exists(output_file):
                raise RuntimeError(
                    f"Parmchk2 completed but output file {output_file} was not created"
                )

        except subprocess.CalledProcessError as e:
            error = f"Parmchk2 command failed: {e}\nStderr: {e.stderr}\nStdout: {e.stdout}"
            logger.error(error)
            raise RuntimeError(error)
        except FileNotFoundError:
            error = "Parmchk2 command not found. Please ensure AmberTools is installed and in PATH."
            logger.error(error)
            raise RuntimeError(error)

        return output_file

    def _rename_mol2_file(self, pre_mol2_file):
        """
        Rename {ligand}_pre.mol2 to {ligand}.mol2.

        This is needed because parmchk2 expects the standard naming convention
        without the "_pre" suffix.

        Args:
            pre_mol2_file (str): Path to the {ligand}_pre.mol2 file

        Returns:
            str: Path to the renamed {ligand}.mol2 file

        Raises:
            ValueError: If ligand is not specified
            RuntimeError: If file operations fail
        """
        if not self.ligand:
            raise ValueError(
                "Ligand name must be specified to rename mol2 file"
            )

        if not os.path.exists(pre_mol2_file):
            raise RuntimeError(f"Source file {pre_mol2_file} does not exist")

        # Target filename
        mol2_file = f"{self.ligand}.mol2"

        try:
            # Copy the file to the new name
            import shutil

            shutil.copy2(pre_mol2_file, mol2_file)

            logger.info(f"Renamed {pre_mol2_file} to {mol2_file}")

            # Verify the new file exists
            if not os.path.exists(mol2_file):
                raise RuntimeError(f"Failed to create {mol2_file}")

        except Exception as e:
            error_msg = f"Failed to rename {pre_mol2_file} to {mol2_file}: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        return mol2_file
