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

    Properties:
        ligand_pdb (str): Path to extracted ligand PDB file
        metal_pdb (str): Path to extracted metal PDB file
        reduced_ligand (str): Path to ligand PDB with added hydrogens
        ligand_mol2_file (str): Path to ligand mol2 file
        ligand_frcmod_file (str): Path to ligand frcmod file
        water_mol2_file (str): Path to water mol2 file
        hpp_pdb_file (str or None): Path to PDB file with AMBER naming scheme
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

    def combine_pdb_files(self, input_files, output_file):
        """
        Combine multiple PDB files into a single PDB file.

        Convenience method to combine PDB files, equivalent to concatenation.

        Args:
            input_files (list): List of input PDB file paths to combine
            output_file (str): Path to the output combined PDB file

        Returns:
            str: Path to the combined PDB file

        Example:
            # Equivalent to: cat 1OKL_Hpp_fixed.pdb ZN.pdb MNS_fixed_H.pdb > 1OKL_H.pdb
            combined_file = settings.combine_pdb_files(
                ['1OKL_Hpp_fixed.pdb', 'ZN.pdb', 'MNS_fixed_H.pdb'],
                '1OKL_H.pdb'
            )
        """
        return self._combine(input_files, output_file)

    def fix_pdb_file(self, input_pdb_file, output_pdb_file=None):
        """
        Fix and renumber PDB file using pdb4amber.

        Convenience method to run pdb4amber for PDB file processing.

        Args:
            input_pdb_file (str): Path to input PDB file to process
            output_pdb_file (str, optional): Path to output PDB file. If not specified,
                                           defaults to <input_basename>_fixed.pdb

        Returns:
            str: Path to the output fixed PDB file

        Example:
            # Equivalent to: pdb4amber -i 1OKL_H.pdb -o 1OKL_fixed_H.pdb
            fixed_file = settings.fix_pdb_file('1OKL_H.pdb', '1OKL_fixed_H.pdb')
        """
        return self._run_pdb4amber(input_pdb_file, output_pdb_file)

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

    def _combine(self, input_files, output_file):
        """
        Combine multiple PDB files into a single PDB file.

        Equivalent to: cat file1.pdb file2.pdb file3.pdb > combined.pdb

        Args:
            input_files (list): List of input PDB file paths to combine
            output_file (str): Path to the output combined PDB file

        Returns:
            str: Path to the combined PDB file

        Raises:
            ValueError: If input_files is empty or contains invalid files
            RuntimeError: If file operations fail
        """
        if not input_files:
            raise ValueError("Input files list cannot be empty")

        if not isinstance(input_files, (list, tuple)):
            raise ValueError("Input files must be provided as a list or tuple")

        # Validate input files exist
        missing_files = []
        for file_path in input_files:
            if not os.path.exists(file_path):
                missing_files.append(file_path)

        if missing_files:
            raise ValueError(
                f"The following input files do not exist: {missing_files}"
            )

        try:
            # Combine files by reading each and writing to output
            with open(output_file, "w") as output:
                for i, file_path in enumerate(input_files):
                    logger.info(
                        f"Adding content from {file_path} to {output_file}"
                    )

                    with open(file_path, "r") as input_file:
                        content = input_file.read()

                        # Write content to output file
                        output.write(content)

                        # Add newline between files if the current file doesn't end with one
                        # and it's not the last file
                        if (
                            not content.endswith("\n")
                            and i < len(input_files) - 1
                        ):
                            output.write("\n")

            # Verify output file was created and has content
            if not os.path.exists(output_file):
                raise RuntimeError(
                    f"Failed to create output file {output_file}"
                )

            if os.path.getsize(output_file) == 0:
                raise RuntimeError(f"Output file {output_file} is empty")

            logger.info(
                f"Successfully combined {len(input_files)} files into {output_file}"
            )

            # Log file sizes for verification
            total_size = sum(os.path.getsize(f) for f in input_files)
            output_size = os.path.getsize(output_file)
            logger.info(
                f"Input files total size: {total_size} bytes, output size: {output_size} bytes"
            )

        except IOError as e:
            error_msg = f"Failed to combine PDB files: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)
        except Exception as e:
            error_msg = f"Unexpected error while combining files: {e}"
            logger.error(error_msg)
            raise RuntimeError(error_msg)

        return output_file

    def _run_pdb4amber(self, input_pdb_file, output_pdb_file=None):
        """
        Execute pdb4amber command to renumber and fix PDB file.

        Runs: pdb4amber -i <input_pdb_file> -o <output_pdb_file>

        Args:
            input_pdb_file (str): Path to input PDB file to process
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
