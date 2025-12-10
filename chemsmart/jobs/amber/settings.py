"""
Amber job settings and configuration classes.

This module provides settings classes for Amber molecular dynamics and
energy calculations, following the same pattern as other computational
chemistry modules in chemsmart.
"""

import logging
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
        email_notifications: bool = True,
        email: str = "",
        simulation_type: str = "md",
        input_file: Optional[str] = None,
        topology_file: Optional[str] = None,
        coordinate_file: Optional[str] = None,
        steps: int = 50000,
        temperature: float = 300.0,
        pressure: Optional[float] = None,
        ensemble: str = "NVT",
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
            email_notifications: Enable email notifications
            email: Email address for notifications
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
        self.email_notifications = email_notifications
        self.email = email
        self.simulation_type = simulation_type
        self.input_file = input_file
        self.topology_file = topology_file
        self.coordinate_file = coordinate_file
        self.steps = steps
        self.temperature = temperature
        self.pressure = pressure
        self.ensemble = ensemble

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
