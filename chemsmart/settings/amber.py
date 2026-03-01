"""
Amber project settings and configuration.

This module provides project-level settings for Amber molecular dynamics
calculations, including default parameters, job configurations, and
YAML-based project management.
"""

import logging
from typing import Optional

from chemsmart.jobs.amber.settings import AmberJobSettings
from chemsmart.utils.mixins import RegistryMixin

logger = logging.getLogger(__name__)
project_settings_registry: list[str] = []


class AmberProjectSettings(RegistryMixin):
    """
    Base class for Amber project settings with default configurations.

    Provides fundamental settings and configurations for Amber molecular
    dynamics calculations. Includes default values for computational
    parameters and methods for creating job-specific settings.

    Attributes:
        PROJECT_NAME (str): Name identifier for the project.
        default_ensemble (str): Default statistical ensemble.
        default_temperature (float): Default temperature in Kelvin.
        default_steps (int): Default number of simulation steps.
        default_timestep (float): Default timestep in ps.
    """

    PROJECT_NAME = "general"
    default_ensemble = "NVT"
    default_temperature = 300.0
    default_steps = 50000
    default_timestep = 0.002

    def main_settings(self):
        """
        Create main Amber job settings with default parameters.

        Returns:
            AmberJobSettings: Default settings for Amber calculations
        """
        return AmberJobSettings(
            ensemble=self.default_ensemble,
            temperature=self.default_temperature,
            steps=self.default_steps,
            simulation_type="md",
        )

    def md_settings(self):
        """
        Create Amber molecular dynamics settings.

        Returns:
            AmberJobSettings: Settings optimized for MD simulations
        """
        settings = self.main_settings()
        settings.simulation_type = "md"
        settings.steps = self.default_steps
        settings.ensemble = self.default_ensemble
        settings.temperature = self.default_temperature
        return settings

    def energy_settings(self):
        """
        Create Amber energy calculation settings.

        Returns:
            AmberJobSettings: Settings for single-point energy calculations
        """
        settings = self.main_settings()
        settings.simulation_type = "energy"
        settings.steps = 0  # No dynamics
        settings.ensemble = "NVE"  # No thermostat needed
        return settings

    def opt_settings(self):
        """
        Create Amber optimization settings.

        Returns:
            AmberJobSettings: Settings for geometry optimization
        """
        settings = self.main_settings()
        settings.simulation_type = "opt"
        settings.steps = 1000  # Optimization cycles
        settings.ensemble = "NVE"  # No thermostat for optimization
        return settings

    @classmethod
    def from_project(cls, project):
        """
        Create project settings from a project directory.

        Args:
            project: Project directory path or name

        Returns:
            AmberProjectSettings: Configured settings instance
        """
        # Try to find YAML settings file
        import os

        if os.path.isdir(project):
            # Look for amber.yml in project directory
            amber_yaml = os.path.join(project, "amber.yml")
            if os.path.exists(amber_yaml):
                try:
                    import yaml

                    with open(amber_yaml, "r") as f:
                        yaml_data = yaml.safe_load(f)
                    return YamlAmberProjectSettings.from_dict(yaml_data)
                except Exception as e:
                    logger.warning(f"Failed to load YAML settings: {e}")

        # Fallback to default settings
        logger.debug(f"Using default Amber project settings for: {project}")
        return cls()


class YamlAmberProjectSettings(AmberProjectSettings):
    """
    YAML-based Amber project settings.

    This class reads Amber settings from YAML files and provides
    job-specific configurations based on the loaded parameters.
    """

    def __init__(
        self,
        md_settings: Optional[AmberJobSettings] = None,
        energy_settings: Optional[AmberJobSettings] = None,
        opt_settings: Optional[AmberJobSettings] = None,
        **kwargs,
    ):
        """
        Initialize YAML-based Amber project settings.

        Args:
            md_settings: Settings for MD simulations
            energy_settings: Settings for energy calculations
            opt_settings: Settings for optimization
            **kwargs: Additional settings
        """
        super().__init__()

        self._md_settings = md_settings or self.main_settings()
        self._energy_settings = energy_settings or self.energy_settings()
        self._opt_settings = opt_settings or self.opt_settings()

        # Store additional configuration
        for key, value in kwargs.items():
            setattr(self, key, value)

    def md_settings(self):
        """
        Get molecular dynamics calculation settings.

        Returns:
            AmberJobSettings: MD simulation settings
        """
        return self._md_settings

    def energy_settings(self):
        """
        Get energy calculation settings.

        Returns:
            AmberJobSettings: Energy calculation settings
        """
        return self._energy_settings

    def opt_settings(self):
        """
        Get optimization calculation settings.

        Returns:
            AmberJobSettings: Optimization settings
        """
        return self._opt_settings

    @classmethod
    def from_dict(cls, config_dict):
        """
        Create settings from dictionary (typically from YAML).

        Args:
            config_dict: Dictionary containing configuration

        Returns:
            YamlAmberProjectSettings: Configured settings
        """
        # Extract different job type configurations
        md_config = config_dict.get("md", {})
        energy_config = config_dict.get("energy", {})
        opt_config = config_dict.get("opt", {})

        # Create settings objects
        md_settings = AmberJobSettings(**md_config) if md_config else None
        energy_settings = (
            AmberJobSettings(**energy_config) if energy_config else None
        )
        opt_settings = AmberJobSettings(**opt_config) if opt_config else None

        return cls(
            md_settings=md_settings,
            energy_settings=energy_settings,
            opt_settings=opt_settings,
            **config_dict.get("general", {}),
        )


# Register the project settings
project_settings_registry.append("AmberProjectSettings")
