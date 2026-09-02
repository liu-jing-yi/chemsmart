import os
from dataclasses import dataclass

import yaml

from chemsmart.settings.user import CHEMSMARTUserSettings

user_settings = CHEMSMARTUserSettings()


@dataclass(frozen=True)
class ChainStep:
    program: str
    job: str


class ChainProjectSettings:
    """YAML chain-project aliases and optional pipeline steps.

    Top-level keys ``crest``, ``xtb``, ``gaussian``, and ``orca`` are aliases
    to existing per-program project names. ``steps`` is an optional list of
    ``{program, job}`` mappings. Aliases-only files are valid.
    """

    PROGRAMS = ("crest", "xtb", "gaussian", "orca")
    PROJECT_NAME = "chain"

    def __init__(self, aliases, steps, project_name):
        self._aliases = dict(aliases)
        self._steps = tuple(steps)
        self.PROJECT_NAME = project_name

    @property
    def steps(self):
        return self._steps

    def project_for(self, program):
        """Return the per-program project alias for ``program``.

        Args:
            program (str): One of ``crest``, ``xtb``, ``gaussian``, ``orca``.

        Returns:
            str: Project name stem, e.g. ``test`` for
            ``~/.chemsmart/gaussian/test.yaml``.

        Raises:
            ValueError: If ``program`` is unknown or has no alias in this
                chain file.
        """
        if program not in self.PROGRAMS:
            allowed = ", ".join(self.PROGRAMS)
            raise ValueError(
                f"Unknown chain program {program!r}. "
                f"Allowed programs: {allowed}."
            )
        if program not in self._aliases:
            raise ValueError(
                f"Chain project {self.PROJECT_NAME!r} has no "
                f"{program} section."
            )
        return self._aliases[program]

    @classmethod
    def from_project(cls, project):
        # Match Gaussian/xTB precedence: user projects override packaged
        # templates, and missing names fail with a discoverable config path.
        user_project_settings = cls._from_user_project_name(project)
        if user_project_settings is not None:
            return user_project_settings

        packaged_project_settings = cls._from_packaged_project_name(project)
        if packaged_project_settings is not None:
            return packaged_project_settings

        templates_path = os.path.join(os.path.dirname(__file__), "templates")
        raise FileNotFoundError(
            f"No chain project settings implemented for {project}.\n\n"
            f"Place new chain project settings .yaml file in "
            f"{user_settings.user_chain_settings_dir}.\n\n"
            f"Templates for such settings.yaml files are available at "
            f"{templates_path}\n\n "
            f"Currently available projects: "
            f"{user_settings.all_available_chain_projects}"
        )

    @classmethod
    def from_yaml(cls, filename):
        filename = os.path.abspath(filename)
        with open(filename) as handle:
            config = yaml.safe_load(handle)
        if config is None:
            config = {}
        if not isinstance(config, dict):
            raise ValueError(
                f"Chain project settings in {filename} must be a mapping."
            )
        project_name = os.path.basename(filename).removesuffix(".yaml")
        return cls._from_config(config, project_name=project_name)

    @classmethod
    def _from_config(cls, config, project_name):
        allowed_keys = set(cls.PROGRAMS) | {"steps"}
        unknown_keys = sorted(set(config) - allowed_keys)
        if unknown_keys:
            allowed = ", ".join(sorted(allowed_keys))
            unknown = ", ".join(unknown_keys)
            raise ValueError(
                f"Unknown chain project keys: {unknown}. "
                f"Allowed keys: {allowed}."
            )
        aliases = cls._parse_aliases(config)
        steps = cls._parse_steps(config, aliases)
        return cls(aliases=aliases, steps=steps, project_name=project_name)

    @classmethod
    def _parse_aliases(cls, config):
        aliases = {}
        for program in cls.PROGRAMS:
            if program not in config:
                continue
            value = config[program]
            if value is None:
                continue
            if not isinstance(value, str) or not value.strip():
                raise ValueError(
                    f"Chain project {program} alias must be a "
                    f"project name string."
                )
            aliases[program] = value.strip()
        return aliases

    @classmethod
    def _parse_steps(cls, config, aliases):
        if "steps" not in config or config["steps"] is None:
            return ()
        raw_steps = config["steps"]
        if not isinstance(raw_steps, list):
            raise ValueError("Chain project 'steps' must be a list.")
        steps = []
        allowed_step_keys = {"program", "job"}
        for index, raw_step in enumerate(raw_steps):
            if not isinstance(raw_step, dict):
                raise ValueError(
                    f"Chain project step {index} must be a mapping."
                )
            unknown_keys = sorted(set(raw_step) - allowed_step_keys)
            if unknown_keys:
                unknown = ", ".join(unknown_keys)
                raise ValueError(f"Unknown chain step keys: {unknown}.")
            program = raw_step.get("program")
            job = raw_step.get("job")
            if not isinstance(program, str) or not program.strip():
                raise ValueError(
                    "Each chain step requires a 'program' string."
                )
            if not isinstance(job, str) or not job.strip():
                raise ValueError("Each chain step requires a 'job' string.")
            program = program.strip()
            job = job.strip()
            if program not in cls.PROGRAMS:
                allowed = ", ".join(cls.PROGRAMS)
                raise ValueError(
                    f"Unknown chain step program {program!r}. "
                    f"Allowed programs: {allowed}."
                )
            if program not in aliases:
                raise ValueError(
                    f"Chain step program {program!r} has no project alias."
                )
            steps.append(ChainStep(program=program, job=job))
        return tuple(steps)

    @classmethod
    def _from_user_project_name(cls, project_name):
        if project_name is None:
            return None
        project_name_yaml_path = os.path.join(
            CHEMSMARTUserSettings().user_chain_settings_dir,
            f"{project_name}.yaml",
        )
        try:
            return cls.from_yaml(project_name_yaml_path)
        except FileNotFoundError:
            return None

    @classmethod
    def _from_packaged_project_name(cls, project_name):
        if project_name is None:
            return None
        current_file_dir = os.path.dirname(os.path.abspath(__file__))
        project_name_yaml_path = os.path.join(
            current_file_dir,
            "templates",
            ".chemsmart",
            "chain",
            f"{project_name}.yaml",
        )
        try:
            return cls.from_yaml(project_name_yaml_path)
        except FileNotFoundError:
            return None
