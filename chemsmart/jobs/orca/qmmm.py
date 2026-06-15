from chemsmart.jobs.orca.job import ORCAJob


class ORCAQMMMJob(ORCAJob):
    TYPE = "orcaqmmm"

    def __init__(
        self, molecule, settings, label, structure_filename=None, **kwargs
    ):
        super().__init__(
            molecule=molecule, settings=settings, label=label, **kwargs
        )
        self.structure_filename = structure_filename
