from chemsmart.cli.orca.orca import orca
from chemsmart.cli.reaction import register_reaction_cli
from chemsmart.jobs.orca.reaction import ORCAReactionJob

reaction = register_reaction_cli(orca, ORCAReactionJob)
