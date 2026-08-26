from chemsmart.cli.gaussian.gaussian import gaussian
from chemsmart.cli.reaction import register_reaction_cli
from chemsmart.jobs.gaussian.reaction import GaussianReactionJob

reaction = register_reaction_cli(gaussian, GaussianReactionJob)
