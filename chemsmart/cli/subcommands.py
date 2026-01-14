from chemsmart.cli.amber import amber
from chemsmart.cli.gaussian import gaussian
from chemsmart.cli.mol import mol
from chemsmart.cli.nciplot import nciplot
from chemsmart.cli.thermochemistry import thermochemistry

subcommands = [
    gaussian,
    amber,
    mol,
    nciplot,
    thermochemistry,
]
