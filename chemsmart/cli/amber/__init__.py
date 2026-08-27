"""
Amber Command Line Interface package.

This package provides the complete CLI interface for Amber molecular dynamics
calculations including MD simulations, energy calculations, and geometry
optimization.
"""

from chemsmart.cli.amber.amber import amber

__all__ = [
    "amber",
]
