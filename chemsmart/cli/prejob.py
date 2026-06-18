"""
Compute-node preparation hooks run before the main chemsmart run pipeline.

Hooks are invoked from generated ``chemsmart_run_*.py`` scripts on the cluster
before ``run(cli_args)`` recreates and executes the job.
"""

from __future__ import annotations

import logging
from collections.abc import Callable, Sequence

logger = logging.getLogger(__name__)

PrejobHook = Callable[[Sequence[str]], None]


def run_prejob_hooks(cli_args: Sequence[str]) -> None:
    """
    Run registered pre-job hooks before the chemsmart run pipeline.

    Args:
        cli_args: Command-line arguments passed to ``chemsmart run``.
    """
    for hook in _registered_hooks():
        logger.debug("Running prejob hook: %s", hook.__name__)
        hook(cli_args)


def _registered_hooks() -> tuple[PrejobHook, ...]:
    """Return hook callables in execution order."""
    return ()
