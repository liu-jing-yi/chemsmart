###################
 Chain CLI Options
###################

This page documents the CLI options for the ``chain`` command group. Use ``chemsmart sub chain --help`` for the complete
list, including nested ``gaussian``, ``orca``, ``xtb``, and ``crest`` subcommands.

*************************
 Basic Command Structure
*************************

**YAML pipeline** (no nested subcommand):

.. code:: bash

   chemsmart sub [OPTIONS] chain [CHAIN_OPTIONS] -f FILE

**Nested program slice**:

.. code:: bash

   chemsmart sub [OPTIONS] chain [CHAIN_OPTIONS] <PROGRAM> [PROGRAM_OPTIONS] <SUBCMD> [SUBCMD_OPTIONS]

where ``<PROGRAM>`` is one of ``gaussian``, ``orca``, ``xtb``, or ``crest``.

***************
 Chain Options
***************

Project and File Options
========================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-p, --project``
      -  string
      -  Chain project settings from ``~/.chemsmart/chain/*.yaml`` (required)

   -  -  ``-f, --filename``
      -  string
      -  Input structure file (required for the YAML pipeline; nested programs also accept ``-f``)

   -  -  ``-l, --label``
      -  string
      -  Custom job label (without extension)

   -  -  ``-a, --append-label``
      -  string
      -  String appended to the base filename for the label

   -  -  ``-i, --index``
      -  string
      -  Structure index (1-based; default: last structure)

   -  -  ``-c, --charge``
      -  int
      -  Molecular charge (pipeline; merged into each child job)

   -  -  ``-m, --multiplicity``
      -  int
      -  Spin multiplicity (pipeline; merged into each child job)

.. note::

   -  ``-p`` uses the chain project name without the ``.yaml`` extension.
   -  For the YAML pipeline, ``-f`` is required. Aliases-only chain files (no ``steps``) require a nested program
      command instead.
   -  ``-l`` and ``-a`` are mutually exclusive.

Execution Control
=================

The chain group accepts the standard job options documented in :doc:`cli-overview`:

-  ``-S/-R, --skip-completed/--no-skip-completed``
-  ``--fake/--no-fake``
-  ``--scratch/--no-scratch``

Nested Program Options
======================

When a nested ``gaussian``, ``orca``, ``xtb``, or ``crest`` command is used, that program's options apply — including
its own ``-p`` override, ``-f``, charge/multiplicity, and subcommand-specific flags. See:

-  :doc:`gaussian-cli-options`
-  :doc:`orca-cli-options`
-  :doc:`xtb-cli-options`
-  :doc:`crest-cli-options`

**********
 Examples
**********

Full pipeline on HPC:

.. code:: bash

   chemsmart sub -s mycluster chain -p combined -f ligand.xyz -c 0 -m 1

Gaussian pKa slice using the chain file's Gaussian alias:

.. code:: bash

   chemsmart sub chain -p combined gaussian -f ligand.xyz -c 0 -m 1 pka submit

Override the Gaussian project while keeping other chain aliases for other nested commands:

.. code:: bash

   chemsmart sub chain -p combined gaussian -p gas_solv -f ligand.xyz opt

************
 Next Steps
************

-  :doc:`chain-jobs` — chain YAML schema, supported pipeline steps, geometry handoff, and HPC submission
-  :doc:`configuration-project-settings` — chain project file location and alias rules
