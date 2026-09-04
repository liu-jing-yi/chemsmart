###################
 Chain CLI Options
###################

This page documents the CLI options for the ``chain`` command group. Use ``chemsmart sub chain --help`` for the complete
list, including the ``custom`` subcommand.

*************************
 Basic Command Structure
*************************

**Custom pipeline** (no nested subcommand, or ``custom``):

.. code:: bash

   chemsmart sub [OPTIONS] chain [CHAIN_OPTIONS] -f FILE -s STEP [-s ...] [custom]

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
      -  Input structure file (required for the CLI pipeline)

   -  -  ``-s, --steps``

      -  string (repeatable)

      -  One pipeline step as ``PROGRAM:JOB`` or quoted ``PROGRAM: [OPTIONS] JOB`` (e.g. ``gaussian:opt`` or
         ``"gaussian: -o maxstep=8,maxsize=12 opt"``). Required for the CLI pipeline; repeat for each phase. Not the
         same as ``run``/``sub`` ``-s/--server``.

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
   -  For the CLI pipeline, ``-f`` and at least one ``-s/--steps`` are required. The ``custom`` subcommand runs the same
      pipeline.
   -  ``-l`` and ``-a`` are mutually exclusive.

Execution Control
=================

The chain group accepts the standard job options documented in :doc:`cli-overview`:

-  ``-S/-R, --skip-completed/--no-skip-completed``
-  ``--fake/--no-fake``
-  ``--scratch/--no-scratch``

**********
 Examples
**********

Full pipeline on HPC:

.. code:: bash

   chemsmart sub -s mycluster chain -p combined -f ligand.xyz -c 0 -m 1 \
     -s crest:conformers -s xtb:opt -s gaussian:opt -s orca:sp

Quoted step with Gaussian opt extras:

.. code:: bash

   chemsmart run chain -p combined -f mol.xyz -c 0 -m 1 \
     -s "gaussian: -o maxstep=8,maxsize=12 opt" \
     -s gaussian:sp

Same pipeline via the ``custom`` subcommand:

.. code:: bash

   chemsmart sub chain -p combined -f ligand.xyz -c 0 -m 1 \
     -s gaussian:opt custom

Ordinary Gaussian opt or pKa uses that program's own project YAML:

.. code:: bash

   chemsmart sub gaussian -p gaussian_project2 -f ligand.xyz -c 0 -m 1 -o maxstep=8,maxsize=12 opt
   chemsmart sub gaussian -p gaussian_project2 -f ligand.xyz -c 0 -m 1 pka submit

************
 Next Steps
************

-  :doc:`chain-jobs` — chain YAML aliases, ``-s/--steps`` pipeline, geometry handoff, and HPC submission
-  :doc:`configuration-project-settings` — chain project file location and alias rules
