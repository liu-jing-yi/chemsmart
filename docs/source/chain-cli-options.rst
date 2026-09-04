###################
 Chain CLI Options
###################

This page documents the CLI options for the ``chain`` command group. Use ``chemsmart sub chain --help`` for the complete
list, including ``custom`` and the workflow subcommands ``pka``, ``fukui``, ``redox``, and ``reaction``.

*************************
 Basic Command Structure
*************************

**Custom pipeline** (no nested subcommand, or ``custom``):

.. code:: bash

   chemsmart sub [OPTIONS] chain [CHAIN_OPTIONS] -f FILE -s STEP [-s ...] [custom]

**Workflow submit** (``pka``, ``fukui``, ``redox``, ``reaction``):

.. code:: bash

   chemsmart sub [OPTIONS] chain -p PROJECT -f FILE -c CHARGE -m MULT \
     pka|fukui|redox|reaction --program {gaussian,orca} [WORKFLOW_OPTIONS]

**Workflow analyze** (pKa, Fukui, and redox only; no ``-p``, ``-f``, or ``--program``):

.. code:: bash

   chemsmart run chain pka|fukui|redox analyze [ANALYZE_OPTIONS]

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
      -  Chain project settings from ``~/.chemsmart/chain/*.yaml``. Required for the ``-s/--steps`` pipeline and for
         ``pka`` / ``fukui`` / ``redox`` / ``reaction`` **submit**. Not required for analyze.

   -  -  ``-f, --filename``
      -  string
      -  Input structure file (required for the CLI pipeline and for workflow submit)

   -  -  ``-s, --steps``

      -  string (repeatable)

      -  One pipeline step as ``PROGRAM:JOB`` or quoted ``PROGRAM: [OPTIONS] JOB`` (e.g. ``gaussian:opt`` or
         ``"gaussian: -o maxstep=8,maxsize=12 opt"``). Required for the CLI pipeline; repeat for each phase. Not the
         same as ``run``/``sub`` ``-s/--server``. Cannot be combined with ``pka``, ``fukui``, ``redox``, or
         ``reaction``.

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
      -  Molecular charge (pipeline and workflow submit; merged into each child job)

   -  -  ``-m, --multiplicity``
      -  int
      -  Spin multiplicity (pipeline and workflow submit; merged into each child job)

Workflow Submit Options
=======================

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``--program``
      -  ``gaussian`` or ``orca``
      -  QC program for workflow **submit**. Required. Theory and solvent come from that program's alias in the chain
         YAML. Not used for analyze.

.. note::

   -  ``-p`` uses the chain project name without the ``.yaml`` extension.

   -  For the CLI pipeline, ``-f`` and at least one ``-s/--steps`` are required. The ``custom`` subcommand runs the same
      pipeline.

   -  For workflow submit, ``-p``, ``-f``, and ``--program`` are required. Analyze does not need ``-p``, ``-f``, or
      ``--program``.

   -  ``-l`` and ``-a`` are mutually exclusive.

   -  Workflow option tables (proton index, Fukui mode, redox reference files, reaction ``--product``, …) match the
      corresponding program CLI. See :doc:`pka-calculations`, :ref:`fukui-jobs`, :doc:`redox-calculations`, and
      :doc:`reaction`.

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

Workflow submit from the chain alias:

.. code:: bash

   chemsmart sub chain -p combined -f ligand.xyz -c 0 -m 1 \
     pka --program gaussian -pi 10 submit
   chemsmart sub chain -p combined -f ligand.xyz -c 0 -m 1 \
     fukui --program orca
   chemsmart sub chain -p combined -f ox.xyz -c 1 -m 1 \
     redox --program gaussian --ref-ox ref_ox.xyz --ref-red ref_red.xyz
   chemsmart sub chain -p combined -f ts.xyz -c 0 -m 1 \
     reaction --program gaussian

Analyze without a chain project:

.. code:: bash

   chemsmart run chain pka analyze -ha acid_HA_opt.log -hr ref_HRef_opt.log -rp 6.75
   chemsmart run chain fukui analyze -n mol_n.log
   chemsmart run chain redox analyze --ox-gas ox_opt.log --red-gas red_opt.log \
     --ref-ox-gas RefOx_opt.log --ref-red-gas RefRed_opt.log \
     --ox-solv ox_sp.log --red-solv red_sp.log \
     --ref-ox-solv RefOx_sp.log --ref-red-solv RefRed_sp.log

Ordinary Gaussian opt or pKa uses that program's own project YAML:

.. code:: bash

   chemsmart sub gaussian -p gaussian_project2 -f ligand.xyz -c 0 -m 1 -o maxstep=8,maxsize=12 opt
   chemsmart sub gaussian -p gaussian_project2 -f ligand.xyz -c 0 -m 1 pka submit

************
 Next Steps
************

-  :doc:`chain-jobs` — chain YAML aliases, ``-s/--steps`` pipeline, workflow subcommands, geometry handoff, and HPC
   submission
-  :doc:`configuration-project-settings` — chain project file location and alias rules
-  :doc:`pka-calculations` — pKa submit and analyze
-  :doc:`redox-calculations` — exchange redox submit and analyze
-  :doc:`reaction` — reaction submit (no analyze)
