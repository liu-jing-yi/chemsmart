.. _gaussian-reaction-calculations:

################################
 Gaussian Reaction Calculations
################################

This page covers **Gaussian reaction job submission** — QST2/QST3 when locating a TS from reactant and product, then
TS/R/P opt+freq and solvent single-points.

.. note::

   Output-file analysis is not part of ``gaussian reaction``. After calculations finish, use ``chemsmart run
   thermochemistry`` on the child outputs. See :ref:`reaction-calculations`. There is no ``chemsmart run reaction
   analyze`` and no ``chain reaction analyze``.

.. contents:: Table of Contents
   :local:
   :depth: 2

*************
 Quick Start
*************

**Case 1 — TS guess**

``-f`` is the TS structure. Guess (QST) is skipped. Optional ``--reactant`` without ``--product`` adds reactant optima
to optimize.

.. code:: bash

   chemsmart sub gaussian -p my_project -f ts_guess.xyz -c 0 -m 1 reaction

Where:

-  ``-p my_project``: Project settings (gas opt/TS theory and solv SP theory)
-  ``-f ts_guess.xyz``: TS guess (XYZ, LOG, COM, …)
-  ``-c 0 -m 1``: Charge and multiplicity of the TS

This runs TS opt+freq (``ts_settings()``) and a solvent single-point (``sp_settings()``).

Chain submit uses the Gaussian alias from the chain YAML:

.. code:: bash

   chemsmart sub chain -p combined -f ts_guess.xyz -c 0 -m 1 \
       reaction --program gaussian

**Case 2 — reactant + product (QST2)**

``-f`` is the reactant. ``--product`` selects path search. Gaussian builds a QST2 ``.com`` (two geometry blocks) and
runs it as a COM job, then characterizes the located TS.

.. code:: bash

   chemsmart sub gaussian -p my_project -f reactant.xyz -c 0 -m 1 reaction \
       --product product.xyz

**Case 2 — QST3**

Pass a TS guess as well. Atom order must match across reactant, product, and TS guess.

.. code:: bash

   chemsmart sub gaussian -p my_project -f reactant.xyz -c 0 -m 1 reaction \
       --product product.xyz --ts-guess ts.xyz

   # Equivalent: -f is the QST3 guess when both --reactant and --product are set
   chemsmart sub gaussian -p my_project -f ts.xyz -c 0 -m 1 reaction \
       --reactant reactant.xyz --product product.xyz

``gaussian ts`` remains the single-structure TS search. See :doc:`gaussian-transition-state`.

************************
 Job Output File Naming
************************

Sub-job labels determine output filenames. For a job with label ``sn2`` (the default when submitting ``sn2.xyz``):

.. code:: text

   sn2_qst.com / sn2_qst.log     # QST2 or QST3 (case 2 only)
   sn2_TS_opt.log                # TS gas-phase opt+freq
   sn2_R_opt.log                 # reactant opt+freq (if given)
   sn2_R1_opt.log / sn2_R2_opt   # extra reactant fragments
   sn2_P_opt.log                 # product opt+freq (if given)
   sn2_TS_sp.log                 # TS solvent single-point
   sn2_R_sp.log / sn2_P_sp.log   # matching SP children

************************************
 Batch Processing with Input Tables
************************************

Pass a ``.csv`` or whitespace-delimited ``.txt`` file via ``-f`` and invoke ``reaction batch`` (or omit the subcommand —
batch is selected automatically when ``-f`` points to a submission table).

.. code:: bash

   chemsmart sub gaussian -p my_project -f reactions.csv reaction batch

.. note::

   When ``-f`` is a submission table, the parent ``gaussian`` command does not require ``-c`` / ``-m``; charge and
   multiplicity are read from each table row.

On HPC clusters, use ``chemsmart sub`` instead of ``chemsmart run``; each ``reaction_id`` receives its own scheduler
script with a reconstructed ``reaction submit`` command. See :ref:`reaction-calculations`.

Table format is documented in :ref:`reaction-calculations`.

************
 Parameters
************

Reaction Options
================

.. list-table::
   :header-rows: 1
   :widths: 15 20 65

   -  -  Short
      -  Long
      -  Description

   -  -  ``-r``
      -  ``--reactant``
      -  Reactant geometry file. Repeatable for extra fragments. With ``--product``, parent ``-f`` is the TS guess.

   -  -
      -  ``--product``
      -  Product geometry file. Repeatable. Presence selects path search (case 2) when no ``--reactant`` is given.

   -  -
      -  ``--ts-guess``
      -  QST3 intermediate when ``-f`` is the reactant. Requires ``--product``.

   -  -  ``-S`` / ``-R``
      -  ``--skip-completed`` / ``--no-skip-completed``
      -  Skip completed child jobs (default) or rerun them.

QST notes
=========

-  QST2 uses two coordinate blocks (reactant, product). QST3 adds a TS-guess block.
-  Structures must have the same number of atoms and the same atom order.
-  The Guess job reuses project TS theory with ``opt=qst2`` or ``opt=qst3``; it does not change the shared Gaussian
   writer or route builder.

**********
 Examples
**********

Example 1: TS Guess Only
========================

.. code:: bash

   chemsmart sub gaussian -p b3lyp_project -f ts_guess.xyz -c 0 -m 1 reaction

Example 2: QST2 from Reactant and Product
=========================================

.. code:: bash

   chemsmart sub gaussian -p b3lyp_project -f reactant.xyz -c 0 -m 1 reaction \
       --product product.xyz

Example 3: QST3
===============

.. code:: bash

   chemsmart sub gaussian -p b3lyp_project -f reactant.xyz -c 0 -m 1 reaction \
       --product product.xyz --ts-guess ts.xyz

Example 4: Batch Submission from CSV
====================================

.. code:: bash

   chemsmart sub gaussian -p b3lyp_project -f reactions.csv reaction batch

Example 5: Thermochemistry on Completed Outputs
===============================================

Until ``chemsmart run reaction analyze`` exists:

.. code:: bash

   chemsmart run thermochemistry -f sn2_TS_opt.log
   chemsmart run thermochemistry -f sn2_R_opt.log
   chemsmart run thermochemistry -f sn2_P_opt.log

**********
 See Also
**********

-  :ref:`reaction-calculations`
-  :ref:`orca-reaction-calculations`
-  :doc:`gaussian-transition-state`
-  :doc:`thermochemistry-analysis`
