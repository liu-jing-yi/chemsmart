.. _reaction-calculations:

###################
 Reaction Workflow
###################

CHEMSMART provides a program-first reaction chain nested under Gaussian and ORCA, the same way pKa and Fukui already
are:

.. code:: bash

   chemsmart sub gaussian [GAUSSIAN_OPTIONS] reaction [REACTION_OPTIONS]
   chemsmart sub orca     [ORCA_OPTIONS]     reaction [REACTION_OPTIONS]

The chain locates a transition state from a guess **or** from reactant + product (Gaussian QST2/QST3 or ORCA NEB-TS),
then runs TS opt+freq, optional reactant/product opt+freq, and dual-level solvent single-points.

.. note::

   Output analysis is **not** part of this command. There is no ``chemsmart run reaction analyze`` yet. After jobs
   finish, use :doc:`thermochemistry-analysis` on the individual opt/SP outputs. A dedicated reaction analysis command
   is planned as a follow-up.

.. toctree::
   :maxdepth: 2
   :caption: Reaction Submission

   gaussian-reaction-calculations
   orca-reaction-calculations

.. contents:: Table of Contents
   :local:
   :depth: 2

************************
 Execution Architecture
************************

**Job submission**

-  ``chemsmart run/sub gaussian ... reaction [submit|batch]`` — Gaussian QST (when locating a TS) then R/TS/P opt and
   SP.
-  ``chemsmart run/sub orca ... reaction [submit|batch]`` — ORCA NEB-TS (when locating a TS) then R/TS/P opt and SP.
-  Use ``chemsmart run`` for local preparation and execution; use ``chemsmart sub`` on HPC clusters to generate
   scheduler scripts.
-  When ``reaction`` is invoked without an explicit subcommand, a submission table triggers ``batch``; otherwise
   ``submit`` runs.

``gaussian ts`` / ``orca ts`` remain single-structure TS searches. ``orca neb`` remains the standalone NEB command.
``gaussian reaction`` / ``orca reaction`` are the R/TS/P workflow.

******************
 Case 1 vs Case 2
******************

The chain always ends in the same TS characterization (opt=ts / OptTS + freq, then optional R/P opt + SP). Case 2 only
**creates** that TS guess.

**Case 1 — TS guess provided.** Skip Guess. Parent ``-f`` is the TS structure. Optional ``--reactant`` (without
``--product``) adds reactant minima to optimize; that is not a path search.

**Case 2 — reactant + product, optional TS guess.** Run Guess, then case 1 with the located TS:

-  **Gaussian:** QST2 if only R+P; QST3 if a TS guess is also given. Atom order must match across coordinate blocks.
-  **ORCA:** project NEB-TS. Reactant is the NEB start, product is the ending geometry, optional TS guess is the
   intermediate.

Failed QST/NEB does not continue to TS opt.

CLI dispatch
============

-  ``-f`` only → case 1, ``-f`` is TS
-  ``-f`` + ``--reactant``, no ``--product`` → case 1, ``-f`` is TS
-  ``-f`` + ``--product``, no ``--reactant`` → case 2, ``-f`` is reactant (same as ``orca neb -f … -e``)
-  ``--reactant`` + ``--product`` → case 2; ``-f`` is the QST3/NEB intermediate
-  ``--ts-guess`` with ``--product`` → QST3 / NEB intermediate when ``-f`` is the reactant

.. code:: bash

   # Case 1
   chemsmart sub gaussian -p proj -f ts_guess.xyz -c 0 -m 1 reaction

   # Case 2 Gaussian QST2 / QST3
   chemsmart sub gaussian -p proj -f reactant.xyz -c 0 -m 1 reaction --product product.xyz
   chemsmart sub gaussian -p proj -f reactant.xyz reaction --product product.xyz --ts-guess ts.xyz

   # Case 2 ORCA NEB then OptTS
   chemsmart sub orca -p proj -f reactant.xyz -c 0 -m 1 reaction --product product.xyz
   chemsmart sub orca -p proj -f reactant.xyz reaction --product product.xyz --ts-guess ts.xyz

**********************************
 What Happens After the TS Exists
**********************************

Roles:

-  **ts** (required): TS job with ``freq=True``, project ``ts_settings()`` (gas-phase)
-  **reactant** / **product** (optional): opt+freq, project ``opt_settings()``
-  **SP**: solvent single-points on optimized geometries, project ``sp_settings()`` (solv) — a different level of theory
   than opt when ``gas`` / ``solv`` differ

Phases:

#. **Guess** — QST or NEB; skipped in case 1
#. **Opt** — R, then TS, then P
#. **SP** — solvent single-points

IRC is not part of this workflow.

Child labels
============

For a job labelled ``sn2``:

.. code:: text

   sn2_qst.com / sn2_neb.inp   # Guess (case 2)
   sn2_TS_opt                 # TS opt+freq
   sn2_R_opt / sn2_R1_opt     # reactant opt+freq (numbered when more than one)
   sn2_P_opt / sn2_P1_opt     # product opt+freq
   sn2_TS_sp / sn2_R_sp / …   # matching solvent single-points

Repeatable ``--reactant`` / ``--product`` cover extra fragments (two reactants, etc.).

*************
 Batch Table
*************

Pass a ``.csv`` or whitespace-delimited ``.txt`` file via ``-f`` and invoke ``reaction batch`` (or omit the subcommand —
batch is selected automatically when ``-f`` points to a submission table).

Required columns:

-  ``reaction_id`` (aliases: ``reaction``, ``id``, ``name``)
-  ``filepath`` (aliases: ``file_path``, ``path``)
-  ``role`` (aliases: ``type``) — ``ts``, ``reactant``, or ``product``
-  ``charge`` (alias: ``q``)
-  ``multiplicity`` (aliases: ``mult``, ``m``)

Role aliases: ``transition_state`` / ``ts_guess`` → ``ts``; ``r`` / ``reactants`` → ``reactant``; ``p`` / ``products`` →
``product``.

Example ``reactions.csv``:

.. code:: text

   reaction_id,filepath,role,charge,multiplicity
   sn2,ts.xyz,ts,0,1
   sn2,nu.xyz,reactant,-1,1
   sn2,rx.xyz,reactant,0,1
   sn2,prod.xyz,product,-1,1

Rows with the same ``reaction_id`` become one job. A group with both reactant and product roles is case 2. A ``ts`` row
without a reactant is case 1 (``-f`` is the TS).

.. note::

   When ``-f`` is a submission table, the parent ``gaussian`` / ``orca`` command does not require ``-c`` / ``-m``;
   charge and multiplicity are read from each table row.

On HPC clusters, ``chemsmart sub ... reaction batch`` writes one scheduler script per ``reaction_id``. Each run wrapper
is rewritten to a single-reaction ``reaction submit`` command (table path replaced by the parent structure, ``batch``
replaced by ``submit``, ``--reactant`` / ``--product`` injected). See :doc:`cli-overview` and
:doc:`configuration-server-settings`.

**********
 Settings
**********

Child jobs copy project ``opt_settings()`` / ``ts_settings()`` / ``sp_settings()`` (and ORCA ``neb_settings()`` for
Guess). Charge and multiplicity are set per structure. There is no new YAML job type: use an existing Gaussian or ORCA
project.

ORCA Hessian/ScanTS flags are not re-exposed on ``reaction``; TS children use project ``ts_settings()``. Guess-phase
ORCA uses project ``neb_settings()`` (default ``NEB-TS``).

No ``xtb reaction`` in this version.

******************
 Analysis (later)
******************

``chemsmart run reaction analyze`` is not implemented. Until that lands, compute thermochemistry on the child outputs
with :doc:`thermochemistry-analysis`.

**********
 See Also
**********

-  :ref:`gaussian-reaction-calculations`
-  :ref:`orca-reaction-calculations`
-  :doc:`gaussian-transition-state`
-  :doc:`orca-transition-state`
-  :doc:`thermochemistry-analysis`
