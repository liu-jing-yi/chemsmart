.. _chain-jobs:

#################
 Chain Workflows
#################

Chain projects tie together CREST, xTB, Gaussian, and ORCA calculations in a single YAML file. Use them in two ways:

#. **Pipeline** — ``chemsmart run/sub chain -p NAME -f mol.xyz`` runs the ``steps`` list in order, passing optimized
   geometries from each phase to the next.
#. **Nested slice** — ``chemsmart run/sub chain -p NAME gaussian pka ...`` (or ``orca``, ``xtb``, ``crest``) runs one
   program's CLI with project aliases from the chain file, without loading other programs.

.. contents:: Table of Contents
   :local:
   :depth: 2

********************
 Chain Project YAML
********************

Chain project files live in ``~/.chemsmart/chain/``. The ``-p`` / ``--project`` option takes the file stem (without
``.yaml``), the same convention as ``gaussian -p test`` → ``~/.chemsmart/gaussian/test.yaml``.

Example (``~/.chemsmart/chain/combined.yaml``):

.. code:: yaml

   crest: crest_project1
   xtb: xtb_project1
   gaussian: gaussian_project2
   orca: orca_project3

   steps:
     - program: crest
       job: conformers
     - program: xtb
       job: opt
     - program: gaussian
       job: opt
     - program: orca
       job: sp

Top-level keys ``crest``, ``xtb``, ``gaussian``, and ``orca`` are aliases to existing per-program project names (the
YAML stems under ``~/.chemsmart/crest/``, ``xtb/``, ``gaussian/``, and ``orca/``). Unknown top-level keys other than
``steps`` are rejected.

``steps`` is optional. An aliases-only file is valid for nested slices. Running ``chain -p NAME`` with no nested command
and no ``steps`` is a usage error.

Each step's ``program`` must have a corresponding alias in the file. The target per-program YAML must exist; otherwise
CHEMSMART raises the same ``from_project`` error as for a missing Gaussian or ORCA project.

Packaged templates are copied to ``~/.chemsmart/chain/`` by ``chemsmart config`` and ``chemsmart update projects``. See
:doc:`configuration-project-settings` for the full schema.

*********************
 YAML Pipeline Steps
*********************

The following ``(program, job)`` pairs are supported in ``steps``:

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   -  -  Program
      -  Job
      -  Description

   -  -  ``crest``
      -  ``conformers``
      -  CREST conformer search; best conformer is passed to the next step

   -  -  ``xtb``
      -  ``opt``
      -  xTB geometry optimization (with frequency by default)

   -  -  ``xtb``
      -  ``sp``
      -  xTB single-point energy

   -  -  ``xtb``
      -  ``hess``
      -  xTB Hessian / frequency calculation

   -  -  ``gaussian``
      -  ``opt``, ``ts``, ``sp``
      -  Gaussian optimization, transition state, or single point

   -  -  ``gaussian``
      -  ``modred``, ``irc``, ``scan``, ``nci``, ``wbi``, ``td``, ``resp``, ``fukui``, ``qrc``
      -  Other Gaussian job types with standard project settings

   -  -  ``orca``
      -  ``opt``, ``ts``, ``sp``, ``modred``, ``irc``, ``scan``, ``fukui``, ``qrc``
      -  ORCA equivalents of the Gaussian rows above

Jobs that need large option surfaces or extra input structures — **pKa**, **reaction**, **QM/MM**, and similar — are
**not** YAML steps. Run them as nested commands, e.g. ``chain -p combined gaussian pka submit``.

Geometry handoff between phases:

-  Optimization, TS, and SP steps use the previous job's optimized structure.
-  CREST passes the best conformer (``crest_best.xyz``).
-  If the previous phase is incomplete, the chain stops so HPC resubmission with ``--skip-completed`` can continue
   later.

Child job labels follow ``{chain_label}_{program}_{job}`` (e.g. ``mol_gaussian_opt``). Charge and multiplicity from
``-c`` / ``-m`` on the chain command are merged into each child.

********************
 Running a Pipeline
********************

.. code:: bash

   chemsmart sub chain -p combined -f molecule.xyz -c 0 -m 1

The whole pipeline is **one** job submission: a single ``ChainJob`` orchestrates all phases. Use ``chemsmart sub`` on
HPC clusters the same way as for Gaussian pKa or other multi-phase workflows — one scheduler script runs the full chain.

Local preparation and dry runs:

.. code:: bash

   chemsmart run chain -p combined -f molecule.xyz -c 0 -m 1
   chemsmart run --fake chain -p combined -f molecule.xyz -c 0 -m 1

``--fake`` uses the chain fake runner for the parent; each phase child receives the appropriate program fake runner when
its job type differs from the parent.

***********************
 Nested Program Slices
***********************

Nested commands reuse the full Gaussian, ORCA, xTB, and CREST CLI groups:

.. code:: bash

   chemsmart sub chain -p combined gaussian -f mol.xyz -c 0 -m 1 pka submit
   chemsmart sub chain -p combined orca -f mol.xyz opt
   chemsmart sub chain -p combined xtb -f mol.xyz opt
   chemsmart sub chain -p combined crest -f mol.xyz conformers

Project resolution for nested commands:

-  If the nested program's own ``-p`` is set, that project is used (explicit override).
-  Otherwise, the alias from the chain YAML is used (e.g. ``gaussian: gaussian_project2``).
-  If the chain file has no section for that program, CHEMSMART raises a usage error.

``chemsmart run gaussian -p test ...`` is unchanged — it still loads ``~/.chemsmart/gaussian/test.yaml`` directly
without a chain file.

***************
 CLI Reference
***************

See :doc:`chain-cli-options` for chain-specific options. Nested ``gaussian``, ``orca``, ``xtb``, and ``crest``
subcommands accept the same options as documented in their respective CLI pages.
