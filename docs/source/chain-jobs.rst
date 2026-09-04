.. _chain-jobs:

#################
 Chain Workflows
#################

Chain projects tie together CREST, xTB, Gaussian, and ORCA calculations. A chain YAML file holds **project aliases**
only; pipeline order comes from repeatable ``-s/--steps`` on the CLI. Combined YAML is used only by ``chain`` jobs.

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

Top-level keys ``crest``, ``xtb``, ``gaussian``, and ``orca`` are aliases to existing per-program project names (the
YAML stems under ``~/.chemsmart/crest/``, ``xtb/``, ``gaussian/``, and ``orca/``). Any other top-level key is rejected.
A legacy ``steps`` key in the file is also rejected — specify the pipeline with ``-s/--steps`` instead.

Each pipeline step's program must have a corresponding alias in the file. The target per-program YAML must exist;
otherwise CHEMSMART raises the same ``from_project`` error as for a missing Gaussian or ORCA project. ``chemsmart sub
gaussian -p test`` still loads ``~/.chemsmart/gaussian/test.yaml`` only — it does not read chain YAML.

Packaged templates are copied to ``~/.chemsmart/chain/`` by ``chemsmart config`` and ``chemsmart update projects``. See
:doc:`configuration-project-settings` for the full schema.

********************
 CLI Pipeline Steps
********************

Repeat ``-s`` on the chain command for each phase. A step is ``PROGRAM:JOB`` or a quoted ``PROGRAM: [OPTIONS] JOB``
string that matches the program CLI after the program name:

.. code:: bash

   chemsmart run chain -p combined -f mol.xyz -c 0 -m 1 \
     -s "gaussian: -o maxstep=8,maxsize=12 opt" \
     -s gaussian:sp

The following ``(program, job)`` pairs are supported:

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

   -  -  ``orca``
      -  ``opt``, ``ts``, ``sp``
      -  ORCA optimization, transition state, or single point

Jobs that need large option surfaces or extra input structures — **pKa**, **Fukui**, **reaction**, QM/MM, IRC, scan, and
similar — are **not** pipeline steps. Run them on the program CLI with that program's own project YAML, e.g. ``chemsmart
sub gaussian -p gaussian_project2 pka submit``.

Geometry handoff between phases:

-  Optimization, TS, and SP steps use the previous job's optimized structure.
-  CREST passes the best conformer (``crest_best.xyz``).
-  If the previous phase is incomplete, the chain stops so HPC resubmission with ``--skip-completed`` can continue
   later.

Child job labels follow ``{chain_label}_{index:02d}_{program}_{job}`` (e.g. ``mol_00_gaussian_opt``,
``mol_01_gaussian_opt`` for two Gaussian optimization steps). Charge and multiplicity from ``-c`` / ``-m`` on the chain
command are merged into each child.

********************
 Running a Pipeline
********************

.. code:: bash

   chemsmart sub chain -p combined -f molecule.xyz -c 0 -m 1 \
     -s crest:conformers -s xtb:opt -s gaussian:opt -s orca:sp

The whole pipeline is **one** job submission: a single ``ChainJob`` orchestrates all phases. Use ``chemsmart sub`` on
HPC clusters the same way as for Gaussian pKa or other multi-phase workflows — one scheduler script runs the full chain.

Local preparation and dry runs:

.. code:: bash

   chemsmart run chain -p combined -f molecule.xyz -c 0 -m 1 \
     -s crest:conformers -s xtb:opt -s gaussian:opt -s orca:sp
   chemsmart run --fake chain -p combined -f molecule.xyz -c 0 -m 1 \
     -s crest:conformers -s xtb:opt -s gaussian:opt -s orca:sp

``--fake`` uses the chain fake runner for the parent; each phase child receives the appropriate program fake runner when
its job type differs from the parent.

*************
 Program CLI
*************

``chemsmart sub chain --help`` lists ``custom``. pKa, Fukui, reaction, and other single-program jobs stay on the program
CLI and use that program's own ``-p`` project YAML:

.. code:: bash

   chemsmart sub chain -p combined -f mol.xyz -s gaussian:opt custom
   chemsmart sub gaussian -p gaussian_project2 -f mol.xyz -c 0 -m 1 pka submit
   chemsmart sub gaussian -p gaussian_project2 -f mol.xyz -c 0 -m 1 -o maxstep=8,maxsize=12 opt

``chemsmart run gaussian -p test ...`` still loads ``~/.chemsmart/gaussian/test.yaml`` only.

***************
 CLI Reference
***************

See :doc:`chain-cli-options` for chain-specific options. Program subcommands accept the same options as documented in
their respective CLI pages.
