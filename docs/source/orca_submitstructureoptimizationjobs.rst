Welcome to the tutorials! We're thrilled to have you here. Please go through the code examples, and don't hesitate to
contact our team if you have questions or feedback.

####################################
 Submit Structure Optimization Jobs
####################################

ChemSmart provides comprehensive tools for structure optimization calculations using ORCA. This section covers geometry
optimization workflows including constrained optimizations.

************************
 Structure Optimization
************************

Geometry optimization is used to find the minimum energy structure of a molecule by adjusting atomic positions until
forces are minimized.

.. code:: console

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] opt [SUBCMD_OPTIONS]

Opt-Specific OPTIONS
====================

.. list-table:: Structure Optimization Options
   :header-rows: 1
   :widths: 30 15 55

   -  -  Option
      -  Type
      -  Description

   -  -  ``-f, --freeze-atoms``
      -  string
      -  Indices of atoms to freeze for constrained optimization. 1-indexed (default=None)

   -  -  ``-i, --invert-constraints/--no-invert-constraints``
      -  bool
      -  Invert the constraints for frozen atoms in optimization (default=False)

Basic Usage
===========

**Basic geometry optimization**:
To submit a regular ORCA optimization job, one can do:

   .. code:: console

      chemsmart sub orca -p project_name -f input.xyz opt

**Optimization with inverted constraints**:
To submit a constrained ORCA optimization with specified atoms frozen, one can do:

   .. code:: console

      chemsmart sub orca -p inverted_opt -f molecule.xyz opt -f 1,2,3 -i

Please note that the ``opt`` subcommand option ``-f, --freeze-atoms`` is different
from the main ``-f, --file`` option used to specify the input file.

Examples
========

**optimize structure directly from an ORCA output file with different charge and multiplicity**

-  After an ORCA optimization job is done, one can use the output file from the previous job to run a
   new optimization with different charge and multiplicity, e.g.:

   .. code:: console

      chemsmart sub -s SLURM ORCA -p project1 -f k_atom_opt.out -c 1 -m 1 -l k_cation_opt opt

   ``-c`` and ``-m`` will override the charge and multiplicity in the original ORCA out file (from K_atom to
   K_cation).


***************************
 Single Point Calculations
***************************

Single point calculations compute the energy and properties of a molecule at a fixed geometry without optimization.

.. code:: console

   chemsmart sub [OPTIONS] orca [ORCA_OPTIONS] sp [SUBCMD_OPTIONS]

Basic Usage
===========

**Basic single point calculation**:

   .. code:: console

      chemsmart sub orca -p project_name -f input.xyz sp
