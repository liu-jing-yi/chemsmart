.. _redox-calculations:

####################
 Redox Calculations
####################

CHEMSMART provides redox workflows in two separate stages:

#. **Job submission** — generate and run Gaussian or ORCA calculations for the oxidized and reduced target and a
   reference couple. See :ref:`gaussian-redox-calculations` and :ref:`orca-redox-calculations`.

#. **Output analysis** — compute the exchange redox potential from completed output files using the backend-independent
   command ``chemsmart run redox analyze`` (also ``chemsmart run chain redox analyze``). Analysis is
   **program-agnostic**: the same workflow reads Gaussian ``.log`` and ORCA ``.out`` files.

The scheme is **exchange-only**. There is no isolated-electron :math:`G(\mathrm{e}^{-})` term and no bundled Fc/Fc+
geometries.

.. toctree::
   :maxdepth: 2
   :caption: Redox Calculations

   gaussian-redox-calculations
   orca-redox-calculations

.. contents:: Table of Contents
   :local:
   :depth: 2

************************
 Execution Architecture
************************

**Job submission**

-  ``chemsmart run/sub gaussian ... redox`` — prepare and run Gaussian redox calculations.
-  ``chemsmart run/sub orca ... redox`` — prepare and run ORCA redox calculations.
-  ``chemsmart run/sub chain -p combined ... redox --program {gaussian,orca}`` — same jobs, with theory and solvent from
   the chain YAML alias for ``--program``. See :ref:`chain-workflow-subcommands`.
-  Use ``chemsmart run`` for local preparation and execution; use ``chemsmart sub`` on HPC clusters to generate
   scheduler scripts.
-  The oxidized target comes from parent ``-f``. The reduced target uses the same geometry (or ``--red``) with charge
   ``ox − n``. Reference geometries come from the registry and/or ``--ref-ox`` / ``--ref-red``.

**Output analysis**

-  ``chemsmart run redox analyze`` — single-system analysis from eight output files. The reference couple is inferred
   from Ref_ox/Ref_red formulas; ``-r`` overrides.
-  ``chemsmart run chain redox analyze`` — the same analysis under ``chain`` (no ``-p``, ``-f``, or ``--program``).
-  Analysis never invokes ``gaussian`` or ``orca`` job submission — only reads completed output files.

********
 Theory
********

The exchange reaction is:

.. math::

   \mathrm{Ox} + \mathrm{Ref_{red}} \rightarrow \mathrm{Red} + \mathrm{Ref_{ox}}

The exchange free-energy change is:

.. math::

   \Delta G_{\mathrm{exchange}} = G(\mathrm{Red}) + G(\mathrm{Ref_{ox}}) - G(\mathrm{Ox}) - G(\mathrm{Ref_{red}})

The target reduction potential on the reference scale is:

.. math::

   E_{\mathrm{target}} = E_{\mathrm{ref}} - \frac{\Delta G_{\mathrm{exchange}}}{n F}

:math:`n` must match the target couple and the reference couple. :math:`F` is Faraday's constant in C/mol. :math:`\Delta
G_{\mathrm{exchange}}` is converted from Hartree to J/mol before dividing by :math:`n F`. The summary reports
:math:`\Delta G` in au, kcal/mol, eV, and J/mol.

*********************
 Dual-Level Approach
*********************

Solution free energies follow the same dual-level construction as pKa:

#. **Thermal corrections** (:math:`G_{\mathrm{corr}}`) from gas-phase frequency calculations using quasi-harmonic Gibbs
   free energy:

   .. math::

      G_{\mathrm{corr}} = G_{\mathrm{qh}}(T) - E_{\mathrm{gas}}

#. **Solvent energies** (:math:`E_{\mathrm{solv}}`) from high-level single-point calculations in implicit solvent.

#. **Total free energy in solution**:

   .. math::

      G_{\mathrm{soln}} = E_{\mathrm{solv}} + G_{\mathrm{corr}}

**********************************
 Job Submission (Gaussian / ORCA)
**********************************

Job submission is backend-specific. Use the dedicated pages for full examples:

-  :ref:`gaussian-redox-calculations`
-  :ref:`orca-redox-calculations`

**Commands**

The built-in ``fc_fc+`` couple has :math:`E_{\mathrm{ref}} = 0.0` V and :math:`n = 1` on the Fc/Fc+ scale. It does not
bundle geometries, so ``--ref-ox`` and ``--ref-red`` are required unless another registered couple supplies them.

.. code:: bash

   chemsmart run gaussian -p my_project -f ox.xyz -c 1 -m 2 redox \
       --ref-ox ref_ox.xyz --ref-red ref_red.xyz

   chemsmart run orca -p my_project -f ox.xyz -c 1 -m 2 redox \
       --ref-ox ref_ox.xyz --ref-red ref_red.xyz

   chemsmart run chain -p combined -f ox.xyz -c 1 -m 2 \
       redox --program gaussian --ref-ox ref_ox.xyz --ref-red ref_red.xyz

Phases: Opt (Ox, Red) → Ref Opt → SP → Ref SP. Child labels for a job labelled ``mol_redox``:

.. code:: text

   mol_redox_ox_opt / mol_redox_red_opt
   mol_redox_RefOx_opt / mol_redox_RefRed_opt
   mol_redox_ox_sp / mol_redox_red_sp
   mol_redox_RefOx_sp / mol_redox_RefRed_sp

***************************
 Reference Couple Registry
***************************

Core jobs consume a ``RedoxReference`` from the public registry. Register additional couples without editing the redox
job class:

.. code:: python

   from chemsmart.cli.redox import (
       RedoxReference,
       get_redox_reference,
       list_redox_references,
       register_redox_reference,
   )

   register_redox_reference(
       RedoxReference(
           name="custom_she",
           E_ref_V=0.40,
           n_electrons=1,
           scale="SHE",
           couple_label="Custom/Custom+",
           ox_file="ref_ox.xyz",
           red_file="ref_red.xyz",
           ox_charge=1,
           ox_multiplicity=2,
           red_charge=0,
           red_multiplicity=1,
       )
   )

``-r/--reference`` selects a registry name. Submit defaults to ``fc_fc+``. Analyze infers the couple from the Hill
formulas of the Ref_ox and Ref_red outputs; ``-r`` overrides. The built-in ``fc_fc+`` couple matches
:math:`\mathrm{C_{10}H_{10}Fe}`. If inference is ambiguous or the formula is unregistered, pass ``-r``.
``-n/--n-electrons`` defaults to the couple and must match it when given.

*******************************************
 Output Analysis (``chemsmart run redox``)
*******************************************

All post-processing lives under ``chemsmart run redox analyze`` or ``chemsmart run chain redox analyze``. No Gaussian or
ORCA backend is invoked during analysis.

Provide all eight outputs. The reference couple is inferred from Ref_ox/Ref_red formulas; pass ``-r`` to override.

.. code:: bash

   chemsmart run redox analyze \
       --ox-gas mol_redox_ox_opt.log \
       --red-gas mol_redox_red_opt.log \
       --ref-ox-gas mol_redox_RefOx_opt.log \
       --ref-red-gas mol_redox_RefRed_opt.log \
       --ox-solv mol_redox_ox_sp.log \
       --red-solv mol_redox_red_sp.log \
       --ref-ox-solv mol_redox_RefOx_sp.log \
       --ref-red-solv mol_redox_RefRed_sp.log \
       -T 333.15 -c 1.0 -csg 100 -ch 100 \
       -o redox.dat

Pass ``-n`` / ``-r`` on the ``redox`` group **before** ``analyze`` when you need to override inference:

.. code:: bash

   chemsmart run redox -n 1 -r fc_fc+ analyze --ox-gas ...
   chemsmart run chain redox -n 1 -r fc_fc+ analyze --ox-gas ...

*************
 CLI Options
*************

Submit Options
==============

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Option
      -  Description
   -  -  ``-r, --reference``
      -  Registry name of the reference couple (default ``fc_fc+``).
   -  -  ``-n, --n-electrons``
      -  Electrons transferred. Defaults to the reference couple; must match it when given.
   -  -  ``-rd, --red``
      -  Reduced target geometry. Defaults to the oxidized structure from parent ``-f`` with charge ``ox − n``.
   -  -  ``-rdc, --red-charge`` / ``-rdm, --red-multiplicity``
      -  Charge and multiplicity of the reduced target.
   -  -  ``-ro, --ref-ox`` / ``-rr, --ref-red``
      -  Oxidized and reduced reference geometries. Required unless the registered couple provides them.
   -  -  ``-roc, --ref-ox-charge`` / ``-rom, --ref-ox-multiplicity``
      -  Charge and multiplicity of the oxidized reference.
   -  -  ``-rrc, --ref-red-charge`` / ``-rrm, --ref-red-multiplicity``
      -  Charge and multiplicity of the reduced reference.
   -  -  ``-T, --temperature`` / ``-c, --concentration`` / ``-csg`` / ``-ch``
      -  Thermochemistry options stored on the job (same defaults as pKa analysis).

Analyze Options
===============

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Option
      -  Description

   -  -  ``-r, --reference``
      -  Registry name of the reference couple. Default: inferred from Ref_ox/Ref_red Hill formulas. Overrides
         inference. Pass it on the ``redox`` group before ``analyze``.

   -  -  ``-n, --n-electrons``
      -  Electrons transferred. Defaults to the reference couple; must match it when given. Pass it on the ``redox``
         group before ``analyze``.

   -  -  ``-oxg, --ox-gas`` / ``-rdg, --red-gas``
      -  Gas-phase opt+freq outputs for the target couple.

   -  -  ``-rog, --ref-ox-gas`` / ``-rrg, --ref-red-gas``
      -  Gas-phase opt+freq outputs for the reference couple.

   -  -  ``-oxs, --ox-solv`` / ``-rds, --red-solv``
      -  Solution-phase SP outputs for the target couple.

   -  -  ``-ros, --ref-ox-solv`` / ``-rrs, --ref-red-solv``
      -  Solution-phase SP outputs for the reference couple.

   -  -  ``-T, --temperature`` / ``-c, --concentration`` / ``-csg`` / ``-ch``
      -  Thermochemistry options for quasi-harmonic G (same defaults as pKa analysis).

   -  -  ``-o, --output``
      -  Write the formatted summary to this file instead of printing it.

**********
 See Also
**********

-  :ref:`gaussian-redox-calculations`
-  :ref:`orca-redox-calculations`
-  :ref:`chain-workflow-subcommands`
-  :ref:`pka-calculations`
-  :doc:`thermochemistry-analysis`
