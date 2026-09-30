.. _pka-calculations:

##################
 pKa Calculations
##################

CHEMSMART provides pKa workflows in two separate stages:

#. **Job submission** — generate and run Gaussian or ORCA calculations for HA, A⁻, and (optionally) a reference acid.
   With ``--pkb``, the input is the free base B; CHEMSMART protonates it and runs the same HA / A⁻ jobs. See
   :ref:`gaussian-pka-calculations`, :ref:`orca-pka-calculations`, and :ref:`pka-pkb`.

#. **Output analysis** — compute pKa values from completed output files using the backend-independent command
   ``chemsmart run pka``. Add ``--pkb`` / ``--pks`` to also report pKb = pKs − pKa. Analysis is **program-agnostic**:
   the same workflow reads Gaussian ``.log`` and ORCA ``.out`` files and extracts the energies and thermal corrections
   needed for the pKa cycle.

.. toctree::
   :maxdepth: 2
   :caption: pKa Calculations

   gaussian-pka-calculations
   orca-pka-calculations

.. contents:: Table of Contents
   :local:
   :depth: 2

************************
 Execution Architecture
************************

**Job submission**

-  ``chemsmart run/sub gaussian ... pka [submit|batch]`` — prepare and run Gaussian pKa calculations.

-  ``chemsmart run/sub orca ... pka [submit|batch]`` — prepare and run ORCA pKa calculations.

-  Add ``--pkb`` to protonate a free-base input and run the same HA / A⁻ jobs; add ``--pks`` for a non-default
   autoprotolysis constant (see :ref:`pka-pkb`).

-  Use ``chemsmart run`` for local preparation and execution; use ``chemsmart sub`` on HPC clusters to generate
   scheduler scripts (see :ref:`pka-hpc-batch-submission`).

-  A single structure yields one job; batch input (CSV table or multi-molecule CDXML) can produce multiple jobs in one
   invocation.

-  When ``pka`` is invoked without an explicit subcommand, a submission table triggers ``batch``; otherwise ``submit``
   runs.

-  Optional CREST conformational sampling is off by default (``--sampling``). Pass ``-N`` / ``--num-conformers`` on the
   ``pka`` group (after ``pka``, not on ``run`` / ``sub``). ``-N`` is the number of CREST conformers that each receive a
   gas-phase opt+freq and a matching solvent single-point. It is not ``-n`` / ``--num-cores``. See
   :ref:`pka-crest-sampling`.

**Output analysis**

-  ``chemsmart run pka analyze`` — single-system analysis from up to eight output files.
-  ``chemsmart run pka batch-analyze`` — table-driven batch analysis.
-  Add ``--pkb`` and/or ``--pks`` on the ``pka`` group to also report pKb = pKs − pKa (see :ref:`pka-pkb`).
-  Both commands use the same pKa analysis workflow. Gaussian ``.log`` and ORCA ``.out`` files can be mixed in the same
   batch table; each file is read and interpreted on its own.
-  Analysis never invokes ``gaussian`` or ``orca`` job submission — only reads completed output files.

Theory
======

The pKa of an acid HA in aqueous solution is defined by the equilibrium:

.. math::

   \text{HA}_{(\text{aq})} \rightleftharpoons \text{A}^{-}_{(\text{aq})} + \text{H}^{+}_{(\text{aq})}

The pKa is related to the standard Gibbs free energy change by:

.. math::

   \text{p}K_{\text{a}} = \frac{\Delta G^{\circ}_{\text{aq}}}{2.303 \cdot R \cdot T}

where :math:`R` is the gas constant and :math:`T` is the temperature.

**********************
 Thermodynamic Cycles
**********************

CHEMSMART supports two thermodynamic cycles for pKa calculations:

**1. Proton Exchange (Isodesmic) Cycle** (Default, Recommended)

This is the default when ``-s`` is omitted for both job submission and output analysis. A reference acid HRef with known
experimental pKa is required to cancel systematic errors:

.. math::

   \text{HA} + \text{Ref}^{-} \rightarrow \text{A}^{-} + \text{HRef}

The pKa is computed as:

.. math::

   \text{p}K_{\text{a}}(\text{HA}) = \text{p}K_{\text{a}}(\text{HRef}) + \frac{\Delta G_{\text{soln}}}{2.303 \cdot R \cdot T}

where:

.. math::

   \Delta G_{\text{soln}} = \left[ G(\text{A}^{-})_{\text{soln}} + G(\text{HRef})_{\text{soln}} \right] - \left[ G(\text{HA})_{\text{soln}} + G(\text{Ref}^{-})_{\text{soln}} \right]

**2. Direct Cycle**

Uses :math:`G_{\text{soln}}(\text{H}^{+})` in the aqueous dissociation cycle:

.. math::

   \text{p}K_{\text{a}} = \frac{G(\text{A}^{-})_{\text{aq}} + G_{\text{soln}}(\text{H}^{+}) - G(\text{HA})_{\text{aq}}}{2.303 \cdot R \cdot T}

Kelly, Cramer, and Truhlar :math:`\Delta G^{*}_{\text{solv}}(\text{H}^{+}) = -265.9` kcal/mol is the solvation free
energy of the proton, not :math:`G_{\text{soln}}(\text{H}^{+})`. If ``-dG`` is omitted,
:math:`G_{\text{soln}}(\text{H}^{+})` is computed as :math:`G^{\circ}_{\text{gas}}(\text{H}^{+}) + RT \ln(RT /
P^{\circ}) + \Delta G^{*}_{\text{solv}}(\text{H}^{+})` (approximately :math:`-270.3` kcal/mol at 298.15 K; aqueous water
only). Pass ``-dG`` to override.

*********************
 Dual-Level Approach
*********************

CHEMSMART implements a dual-level approach:

#. **Thermal corrections** (:math:`G_{\text{corr}}`) from gas-phase frequency calculations using quasi-harmonic Gibbs
   free energy:

   .. math::

      G_{\text{corr}} = G_{\text{qh}}(T) - E_{\text{gas}}

#. **Solvent energies** (:math:`E_{\text{solv}}`) from high-level single-point calculations in implicit solvent (e.g.,
   SMD or CPCM).

#. **Total free energy in solution**:

   .. math::

      G_{\text{soln}} = E_{\text{solv}} + G_{\text{corr}}

#. **Ensemble free energy** (when a species has more than one conformer):

   .. math::

      G_{\text{eff}} = -RT \ln \sum_i \exp(-G_{\text{soln},i}/RT)

   Each conformer still has :math:`G_{\text{soln},i} = E_{\text{solv},i} + G_{\text{corr},i}`. Analysis replaces that
   species’ :math:`G_{\text{soln}}` with :math:`G_{\text{eff}}` (evaluated with a log-sum-exp shift). A single conformer
   is unchanged.

.. note::

   All internal energies are stored in Hartree (au). :math:`\Delta G_{\text{soln}}` (proton exchange) and :math:`\Delta
   G_{\text{diss}}` (direct dissociation) are converted to kcal/mol for the pKa formula (1 Hartree = 627.5094740631
   kcal/mol).

**********************************
 Job Submission (Gaussian / ORCA)
**********************************

Job submission is backend-specific. Use the dedicated pages for full examples and parameter tables:

-  :ref:`gaussian-pka-calculations`
-  :ref:`orca-pka-calculations`

**Commands**

The default scheme is **proton exchange**, which requires a reference acid (``-r``, ``-rpi``, ``-rc``, ``-rm``). Use
``-s direct`` only when you want the direct dissociation cycle without a reference acid.

.. code:: bash

   # Proton exchange (default) — reference acid required
   chemsmart run gaussian -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   chemsmart run orca -p my_project -f acid.xyz -c 0 -m 1 pka \
       -pi 10 -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   # Direct cycle — no reference acid; must set -s direct explicitly
   chemsmart run gaussian -p my_project -f acid.xyz -c 0 -m 1 pka -pi 10 -s direct
   chemsmart run orca -p my_project -f acid.xyz -c 0 -m 1 pka -pi 10 -s direct

   # Opt-in CREST sampling (N=1 uses crest_best.xyz; -N must follow pka)
   chemsmart run gaussian -p my_project -f acid.xyz -c 0 -m 1 pka \
       --sampling -N 3 -pi 10 -s direct

   # pKb submit: input is the free base B; omit -pi when SMARTS finds one site
   chemsmart run gaussian -p my_project -f pyridine.xyz -c 0 -m 1 pka \
       --pkb -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   # Batch submission (proton exchange requires reference options on the pka group)
   chemsmart run gaussian -p my_project -f pka_input.csv pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1 batch

   # Batch submission (direct cycle)
   chemsmart run gaussian -p my_project -f pka_input.csv pka -s direct batch

**Submission input table** (``pka batch``)

Comma- or whitespace-delimited table with columns ``filepath``, ``proton_index``, ``charge``, ``multiplicity``. With
``--pkb``, ``proton_index`` is the **basic-atom** index to protonate. Leave ``proton_index`` blank to use ChemDraw
colour (CDXML) or a unique SMARTS match (see :ref:`pka-site-resolution`).

.. _pka-crest-sampling:

****************************************
 Optional CREST Conformational Sampling
****************************************

CREST sampling is **opt-in**. Without ``--sampling``, or with ``-N 1``, each species gets one gas-phase opt+freq and one
solvent single-point. Those jobs keep legacy filenames and do not add ``_c1``.

**Flags** (on the ``pka`` group, after ``pka``):

-  ``--sampling`` / ``--no-sampling`` — run CREST on HA and A⁻ before DFT. Default: off. When a reference acid is
   configured, HRef and Ref⁻ are sampled as well.
-  ``-N`` / ``--num-conformers`` — number of lowest-energy CREST conformers to take into DFT (must be ``>= 1``; default
   ``1``). Values greater than 1 require ``--sampling``.

Place ``-N`` **after** ``pka``.

**Geometry selection**

-  ``N = 1``: use ``crest_best.xyz`` (else the first frame of energy-sorted ``crest_conformers.xyz``). DFT filenames
   stay the legacy names without ``_c1``.
-  ``N > 1``: use the first N frames of energy-sorted ``crest_conformers.xyz``.

CREST is gas-phase unless a CREST project YAML (the same ``-p`` name as Gaussian/ORCA, when present) already sets a
solvent. The HA / A⁻ pair is built **before** sampling, so atom order after CREST does not matter.

**Workflow**

For each sampled species (HA, A⁻, and HRef / Ref⁻ when a reference acid is set):

.. code:: text

   CREST conformers
       -> N gas-phase optimization/frequency calculations
       -> N solvent single-point calculations
       -> N conformer solution free energies
       -> one ensemble effective free energy

Each selected conformer receives its own gas-phase opt+freq job and a solvent single-point on that optimized geometry.
Gas-phase and solvent outputs must form a matching pair: the same count, in conformer order (``c1``, ``c2``, …,
``c10``). A mismatch raises an error that names the species and both counts. Analysis does not drop later conformers
when ensemble files exist.

Parent completion requires DFT opt+SP only; CREST is best-effort.

If CREST has not finished, the serial pKa job waits (same as opt/SP phases) so HPC resubmits can continue. After CREST
terminates abnormally or completes without usable geometries, CHEMSMART logs a warning and falls back to the input
geometry (or the shorter available conformer set).

**Analysis**

For each conformer, analysis computes :math:`G_{\text{soln},i} = E_{\text{solv},i} + G_{\text{corr},i}`. When a species
has more than one conformer, that species’ free energy in the pKa cycle is the ensemble :math:`G_{\text{eff}}`. A single
conformer is unchanged. See Dual-Level Approach above.

.. _pka-site-resolution:

***********************************************
 Ionizable Site Resolution (pKa and ``--pkb``)
***********************************************

Job submission locates a single ionizable site in this order. The same sequence applies to ordinary pKa and to
``--pkb``:

#. **Explicit** ``-pi`` / ``--proton-index`` (always wins).

   -  pKa: 1-based index of the **hydrogen to remove**.
   -  ``--pkb``: 1-based index of the **heavy atom to protonate**. A hydrogen index is rejected.

#. **ChemDraw colour** (``.cdxml`` / ``.cdx`` only).

   -  pKa: uniquely coloured **hydrogen**.
   -  ``--pkb``: uniquely coloured **non-hydrogen** (the basic atom).
   -  Use ``-cc`` / ``--color-code`` when more than one non-majority colour exists.
   -  This step is skipped when the file is not CDXML, or when every atom shares one colour (no colour markup).
   -  Ambiguous colour markup (several uniquely coloured candidates) is an error. SMARTS is **not** used as a fallback
      in that case.

#. **RDKit SMARTS** — exactly one matching site. Zero or several matches raise an error asking for ``-pi`` or ``-cc``.

SMARTS patterns are conservative:

-  **pKa:** carboxylic acid H, phenol H, thiol H, and ammonium/iminium H. Generic alcohols are not selected.
-  **``--pkb``:** neutral nitrogen bases (aliphatic and aromatic amines, pyridine-like ring nitrogen, anilines). Amides,
   nitro groups, nitriles, and quaternary nitrogen are excluded.

.. note::

   Colour a hydrogen for pKa. For ``--pkb``, colour the **basic heavy atom** (for example the pyridine nitrogen) the
   same way. Colouring a hydrogen while ``--pkb`` is set is an error: colour the basic atom, or drop ``--pkb`` and run
   pKa.

***************************************
 ChemDraw CDXML / CDX Input (pKa Jobs)
***************************************

pKa job submission can read structures directly from ChemDraw ``.cdxml`` and ``.cdx`` files. CHEMSMART reads atom
colours in the drawing to identify the ionizable site (see :ref:`pka-site-resolution`). For pKa, colour the **acidic
proton** (or the ``H`` in a functional group such as –OH) with a distinct colour. For ``--pkb``, colour the **basic
heavy atom** instead. CHEMSMART auto-detects a uniquely coloured site so ``-pi`` is often unnecessary.

**Single-molecule submit**

When ``-f`` points to one CDXML structure with a single fragment, omit ``-pi`` if the coloured site is unique:

.. code:: bash

   chemsmart run gaussian -p my_project -f phenol.cdxml -c 0 -m 1 pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   chemsmart run gaussian -p my_project -f phenol.cdxml -c 0 -m 1 pka -s direct

   # pKb: colour the pyridine nitrogen (not a hydrogen)
   chemsmart run gaussian -p my_project -f pyridine.cdxml -c 0 -m 1 pka \
       --pkb -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

Use ``-cc`` / ``--color-code`` when several atoms share similar styling and you need to select a specific ChemDraw
colour-table index. The reference acid may also be a CDXML file; in that case ``-rpi`` can be omitted when the reference
proton is uniquely coloured (or use ``-rcc`` / ``--reference-color-code``). If the CDXML file has **no colour markup**
(all atoms the same colour), site resolution falls through to SMARTS, the same as for XYZ.

**Multi-molecule CDXML (one job per fragment)**

A single ``.cdxml`` / ``.cdx`` file may contain **multiple molecules** (multiple ChemDraw fragments). CHEMSMART performs
**per-fragment** colour detection and creates **one pKa job per fragment**. With ``--pkb``, each fragment is protonated
at its coloured (or SMARTS) basic atom.

Pass the file with ``pka batch`` (or ``pka submit`` for a single-fragment file):

.. code:: bash

   chemsmart run gaussian -p my_project -f acids.cdxml -c 0 -m 1 pka \
       -r ref_acid.xyz -rpi 21 -rc 1 -rm 1 batch

   chemsmart run orca -p my_project -f acids.cdxml -c 0 -m 1 pka -s direct batch

Job labels are derived from the filename, e.g. ``acids_frag1_pka`` (Gaussian) or ``acids_frag1_pka`` (ORCA).

**Charge and multiplicity**

How charge and multiplicity are resolved depends on the input mode:

.. list-table::
   :header-rows: 1
   :widths: 35 65

   -  -  Input mode
      -  Charge / multiplicity source

   -  -  CSV batch table
      -  **Required columns** on each row (``charge``, ``multiplicity``). The parent ``gaussian`` / ``orca`` command
         does not need ``-c`` / ``-m`` when ``-f`` is a table.

   -  -  Multi-fragment CDXML (``-f`` is ``.cdxml`` / ``.cdx``)
      -  Parent ``-c`` / ``-m`` apply to every fragment by default. If either is omitted on the backend command,
         CHEMSMART may copy values from the parsed CDXML structure when the drawing supplies them.

   -  -  Single-molecule submit (XYZ, LOG, CDXML, …)
      -  ``-c`` and ``-m`` on the backend command are required unless already present on merged project/job settings.

Blank ``proton_index`` in a table row uses the site order in :ref:`pka-site-resolution` (ChemDraw colour, then SMARTS).
It does **not** remove the requirement for ``charge`` and ``multiplicity`` columns in CSV tables.

**CDXML paths inside a CSV batch table**

You can mix XYZ and CDXML inputs in the same submission table. Each row becomes one pKa job. Leave ``proton_index``
blank to auto-detect the site (coloured atom on CDXML; unique SMARTS match on XYZ and uncoloured CDXML). An explicit
value overrides detection. With ``--pkb``, a blank or explicit ``proton_index`` is the **basic-atom** index.

.. code:: text

   filepath,proton_index,charge,multiplicity
   /path/to/acid1.xyz,12,0,1
   /path/to/acid2.cdxml,,0,1

Multi-molecule CDXML files cannot be expanded from a table row. Pass them directly as ``-f`` with ``pka batch`` (see
above) to create one job per ChemDraw fragment.

.. note::

   If ``-f`` is a CDXML file (not a CSV table), CHEMSMART routes to per-fragment site detection automatically. For
   general CDXML structure handling outside pKa, see :doc:`chemdraw-organometallic`.

.. _pka-chemdraw-preview:

******************
 ChemDraw preview
******************

``--preview`` reads a structure file, prints one row per ChemDraw fragment, and stops before creating or submitting
Gaussian, ORCA, or CREST jobs. It does not write job directories.

.. code:: bash

   chemsmart run gaussian -p my_project -f acids.cdxml -c 0 -m 1 \
       pka --preview -s direct batch

   chemsmart run orca -p my_project -f acids.cdxml -c 0 -m 1 \
       pka --preview -s direct batch

The table columns are fragment number, job label, mode (``pKa`` or ``pKb``), selected site (atom index and element),
selection source, input charge, input multiplicity, and status. Fragment order follows the ChemDraw document.

-  ``-pi`` / ``--proton-index`` takes precedence and is reported as ``explicit``.
-  A uniquely coloured site is reported as ``ChemDraw colour``.
-  An uncoloured drawing with one SMARTS match is reported as ``SMARTS``.

For ``pKa``, the site is the hydrogen that will be removed. For ``--pkb``, the site is the heavy atom that will be
protonated, and the status also gives the index of the hydrogen added to the free base. Ambiguous colour markup or more
than one SMARTS site stops the command before any job is created. The charge and multiplicity columns are the values
submission would use: parent ``-c`` / ``-m`` when they are set, otherwise the values read from the structure.

**Proton and reference options for CDXML**

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-pi``
      -  ``--proton-index``
      -  Optional when a uniquely coloured ChemDraw site or a unique SMARTS match is present. pKa: hydrogen to remove.
         ``--pkb``: heavy atom to protonate.

   -  -  ``-cc``
      -  ``--color-code``
      -  ChemDraw colour-table index for the target site (``.cdxml`` / ``.cdx`` only). pKa: acidic proton. ``--pkb``:
         basic heavy atom.

   -  -  ``-rpi``
      -  ``--reference-proton-index``
      -  Optional when ``-r`` is a CDXML file with a uniquely coloured reference proton.

   -  -  ``-rcc``
      -  ``--reference-color-code``
      -  ChemDraw colour-table index for the reference proton (``.cdxml`` / ``.cdx`` reference only).

**Job output file naming**

Each pKa job creates gas-phase opt+freq and solvent single-point sub-jobs. Output filenames follow the sub-job label:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   -  -  Sub-job label suffix
      -  Species
   -  -  ``_HA_opt``
      -  Target acid HA (gas-phase opt+freq)
   -  -  ``_A_opt``
      -  Target conjugate base A⁻ (gas-phase opt+freq)
   -  -  ``_HA_sp``
      -  HA solvent single-point
   -  -  ``_A_sp``
      -  A⁻ solvent single-point
   -  -  ``_HRef_opt`` / ``_Ref_opt``
      -  Reference acid / conjugate base (proton exchange only)
   -  -  ``_HRef_sp`` / ``_Ref_sp``
      -  Reference solvent single-points (proton exchange only)
   -  -  ``_HA_crest`` / ``_A_crest``
      -  CREST sampling of HA / A⁻ (``--sampling`` only)
   -  -  ``_HA_opt_cN`` / ``_A_opt_cN``
      -  N-th CREST conformer gas-phase opt+freq when ``-N`` > 1 (1-based)
   -  -  ``_HA_sp_cN`` / ``_A_sp_cN``
      -  N-th CREST conformer solvent single-point when ``-N`` > 1 (1-based)
   -  -  ``_HRef_opt_cN`` / ``_Ref_opt_cN``
      -  N-th reference conformer gas-phase opt+freq when ``-N`` > 1
   -  -  ``_HRef_sp_cN`` / ``_Ref_sp_cN``
      -  N-th reference conformer solvent single-point when ``-N`` > 1

Gaussian batch jobs use the input stem as the job label (e.g. ``acid1_HA_opt.log``). ORCA batch jobs append ``_pka`` to
the stem (e.g. ``acid1_pka_HA_opt.out``). The output-analysis autodiscovery convention below is aligned with the
``{basename}_pka_*`` pattern used by ORCA submission and by typical batch output tables. With ``--pkb``, HA is BH⁺ and
A⁻ is B; file names are unchanged. With ``--sampling`` and ``-N 1``, DFT names stay the legacy forms above and do not
gain ``_c1``. With ``-N`` greater than 1, every gas-phase opt and every solvent single-point label gains ``_c1``,
``_c2``, … so each conformer has a matching gas/solvent pair. The same suffixes apply to HRef and Ref⁻.

.. _pka-pkb:

**************************
 pKb from pKa (``--pkb``)
**************************

There is no separate ``pkb`` command. Internally HA is always BH⁺ and A⁻ is always B. ``--pkb`` has two independent
uses:

**Submit** — the input file is the **free base** B. CHEMSMART protonates the resolved basic site, then runs the existing
pKa jobs. ``-c`` / ``-m`` remain the charge and multiplicity of the file passed in. After protonation, the HA (BH⁺)
charge is the input charge plus one; the conjugate-base (B) charge still defaults to the original ``-c``.

Proton-exchange submit still requires a reference acid (``-r``, ``-rpi``, ``-rc``, ``-rm``, …). The experimental value
is the **pKa** of HRef in that solvent. For amines, HRef should be a similar cationic acid. Direct-cycle submit (``-s
direct``) is still allowed; pKb conversion does not depend on how pKa was obtained.

**Analyze** — ``--pkb`` and/or ``--pks`` convert a computed pKa after the usual analysis step:

.. math::

   \mathrm{p}K_{\mathrm{b}} = \mathrm{p}K_{\mathrm{s}} - \mathrm{p}K_{\mathrm{a}}

This works on any completed HA / A⁻ outputs, including jobs that were not submitted with ``--pkb``.

**pKs**

``--pks`` is the solvent autoprotolysis constant. If ``--pkb`` is set and ``--pks`` is omitted, ``14.0`` is used. That
default is the conventional aqueous value (approximately :math:`\mathrm{p}K_{\mathrm{w}}` near 25 °C). It is **not**
valid for other solvents, which have different autoprotolysis constants. Supplying ``--pks`` without ``--pkb`` on
analysis also enables pKb reporting.

.. warning::

   Default pKs = 14.0 is for aqueous water only. When that default is applied and ``-si`` / ``--solvent-id`` is not
   water (or ``h2o``), CHEMSMART warns: pass ``--pks`` for a literature non-aqueous autoprotolysis constant.

.. code:: bash

   # Submit free base; default pKs = 14 (water)
   chemsmart run gaussian -p my_project -f pyridine.xyz -c 0 -m 1 pka \
       --pkb -pi 1 -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   # Non-aqueous: supply literature pKs (warning if default 14 is used)
   chemsmart run gaussian -p my_project -f pyridine.xyz -c 0 -m 1 -si acetonitrile pka \
       --pkb --pks 33.3 -r ref_acid.xyz -rpi 21 -rc 1 -rm 1

   # Analysis: print pKb from existing HA/A- outputs
   chemsmart run pka --pkb analyze \
       -ha acid1_pka_HA_opt.log -hr ref_acid_pka_HRef_opt.log -rp 6.75

   chemsmart run pka --pkb --pks 16.7 analyze \
       -ha acid1_pka_HA_opt.log -hr ref_acid_pka_HRef_opt.log -rp 6.75

.. _pka-hpc-batch-submission:

********************************************
 HPC Cluster Submission (``chemsmart sub``)
********************************************

On a cluster, use ``chemsmart sub`` instead of ``chemsmart run`` to write scheduler scripts and per-job run wrappers.
The pKa workflow is unchanged at the chemistry level; only the launch path differs.

**Single job**

.. code:: bash

   chemsmart sub gaussian -p my_project -f acid.xyz -c 0 -m 1 pka -pi 10 -s direct submit

**Batch table or multi-fragment CDXML**

One scheduler submission is created **per table row** or **per ChemDraw fragment**. CHEMSMART expands ``pka batch`` into
multiple jobs locally, then writes a separate ``chemsmart_sub_<label>.sh`` and ``chemsmart_run_<label>.py`` for each
job.

.. code:: bash

   chemsmart sub gaussian -p my_project -f pka_input.csv pka -s direct batch

   chemsmart sub gaussian -p my_project -f acids.cdxml -c 0 -m 1 pka -s direct batch

Per-job script reconstruction
=============================

Each cluster run wrapper must replay **one** pKa submission, not the entire batch table or full multi-fragment CDXML
file. When a job is created from ``pka batch``, CHEMSMART stores row- or fragment-level metadata and rewrites the CLI
inside ``chemsmart_run_<label>.py`` before submission:

#. **CSV batch rows** — replace the table path in ``-f`` / ``--filename`` with that row's ``filepath``; change ``batch``
   to ``submit``.
#. **Multi-fragment CDXML** — point ``-f`` at the same CDXML file but add ``--index`` / ``-i`` so only one fragment is
   processed; change ``batch`` to ``submit``.
#. **Explicit per-job options** — inject or update ``--proton-index``, ``--charge``, ``--multiplicity``, and ``--label``
   so the reconstructed command is self-contained and passes Click validation on the cluster node.
#. **Proton exchange tables** — rows after the first may run with ``-s direct``; reference-acid flags are dropped from
   later rows automatically (same behaviour as local ``pka batch``).

Example: a two-row CSV batch submitted with ``chemsmart sub ... pka batch`` yields two run scripts. The script for row
two might equivalent to:

.. code:: bash

   chemsmart run gaussian -p my_project -f /path/to/acid2.xyz -c 0 -m 1 \
       pka -pi 8 -s direct submit

Example: a five-fragment CDXML file ``pka_scale.cdxml`` yields labels such as ``pka_scale_frag1_pka``, …,
``pka_scale_frag5_pka``. Each run script targets one fragment via ``--index`` and the matching ``--proton-index``,
``--label``, ``-c``, and ``-m``.

This reconstruction is what allows ``chemsmart sub ... pka batch`` on a cluster to behave like five independent ``pka
submit`` calls while you only maintain one top-level submission command locally.

See also :doc:`cli-overview` for general ``chemsmart sub`` usage and :doc:`configuration-server-settings` for scheduler
configuration.

*****************************************
 Output Analysis (``chemsmart run pka``)
*****************************************

All post-processing lives under ``chemsmart run pka``. No Gaussian or ORCA backend is invoked during analysis.

Thermochemistry extraction
==========================

For each output file, analysis:

#. Open the file and detect the program (Gaussian or ORCA).
#. Reads the gas-phase SCF energy and quasi-harmonic Gibbs free energy (for opt+freq outputs).
#. Reads the solvent-phase SCF energy (for single-point outputs).
#. Raises a clear error if a required quantity cannot be extracted.

Computing pKa from Output Files (``analyze``)
=============================================

If you have completed output files for a single acid, compute pKa with ``analyze``.

**Proton exchange (default)**

Only ``-ha`` and ``-hr`` are strictly required; the remaining six companion files are auto-discovered when they follow
the naming convention below. ``-rp`` / ``--reference-pka`` is required.

.. code:: bash

   chemsmart run pka analyze \
       -ha acid1_pka_HA_opt.log \
       -hr ref_acid_pka_HRef_opt.log \
       -rp 6.75 \
       -T 333.15 -c 1.0 -csg 100 -ch 100

Provide all eight files explicitly when auto-discovery is not appropriate:

.. code:: bash

   chemsmart run pka analyze \
       -ha acid1_pka_HA_opt.log \
       -a acid1_pka_A_opt.log \
       -hr ref_acid_pka_HRef_opt.log \
       -r ref_acid_pka_Ref_opt.log \
       -has acid1_pka_HA_sp.log \
       -as acid1_pka_A_sp.log \
       -hrs ref_acid_pka_HRef_sp.log \
       -rs ref_acid_pka_Ref_sp.log \
       -rp 6.75 \
       -T 298.15

**Direct dissociation**

Four output files are required (HA, A⁻, and their solvent single-points). Specify ``-s direct`` on the ``pka`` group
**before** the ``analyze`` subcommand. Omit ``-dG`` to use the computed aqueous :math:`G_{\text{soln}}(\text{H}^{+})`
default, or pass ``-dG`` to override:

.. code:: bash

   chemsmart run pka -s direct analyze \
       -ha acid1_pka_HA_opt.log \
       -T 298.15

Only ``-ha`` is strictly required; ``-a``, ``-has``, and ``-as`` are auto-discovered from the target-acid suffix
convention when omitted.

**pKb conversion**

``--pkb`` (and/or ``--pks``) reports :math:`\mathrm{p}K_{\mathrm{b}} = \mathrm{p}K_{\mathrm{s}} -
\mathrm{p}K_{\mathrm{a}}` from the same HA / A⁻ files. See :ref:`pka-pkb`.

.. code:: bash

   chemsmart run pka --pkb analyze \
       -ha acid1_pka_HA_opt.log \
       -hr ref_acid_pka_HRef_opt.log \
       -rp 6.75

   chemsmart run pka --pkb --pks 16.7 analyze \
       -ha acid1_pka_HA_opt.log \
       -hr ref_acid_pka_HRef_opt.log \
       -rp 6.75

File autodetection (``analyze``)
================================

When companion paths are omitted, CHEMSMART derives them from the HA and HRef gas-phase files using the same suffix
patterns as ``batch-analyze``:

**From the HA gas-phase file (``-ha``)**

-  ``<basename>_pka_A_opt.<ext>`` — conjugate base gas-phase
-  ``<basename>_pka_HA_sp.<ext>`` — HA solvent single-point
-  ``<basename>_pka_A_sp.<ext>`` — conjugate base solvent SP

If ``<basename>_pka_HA_opt_c*.<ext>`` (or ``.out``) files exist, auto-discovery uses that sorted ensemble and the
matching ``_A_opt_c*``, ``_HA_sp_c*``, and ``_A_sp_c*`` files instead of the single-file suffixes. Otherwise the
single-file suffixes above are used. ``N = 1`` keeps those legacy suffixes and does not add ``_c1``.

Gas-phase and solvent outputs for each species must form a matching pair: the same number of files, in conformer order
(``c1``, ``c2``, …, ``c10``). Analysis rejects unequal counts and names the species and both counts.

**From the HRef gas-phase file (``-hr``)**

-  ``<basename>_pka_Ref_opt.<ext>`` — reference conjugate base
-  ``<basename>_pka_HRef_sp.<ext>`` — reference acid solvent SP
-  ``<basename>_pka_Ref_sp.<ext>`` — reference conjugate base solvent SP

Alternative suffixes (``_HRef_opt``, ``_pka_cb``, etc.) are also recognised. The file extension (``.log`` or ``.out``)
is chosen from the detected program. Override any auto-discovered path with the corresponding flag.

If a required file is missing or cannot be parsed, analysis stops with a clear error (missing paths or missing
thermochemistry data).

Batch Processing of Output Files (``batch-analyze``)
====================================================

Parse a table of pre-computed output file paths to calculate pKa values in batch.

**Proton exchange (default)**

.. code:: bash

   chemsmart run pka -T 333.15 -c 1.0 -csg 100 -ch 100 batch-analyze \
       -o pka_output_table.csv \
       -O results.dat

**Direct dissociation**

Specify ``-s direct``. Omit ``-dG`` to use the computed aqueous default:

.. code:: bash

   chemsmart run pka -s direct batch-analyze \
       -o pka_output_table_direct.csv \
       -O results_direct.dat

**pKb conversion**

Add ``--pkb`` and optionally ``--pks`` on the ``pka`` group. The summary table gains ``pKb`` and ``pKs`` columns.

.. code:: bash

   chemsmart run pka --pkb batch-analyze -o pka_output_table.csv -O results.dat
   chemsmart run pka --pkb --pks 16.7 batch-analyze -o pka_output_table.csv

The formatted batch summary table is printed to stdout. When ``-O`` / ``--output-results`` is given, the same formatted
report is written to that file (not a wide CSV of input columns).

Output table format
-------------------

The output table (``-o`` / ``--output-table``) must contain at least a ``basename`` column. Other file paths may be
omitted and are auto-discovered from ``basename`` when blank.

**Required column**

-  ``basename``: Unique identifier for each acid. Used for file auto-discovery.

**Target-acid columns (both schemes)**

Accepted header aliases include ``ha_opt``, ``a_opt``, ``ha_solv``, ``a_solv``, etc. (see column list below).

When blank, CHEMSMART searches for ``<basename><suffix>.<ext>`` in the current working directory. Suffixes are tried in
order; both ``.log`` and ``.out`` are tested.

.. list-table::
   :header-rows: 1
   :widths: 20 50 30

   -  -  Column
      -  Description
      -  Auto-discovery suffixes (first match wins)

   -  -  ``ha_gas``
      -  HA gas-phase opt+freq output
      -  ``_pka_HA_opt_c*``, ``_pka_HA_opt``, ``_pka_HA``, ``_pka``

   -  -  ``a_gas``
      -  A⁻ gas-phase opt+freq output
      -  ``_pka_A_opt_c*``, ``_pka_A_opt``, ``_pka_A``, ``_pka_cb``

   -  -  ``ha_sp``
      -  HA solvent single-point output
      -  ``_pka_HA_sp_c*``, ``_pka_HA_sp``, ``_pka_sp``

   -  -  ``a_sp``
      -  A⁻ solvent single-point output
      -  ``_pka_A_sp_c*``, ``_pka_A_sp``, ``_pka_cb_sp``

**Reference-acid columns (proton exchange only)**

Ignored for direct dissociation. Blank reference columns inherit values from the previous row.

.. list-table::
   :header-rows: 1
   :widths: 20 80

   -  -  Column
      -  Description
   -  -  ``href_gas``
      -  HRef gas-phase opt+freq output (not auto-discovered from ``basename``; provide explicitly or inherit)
   -  -  ``ref_gas``
      -  Ref⁻ gas-phase opt+freq output
   -  -  ``href_sp``
      -  HRef solvent single-point output
   -  -  ``ref_sp``
      -  Ref⁻ solvent single-point output
   -  -  ``pka_ref``
      -  Experimental pKa of the reference acid

**Example** ``pka_output_table.csv``:

.. code:: text

   basename,ha_gas,a_gas,ha_sp,a_sp,href_gas,ref_gas,href_sp,ref_sp,pka_ref
   phenol,,,,,ref_acid_pka_HRef_opt.log,ref_acid_pka_Ref_opt.log,ref_acid_pka_HRef_sp.log,ref_acid_pka_Ref_sp.log,6.75
   benzoic_acid,,,,,,,,,6.75

``batch-analyze`` options
-------------------------

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-o``
      -  ``--output-table``
      -  **Required.** Path to the output-file table.

   -  -  ``-O``
      -  ``--output-results``
      -  Optional path for the formatted results report. Stdout always receives the summary table.

   -  -  ``-p``
      -  ``--program``
      -  Require every populated output path to match ``gaussian`` or ``orca``. Default: ``auto`` (per-file detection;
         supports mixed tables).

Mixed Gaussian / ORCA tables
----------------------------

With ``-p auto`` (default), the program behind each output file is detected automatically. A batch table may contain
Gaussian ``.log`` targets and ORCA ``.out`` reference files in the same run. Use ``-p gaussian`` or ``-p orca`` only
when you want to **validate** that all populated paths belong to one backend.

*************************
 Analysis Scheme Options
*************************

These options apply to ``chemsmart run pka`` (``analyze`` and ``batch-analyze``). They are separate from submission
options on ``chemsmart run/sub gaussian ... pka`` and ``chemsmart run/sub orca ... pka``.

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-s``
      -  ``--scheme``
      -  Thermodynamic cycle: ``direct`` or ``proton exchange``. Default: ``proton exchange``.

   -  -  ``-dG``

      -  ``--delta-g-proton``

      -  :math:`G_{\text{soln}}(\text{H}^{+})` override in kcal/mol for the direct cycle. If omitted, a T-dependent
         aqueous default is computed from Kelly, Cramer, and Truhlar :math:`\Delta G^{*}_{\text{solv}}(\text{H}^{+}) =
         -265.9` kcal/mol.

   -  -
      -  ``--pkb``
      -  Also print pKb = pKs − pKa. If ``--pks`` is omitted, 14.0 is used.

   -  -
      -  ``--pks``
      -  Solvent autoprotolysis constant. Supplying ``--pks`` also enables pKb reporting.

.. note::

   If ``-dG`` is omitted for the direct cycle, :math:`G_{\text{soln}}(\text{H}^{+})` is computed from Kelly, Cramer, and
   Truhlar :math:`\Delta G^{*}_{\text{solv}}(\text{H}^{+}) = -265.9` kcal/mol (aqueous water at the requested
   temperature). Pass ``-dG`` to override.

.. warning::

   Default pKs = 14.0 is for aqueous water only. When ``--pkb`` is used without ``--pks`` and the solvent is not water,
   CHEMSMART warns and you should pass a literature ``--pks``. See :ref:`pka-pkb`.

*********************
 Output File Options
*********************

Used by ``analyze`` (not ``batch-analyze``, which reads paths from the table).

**Gas-phase optimization + frequency files**

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-ha``
      -  ``--ha``
      -  HA gas-phase opt+freq output.

   -  -  ``-a``
      -  ``--a``
      -  A⁻ gas-phase opt+freq output.

   -  -  ``-hr``
      -  ``--href``
      -  HRef gas-phase opt+freq output.

   -  -  ``-r``
      -  ``--ref``
      -  Ref⁻ gas-phase opt+freq output.

**Solvent single-point files**

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-has``
      -  ``--ha-solv``
      -  HA solvent single-point output.

   -  -  ``-as``
      -  ``--a-solv``
      -  A⁻ solvent single-point output.

   -  -  ``-hrs``
      -  ``--href-solv``
      -  HRef solvent single-point output.

   -  -  ``-rs``
      -  ``--ref-solv``
      -  Ref⁻ solvent single-point output.

*************************
 Thermochemistry Options
*************************

Shared by ``analyze`` and ``batch-analyze``.

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   -  -  Short
      -  Long
      -  Description

   -  -  ``-T``
      -  ``--temperature``
      -  Temperature in Kelvin. Default: ``298.15`` K.

   -  -  ``-c``
      -  ``--concentration``
      -  Concentration in mol/L. Default: ``1.0`` mol/L.

   -  -  ``-P``
      -  ``--pressure``
      -  Pressure in atm. Default: ``1.0`` atm.

   -  -  ``-csg``
      -  ``--cutoff-entropy-grimme``
      -  Cutoff frequency (cm⁻¹) for entropy using Grimme's quasi-RRHO. Default: ``100.0``.

   -  -  ``-cst``
      -  ``--cutoff-entropy-truhlar``
      -  Cutoff frequency (cm⁻¹) for entropy using Truhlar's quasi-RRHO. Mutually exclusive with ``-csg``.

   -  -  ``-ch``
      -  ``--cutoff-enthalpy``
      -  Cutoff frequency (cm⁻¹) for enthalpy using Head-Gordon's method. Default: ``100.0``.

   -  -  ``-rp``
      -  ``--reference-pka``
      -  Experimental pKa of HRef. Required for proton exchange analysis.

Output Format
=============

When computing pKa from output files, CHEMSMART prints a detailed summary. The format depends on the analysis scheme.

**Proton exchange**

.. code:: text

   ==============================================================================
   pKa Calculation - Dual-level Proton Exchange Scheme
   ==============================================================================
   Reaction: HA + Ref- -> A- + HRef
   Temperature: 373.15 K

   Method:
     G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)
     G_soln = E_solv + G_corr  (solution free energy)
     DG_soln = [G(A-)_soln + G(HRef)_soln] - [G(HA)_soln + G(Ref-)_soln]
     pKa = pKa_ref + DG_soln / (RT * ln10)
   ------------------------------------------------------------------------------

   Gas-Phase Electronic Energies (E_gas, au):
     HA:    -345.7419436500
     A-:    -344.9153986020
     HRef:  -365.8436493070
     Ref-:  -365.4561783660

   Thermal Corrections (G_corr = qh-G - E_gas, au):
     HA:    0.0931931305
     A-:    0.0758935969
     HRef:  0.1404528844
     Ref-:  0.1267467582

   Solvent Single-Point Energies (E_solv, au):
     HA:    -346.4882221850
     A-:    -345.8989956310
     HRef:  -366.5974351550
     Ref-:  -366.1368369100

   Solution Free Energies (G_soln = E_solv + G_corr, au):
     HA:    -346.3950290545
     A-:    -345.8231020341
     HRef:  -366.4569822706
     Ref-:  -366.0100901518
   ------------------------------------------------------------------------------

   pKa Calculation:
     DG_soln = 0.1250349015 au
             = 78.4606 kcal/mol
     pKa(HRef)_ref = 6.75

     *** Computed pKa(HA) = 52.70 ***
   ==============================================================================

When ``--pkb`` or ``--pks`` is set, the summary adds pKs and pKb lines:

.. code:: text

   pKs = 14.00 (default aqueous)

   *** Computed pKa(HA) = 52.70 ***
   *** Computed pKb(B)  = -38.70 ***

With an explicit ``--pks``, the source reads ``(user-supplied)`` instead of ``(default aqueous)``.

When a species has more than one conformer, the solution-free-energy lines read ``HA (2 conformers, G_eff):`` and the
method block notes :math:`G_{\text{eff}} = -RT \ln \sum \exp(-G_i/RT)`.

**Direct dissociation**

.. code:: text

   ==============================================================================
   pKa Calculation - Direct Dissociation Scheme
   ==============================================================================
   Reaction: HA -> A- + H+
   Temperature: 298.15 K

   Method:
     G_corr = qh-G(T) - E_gas  (from gas-phase freq calculation)
     G_soln = E_solv + G_corr  (solution free energy)
     DG_diss = G_soln(A-) + G_soln(H+) - G_soln(HA)
     pKa = DG_diss / (2.303 * R * T)
   ------------------------------------------------------------------------------

   Gas-Phase Electronic Energies (E_gas, au):
     HA:  -345.7419436500
     A-:  -344.9153986020

   Thermal Corrections (G_corr = qh-G - E_gas, au):
     HA:  0.0931931305
     A-:  0.0758935969

   Solvent Single-Point Energies (E_solv, au):
     HA:  -346.4882221850
     A-:  -345.8989956310

   Solution Free Energies (G_soln = E_solv + G_corr, au):
     HA:  -346.3950290545
     A-:  -345.8231020341
   ------------------------------------------------------------------------------

   pKa Calculation:
     G_soln(H+) = -270.2811 kcal/mol
                  (computed aqueous default for water at 298.15 K)
     DG_diss = 0.1412066215 au
             = 88.6085 kcal/mol

     *** Computed pKa(HA) = 64.95 ***
   ==============================================================================

**Batch analyze output**

``batch-analyze`` prints a compact table whose ΔG column header matches the scheme. The same formatted report is written
to ``-O`` when provided.

**Typical table** (single-acid basenames):

.. code:: text

   ==============================================================================
   Batch pKa Results (Dual-level Proton Exchange)
   ==============================================================================
   Temperature: 298.15 K
   Pressure: 1.0 atm
   basename                              pKa   ΔG_soln (kcal/mol)
   ------------------------------------------------------------------------------
   phenol                               10.12             13.4567
   benzoic_acid                          4.20              5.7890
   ==============================================================================

**Multi-fragment CDXML workflow** — after ``chemsmart sub ... -f pka_scale.cdxml pka batch`` and ``chemsmart run pka
batch-analyze``, basenames match the fragment labels (``<stem>_frag<N>_pka`` or the Gaussian stem without ``_pka``
suffix, depending on backend). Example (values from a test reference acid; not physically meaningful):

.. code:: text

   ==============================================================================
   Batch pKa Results (Dual-level Proton Exchange)
   ==============================================================================
   Temperature: 278.15 K
   Pressure: 1.0 atm
   basename                              pKa   ΔG_soln (kcal/mol)
   ------------------------------------------------------------------------------
   pka_scale_frag1                      -8.17            -23.8913
   pka_scale_frag4                      -1.24            -15.0654
   pka_scale_frag3                       9.22             -1.7521
   pka_scale_frag5                      -2.68            -16.8958
   ==============================================================================

For direct dissociation, the header reads ``Batch pKa Results (Direct Dissociation)`` and the column is labeled
``ΔG_diss (kcal/mol)``. With ``--pkb`` or ``--pks``, the table also includes ``pKb`` and ``pKs`` columns.

References
==========

#. Kelly, C. P.; Cramer, C. J.; Truhlar, D. G. (2006). *J. Phys. Chem. B*, 110, 16066. (Absolute proton solvation
   energy)
#. Grimme, S. (2012). *Chem. Eur. J.*, 18, 9955. (Quasi-RRHO method)
#. Marenich, A. V.; Cramer, C. J.; Truhlar, D. G. (2009). *J. Phys. Chem. B*, 113, 6378. (SMD solvation model)

See Also
========

-  :doc:`thermochemistry-analysis`
