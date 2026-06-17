# Implementation Plan: ORCA Ionic-Crystal QM/MM Input Generation

## Summary of Proposed Design

Add a new preparation path for ionic-crystal QM/MM jobs, driven by `--geninput` on the ORCA QMMM subcommand. When enabled and `--jobtype ionic-crystal-qmmm` is used, CHEMSMART will:

1. Generate a CrystalPrep template input (`ORCA_crystalprep`) from a `.cif`.
2. Generate a submission script that runs ORCA utilities sequentially (CrystalPrep + ORCA_mm), validates expected files, then runs the final ORCA QM/MM input.
3. Preserve the existing final `.inp` generation and submission behavior for all jobs, only adding the prep path when `--geninput` is present.

## Expected CLI Interface

### New flag (ionic crystal only)
- `--geninput`: enable CrystalPrep/ORCA_mm preparation steps.

### CrystalPrep template options (new group)
Suggested flag names (prefix `--cp-` to avoid collisions):

- `--cp-input-cif` (default: use `-f` if it is a `.cif`)
- `--cp-docif` (bool)
- `--cp-dosupercell` (bool)
- `--cp-scdimension` (e.g., `15x15x15`)
- `--cp-doembedding` (bool)
- `--cp-dolayers` (bool)
- `--cp-doicqmmminput` (bool)
- `--cp-atomtype` (repeatable, e.g. `--cp-atomtype Na 0 1.0 0.0`)
- `--cp-qccharge`, `--cp-qcmult`, `--cp-neutralize` (optional)
- `--cp-template-out` (optional template output filename override)

### Final ORCA QM/MM input
Continue to use existing QMMM flags for final `.inp` generation (already implemented):

- `--high-level-functional`, `--high-level-basis`, `--low-level-method`
- `--use-qm-info-from-pdb`, `--use-qm3-info-from-pdb`
- `--ecp-layer-ecp`, `--conv-charges`, `--enforce-total-charge`
- `--charge-total`, `--print-level`, `--charge-high`, `--mult-high`

## Files / Classes to Modify

### CLI
- `chemsmart/cli/orca/qmmm.py`
  - Add `--geninput` and CrystalPrep flags
  - Gather CrystalPrep options in a separate config object
  - Enforce ionic-crystal + `.cif` gating

### ORCA QM/MM input
- `chemsmart/jobs/orca/settings.py`
  - No major changes to final `.inp` generation (keep stable)

### Submission script
- `chemsmart/settings/submitters.py` or `chemsmart/settings/server.py`
  - Extend script template to optionally insert a CrystalPrep/ORCA_mm section

### New module (small, isolated)
- `chemsmart/jobs/orca/crystalprep.py`
  - `CrystalPrepInputWriter` for template generation
  - Simple dataclass for CrystalPrep options

## Control Flow for `--geninput`

1. **Parse CLI**
   - Require `--jobtype ionic-crystal-qmmm`
   - Require `.cif` input (from `-f` or `--cp-input-cif`)
2. **Generate CrystalPrep template**
   - Write `{basename}.cp.inp` (or `--cp-template-out`)
3. **Generate submission script**
   - Insert prep steps before final ORCA run
   - Validate expected files after each step
4. **Final ORCA run**
   - Use existing QMMM writer to produce the final `.inp`
   - Run as normal once prep steps pass

## Shell Script Structure

The generated script should extend the existing `chemsmart_sub_*.sh` template with a prep block:

```bash
ORCA_crystalprep "$CP_TEMPLATE" -geninput
# check template exists
ORCA_crystalprep "$CP_TEMPLATE"
# check PDB / XYZ outputs
ORCA_mm -makeff "$SUPERCELL_XYZ" -CEL Na 1.0 -CEL Cl -1.0
# check ORCAFF file

# final ORCA run (existing behavior)
"$ORCA_DIR/orca" "$INPUT" > "$OUTPUT"
```

Add explicit `if [ ! -f ... ]` checks after each step.

## Backward Compatibility

- No behavior change for non-ionic-crystal jobs.
- No behavior change for ionic-crystal jobs unless `--geninput` is passed.
- The current final `.inp` generation remains the default path.

## Validation & Failure Modes

- Reject `--geninput` unless `--jobtype ionic-crystal-qmmm` is set.
- Reject if `.cif` input is missing.
- Validate `ORCA_crystalprep` and `ORCA_mm` outputs are present before ORCA run.
- Fail with explicit error messages on missing files.

## Testing Strategy

### Unit tests
- `tests/test_orca_cli.py`: flag gating + `.cif` validation
- `tests/test_orca_crystalprep.py`: template generation snapshot

### Script generation tests
- Add tests for the submission script builder to verify:
  - prep commands included
  - file checks present
  - final ORCA run still present

### Regression tests
- Keep existing ionic-crystal `.inp` tests intact
- Ensure non-ionic jobs unchanged

## Risks / Edge Cases

- ORCA_crystalprep file naming may vary by version; keep filenames overridable.
- Prep steps must run once per submission, not per node.
- ORCA_mm flags may evolve; keep command generation centralized in one helper.

## Next Steps

1. Add CrystalPrep template writer + data model
2. Extend QMMM CLI with `--geninput` and CrystalPrep flags
3. Add optional prep section to submission script generator
4. Add tests for CrystalPrep template and script output
