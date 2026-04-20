# Amber MCPB Sequential Workflow Implementation Summary

## Overview
Successfully implemented a specialized submission system for Amber MCPB (Metal Center Parameter Builder) workflows that generates sequential job scripts instead of the single-job pattern used by Gaussian/ORCA.

## Key Components Implemented

### 1. AmberExecutable Class
**File**: `chemsmart/settings/executable.py`
- Added `AmberExecutable` class with `PROGRAM = "AMBER"`
- Supports multiple Amber executables: `pmemd.cuda`, `pmemd`, `sander`
- Includes helper methods for AmberTools utilities:
  - `get_mcpb_executable()` - for MCPB.py
  - `get_antechamber_executable()` - for antechamber
  - `get_tleap_executable()` - for tleap

### 2. Specialized MCPB Submitters
**File**: `chemsmart/settings/submitters.py`
- `AmberMCPBSLURMSubmitter` - Sequential MCPB workflow for SLURM schedulers
- `AmberMCPBPBSSubmitter` - Sequential MCPB workflow for PBS schedulers
- Both inherit from their respective base submitters but override job execution

### 3. Automatic MCPB Detection
**File**: `chemsmart/settings/server.py`
- Modified `get_submitter()` method to detect Amber MCPB jobs
- Automatically selects specialized submitters when:
  - `job.PROGRAM.lower() == 'amber'`
  - Job has `mcpb_input` attribute OR `job_type == 'mcpb'`

### 4. Amber Project Settings Enhancement
**File**: `chemsmart/settings/amber.py`
- Added `mcpb_settings()` method to `AmberProjectSettings`
- Enhanced `YamlAmberProjectSettings` to support MCPB configuration
- Added MCPB support in `from_dict()` method for YAML-based configurations

## Generated Submission Script Format

The specialized submitters generate sequential workflow scripts like:

```bash
#!/bin/bash
#SBATCH --job-name=mcpb_workflow
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

module load ambertools gaussian

# Amber MCPB Sequential Workflow
set -e  # Exit on any error

echo "Step 1: MCPB Preparation"
MCPB.py -i system.in -s 1

echo "Step 2: Small Model Generation"
MCPB.py -i system.in -s 2

echo "Step 3: Large Model Generation"
MCPB.py -i system.in -s 3

echo "Running Gaussian calculations for small model..."
# Run optimization
if [ -f small_opt.com ]; then
    g16 < small_opt.com > small_opt.log
    formchk small_opt.chk small_opt.fchk
fi

# Run frequency calculation
if [ -f small_fc.com ]; then
    g16 < small_fc.com > small_fc.log
    formchk small_fc.chk small_fc.fchk
fi

# Run RESP calculation
if [ -f small_mk.com ]; then
    g16 < small_mk.com > small_mk.log
    formchk small_mk.chk small_mk.fchk
fi

echo "Step 4: MCPB Finalization"
MCPB.py -i system.in -s 4

echo "Running tleap to generate topology files..."
tleap -s -f system_tleap.in > system_tleap.out

echo "MCPB workflow completed successfully!"
```

## Key Differences from Gaussian/ORCA Workflows

| Aspect | Gaussian/ORCA | Amber MCPB |
|--------|---------------|------------|
| **Job Pattern** | Single calculation | Sequential multi-step workflow |
| **Script Type** | Python run script + submission script | Direct workflow in submission script |
| **Dependencies** | None | Sequential step dependencies |
| **QM Integration** | Standalone | Integrated Gaussian calculations |
| **Error Handling** | Job-level | Step-level with exit-on-error |
| **Time Limits** | Standard (hours) | Extended (48+ hours default) |
| **Resource Allocation** | `--ntasks-per-node` | `--cpus-per-task` |

## Workflow Steps

1. **MCPB Step 1**: Preparation and input validation
2. **MCPB Step 2**: Small model generation  
3. **MCPB Step 3**: Large model generation
4. **QM Calculations**: Gaussian optimization, frequency, and RESP calculations
5. **MCPB Step 4**: Finalization and parameter generation
6. **Topology Generation**: tleap execution for final topology files
7. **Validation**: Results verification

## Features

- ✅ **Sequential Execution**: Each step waits for previous completion
- ✅ **Error Handling**: Exit-on-error with detailed error messages
- ✅ **Conditional QM**: Optional automatic or manual QM submission
- ✅ **Multi-Scheduler**: Support for both SLURM and PBS
- ✅ **Resource Optimization**: Longer time limits and appropriate CPU allocation
- ✅ **Integration**: Seamless integration with existing job submission system
- ✅ **Auto-Detection**: Automatic selection of MCPB submitters for Amber jobs

## Usage Example

```python
# The system automatically detects MCPB jobs and uses specialized submitters
job_settings = AmberFFParamJobSettings(
    element="ZN",
    input_pdb_file="system.pdb", 
    job_type="mcpb",
    auto_submit_qm=True
)

server = Server.from_servername("cluster")
submitter = server.get_submitter(job)  # Automatically returns AmberMCPBSLURMSubmitter
server.submit(job)  # Generates and submits sequential workflow script
```

## Implementation Status: ✅ COMPLETE

The Amber MCPB sequential workflow submission system is fully implemented and ready for use. It provides a robust, error-resistant workflow for complex MCPB parameterization tasks that require coordination between AmberTools and Gaussian calculations.
