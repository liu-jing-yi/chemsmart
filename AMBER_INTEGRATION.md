# Amber Integration for chemsmart

## Overview

The Amber molecular dynamics package has been successfully integrated into chemsmart, providing comprehensive support for MD simulations, energy calculations, and geometry optimizations using the Amber force fields.

## Implementation Summary

### Core Components Implemented

1. **AmberJobSettings** (`chemsmart/jobs/amber/settings.py`)
   - Comprehensive settings class for Amber calculations
   - Supports MD, energy, and optimization simulations
   - Includes validation, merging, and configuration management
   - Compatible with existing chemsmart settings framework

2. **AmberJob Classes** (`chemsmart/jobs/amber/job.py`)
   - `AmberJob`: Base class for all Amber calculations
   - `AmberMDJob`: Specialized for molecular dynamics simulations
   - `AmberEnergyJob`: Specialized for single-point energy calculations
   - `AmberOptJob`: Specialized for geometry optimization
   - Full integration with chemsmart job framework

3. **AmberInputWriter** (`chemsmart/jobs/amber/writer.py`)
   - Generates proper Amber input files
   - Supports different simulation types (MD, energy, optimization)
   - Handles temperature and pressure control
   - Compatible with Amber input format

4. **AmberJobRunner** (`chemsmart/jobs/amber/runner.py`)
   - Specialized runner for Amber calculations
   - Supports both SANDER and PMEMD executables
   - Handles MPI and OpenMP parallelization
   - Environment setup and command execution

5. **CLI Integration** (`chemsmart/cli/amber/`)
   - Main `amber` command group
   - Subcommands: `md`, `energy`, `opt`
   - Full integration with chemsmart CLI framework
   - Support for all standard chemsmart options

6. **Project Settings** (`chemsmart/settings/amber.py`)
   - `AmberProjectSettings`: Base project configuration
   - `YamlAmberProjectSettings`: YAML-based configuration
   - Job-specific settings (MD, energy, optimization)
   - Compatible with existing project system

### Integration Points

1. **Main CLI Integration**
   - Added to `chemsmart/cli/subcommands.py`
   - Available via `chemsmart run amber ...` commands
   - Compatible with all chemsmart CLI features

2. **Job Framework Integration**
   - Inherits from `Job` base class
   - Compatible with existing jobrunner system
   - Supports skip_completed, local execution, etc.

3. **Settings Framework Integration**
   - Inherits from `MolecularJobSettings`
   - Implements required methods: `merge()`, `from_dict()`, `copy()`
   - Compatible with project-based configuration

## Usage Examples

### Basic Molecular Dynamics Simulation
```bash
chemsmart run amber -f system.pdb md -s 100000 -T 300 -e NVT
```

### NPT Ensemble Simulation
```bash
chemsmart run amber -f molecule.pdb md -s 50000 -T 310 -P 1.0 -e NPT
```

### Single-Point Energy Calculation
```bash
chemsmart run amber -f structure.pdb energy
```

### Geometry Optimization
```bash
chemsmart run amber -f initial.pdb opt -mc 1000
```

### With Custom Topology Files
```bash
chemsmart run amber md -top system.prmtop -crd system.inpcrd -s 200000
```

### Using Project Settings
```bash
chemsmart run amber -p my_project -f system.pdb md
```

### Submission to Cluster
```bash
chemsmart sub -s cluster_name amber -f system.pdb md -s 1000000
```

## Features

### Supported Ensembles
- **NVE**: Constant energy, volume, particles
- **NVT**: Constant temperature, volume, particles (Langevin thermostat)
- **NPT**: Constant temperature, pressure, particles (pressure coupling)

### Simulation Types
- **MD**: Molecular dynamics simulation
- **Energy**: Single-point energy calculation
- **Opt**: Geometry optimization using Amber force fields

### Advanced Features
- Pre-optimization of endpoint geometries
- Custom topology and coordinate file support
- Restart capability from previous calculations
- Integration with job scheduling systems
- Full compatibility with chemsmart molecule objects

### Configuration Options
- Temperature control (Langevin thermostat)
- Pressure control (isotropic pressure scaling)
- Timestep configuration
- Output frequency control
- Resource allocation (nodes, processors)
- Runtime limits and queue management

## File Structure

```
chemsmart/
├── jobs/amber/
│   ├── __init__.py          # Package exports
│   ├── job.py              # Job classes
│   ├── runner.py           # Job runner
│   ├── settings.py         # Settings classes
│   └── writer.py           # Input file writer
├── cli/amber/
│   ├── __init__.py         # CLI package exports
│   ├── amber.py            # Main amber command
│   ├── md.py               # MD subcommand
│   ├── energy.py           # Energy subcommand
│   └── opt.py              # Optimization subcommand
└── settings/
    └── amber.py            # Project settings
```

## Integration Status

✅ **Core Functionality**: Complete
✅ **Job Framework**: Complete  
✅ **CLI Integration**: Complete
✅ **Settings Framework**: Complete
✅ **Project Integration**: Complete
✅ **Documentation**: Complete

## Next Steps

The Amber integration is ready for use. Users can now:

1. Run Amber calculations using standard chemsmart commands
2. Use project-based configuration for Amber jobs
3. Submit Amber jobs to computing clusters
4. Leverage all existing chemsmart features with Amber

## Compatibility

The Amber integration maintains full compatibility with:
- Existing chemsmart job framework
- CLI command structure and options
- Project configuration system
- Job runners and scheduling
- Molecule handling and file I/O
- All other computational chemistry packages in chemsmart

The implementation follows the exact same patterns as the existing Gaussian and ORCA integrations, ensuring consistency and maintainability.
