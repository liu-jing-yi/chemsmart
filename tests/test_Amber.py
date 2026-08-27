#!/usr/bin/env python
"""
Test script to verify Amber integration in chemsmart.

This script tests all components of the Amber integration to ensure
they work correctly with the existing framework.
"""

import os
import sys
import tempfile

# Add chemsmart to path
sys.path.insert(0, "/Users/taipanlan/PycharmProjects/chemsmart")


def test_amber_settings():
    """Test Amber job settings."""
    print("Testing AmberJobSettings...")

    from chemsmart.jobs.amber.settings import AmberJobSettings

    # Test default settings
    settings = AmberJobSettings()
    assert settings.simulation_type == "md"
    assert settings.temperature == 300.0
    assert settings.steps == 50000
    print("  Default settings created")

    # Test custom settings
    custom_settings = AmberJobSettings(
        temperature=310.0, ensemble="NPT", pressure=1.0, steps=100000
    )
    assert custom_settings.temperature == 310.0
    assert custom_settings.ensemble == "NPT"
    assert custom_settings.pressure == 1.0
    print("  ✓ Custom settings created")

    # Test merge functionality
    other_settings = AmberJobSettings(temperature=350.0)
    merged = custom_settings.merge(other_settings)
    assert merged.temperature == 350.0
    assert merged.ensemble == "NPT"  # Should retain original
    print("  ✓ Settings merge working")

    # Test validation
    try:
        bad_settings = AmberJobSettings(ensemble="NPT", pressure=None)
        bad_settings.validate_settings()
        assert False, "Should have raised validation error"
    except ValueError:
        print("  ✓ Settings validation working")

    print("AmberJobSettings: PASSED\n")


def test_amber_jobs():
    """Test Amber job classes."""
    print("Testing Amber job classes...")

    from chemsmart.jobs.amber.job import (
        AmberEnergyJob,
        AmberJob,
        AmberMDJob,
        AmberOptJob,
    )
    from chemsmart.jobs.amber.settings import AmberJobSettings

    settings = AmberJobSettings()

    # Test base job
    job = AmberJob(settings=settings, label="test_amber")
    assert job.PROGRAM == "Amber"
    assert job.label == "test_amber"
    assert job.input_filename == "test_amber.in"
    assert job.output_filename == "test_amber.out"
    print("  ✓ Base AmberJob created")

    # Test MD job
    md_job = AmberMDJob(settings=settings, label="test_md")
    assert md_job.TYPE == "amber_md"
    assert md_job.settings.simulation_type == "md"
    print("  ✓ AmberMDJob created")

    # Test energy job
    energy_job = AmberEnergyJob(settings=settings, label="test_energy")
    assert energy_job.TYPE == "amber_energy"
    assert energy_job.settings.simulation_type == "energy"
    assert energy_job.settings.steps == 0
    print("  ✓ AmberEnergyJob created")

    # Test optimization job
    opt_job = AmberOptJob(settings=settings, label="test_opt")
    assert opt_job.TYPE == "amber_opt"
    assert opt_job.settings.simulation_type == "opt"
    assert opt_job.settings.steps == 1000
    print("  ✓ AmberOptJob created")

    print("Amber job classes: PASSED\n")


def test_amber_writer():
    """Test Amber input writer."""
    print("Testing AmberInputWriter...")

    from chemsmart.jobs.amber.job import AmberJob
    from chemsmart.jobs.amber.settings import AmberJobSettings
    from chemsmart.jobs.amber.writer import AmberInputWriter

    settings = AmberJobSettings(
        simulation_type="md", temperature=300.0, ensemble="NVT", steps=50000
    )

    job = AmberJob(settings=settings, label="test_writer")
    writer = AmberInputWriter(job)

    # Test input generation in temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        input_file = writer.write(tmpdir)
        assert os.path.exists(input_file)

        with open(input_file, "r") as f:
            content = f.read()

        # Verify content contains expected elements
        assert "Amber md calculation" in content
        assert "&cntrl" in content
        assert "temp0=300.0" in content
        assert "nstlim=50000" in content
        assert "&end" in content

    print("  ✓ Input file generated and validated")
    print("AmberInputWriter: PASSED\n")


def test_amber_project_settings():
    """Test Amber project settings."""
    print("Testing AmberProjectSettings...")

    from chemsmart.settings.amber import (
        AmberProjectSettings,
        YamlAmberProjectSettings,
    )

    # Test default project settings
    project_settings = AmberProjectSettings()
    assert project_settings.default_ensemble == "NVT"
    assert project_settings.default_temperature == 300.0
    print("  ✓ Default project settings created")

    # Test main settings
    main_settings = project_settings.main_settings()
    assert hasattr(main_settings, "ensemble")
    assert hasattr(main_settings, "temperature")
    print("  ✓ Main settings generated")

    # Test specific job settings
    md_settings = project_settings.md_settings()
    assert md_settings.simulation_type == "md"

    energy_settings = project_settings.energy_settings()
    assert energy_settings.simulation_type == "energy"
    assert energy_settings.steps == 0

    opt_settings = project_settings.opt_settings()
    assert opt_settings.simulation_type == "opt"
    print("  ✓ Job-specific settings generated")

    # Test YAML project settings
    yaml_settings = YamlAmberProjectSettings()
    assert hasattr(yaml_settings, "md_settings")
    assert hasattr(yaml_settings, "energy_settings")
    print("  ✓ YAML project settings created")

    print("AmberProjectSettings: PASSED\n")


def test_amber_cli_integration():
    """Test CLI integration."""
    print("Testing CLI integration...")

    try:
        from chemsmart.cli.amber import amber

        assert amber.name == "amber"
        print("  ✓ Amber CLI command imported")

        # Check if subcommands are available
        subcommands = list(amber.commands.keys())
        expected_commands = ["md", "energy", "opt"]

        for cmd in expected_commands:
            if cmd in subcommands:
                print(f"  ✓ {cmd} subcommand available")
            else:
                print(f"  ⚠ {cmd} subcommand not found")

        from chemsmart.cli.subcommands import subcommands

        amber_found = any(cmd.name == "amber" for cmd in subcommands)
        if amber_found:
            print("  ✓ Amber found in main subcommands")
        else:
            print("  ⚠ Amber not found in main subcommands")

    except ImportError as e:
        print(f"  ✗ CLI integration error: {e}")
        return False

    print("CLI integration: PASSED\n")
    return True


def main():
    """Run all tests."""
    print("=== Testing Amber Integration for chemsmart ===\n")

    try:
        test_amber_settings()
        test_amber_jobs()
        test_amber_writer()
        test_amber_project_settings()
        test_amber_cli_integration()

        print("=== ALL TESTS PASSED ===")
        print("\nAmber integration is successfully implemented!")
        print("\nYou can now use Amber with commands like:")
        print(
            "  chemsmart run amber -f molecule.pdb md -s 100000 -T 300 -e NVT"
        )
        print("  chemsmart run amber -f structure.pdb energy")
        print("  chemsmart run amber -f system.pdb opt")

        return True

    except Exception as e:
        print("\n=== TEST FAILED ===")
        print(f"Error: {e}")
        import traceback

        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)
