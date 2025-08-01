import subprocess
import os
import filecmp
import pytest

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(TEST_DIR)
SCRIPT = os.path.join(PROJECT_DIR, "extract_relationships_from_ped.py")
INPUT_PED = os.path.join(TEST_DIR, "example_no_consanguinity.ped")
EXPECTED_OUTPUT = os.path.join(TEST_DIR, "example_no_consanguinity_with_relationships.ped")
ACTUAL_OUTPUT = os.path.join(TEST_DIR, "test_output.ped")

def test_relationship_extraction():
    # Run the script
    result = subprocess.run([
        "python", SCRIPT, INPUT_PED, ACTUAL_OUTPUT
    ], capture_output=True, text=True)
    assert result.returncode == 0, f"Script failed: {result.stderr}"
    # Compare output files (ignoring comment/header lines)
    with open(EXPECTED_OUTPUT, 'r') as f:
        expected_lines = [line for line in f if not line.startswith('#') and line.strip()]
    with open(ACTUAL_OUTPUT, 'r') as f:
        actual_lines = [line for line in f if not line.startswith('#') and line.strip()]
    assert expected_lines == actual_lines, "Output does not match expected relationships."
