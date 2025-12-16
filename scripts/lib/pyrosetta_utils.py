"""PyRosetta utilities for MCP Rosetta scripts.

This module provides a common interface for PyRosetta functionality
with graceful handling of the dependency.
"""

import logging
from typing import Tuple, Optional, Any

logger = logging.getLogger(__name__)

def setup_pyrosetta() -> Tuple[bool, Optional[Any], Optional[Any]]:
    """Initialize PyRosetta with standard settings.

    Returns:
        Tuple of (success, pyrosetta_module, rosetta_module)
        If PyRosetta is not available, returns (False, None, None)
    """
    try:
        import pyrosetta
        import pyrosetta.rosetta as rosetta

        # Initialize PyRosetta with constant seed for reproducibility
        pyrosetta.init(extra_options="-constant_seed")
        logger.info("PyRosetta initialized successfully")
        return True, pyrosetta, rosetta
    except ImportError as e:
        logger.error(f"PyRosetta import failed: {e}")
        logger.error("Please install PyRosetta: conda install -c rosettacommons pyrosetta")
        return False, None, None

def create_demo_pdb(file_path: str, structure_type: str = "single") -> None:
    """Create a minimal demo PDB structure.

    Args:
        file_path: Path where to save the demo PDB
        structure_type: Type of structure to create ("single", "complex", "loop")
    """
    import os

    if structure_type == "single":
        # Simple alanine residue
        content = """ATOM      1  N   ALA A   1      20.154  11.200  19.962  1.00 20.00           N
ATOM      2  CA  ALA A   1      19.030  11.200  20.875  1.00 20.00           C
ATOM      3  C   ALA A   1      17.662  11.200  20.146  1.00 20.00           C
ATOM      4  O   ALA A   1      17.662  11.200  18.921  1.00 20.00           O
ATOM      5  CB  ALA A   1      19.030  12.424  21.788  1.00 20.00           C
TER
END
"""
    elif structure_type == "complex":
        # Two-chain complex
        content = """ATOM      1  N   ALA A   1      20.154  11.200  19.962  1.00 20.00           N
ATOM      2  CA  ALA A   1      19.030  11.200  20.875  1.00 20.00           C
ATOM      3  C   ALA A   1      17.662  11.200  20.146  1.00 20.00           C
ATOM      4  O   ALA A   1      17.662  11.200  18.921  1.00 20.00           O
ATOM      5  CB  ALA A   1      19.030  12.424  21.788  1.00 20.00           C
TER
ATOM      6  N   GLY B   1      25.154  15.200  25.962  1.00 20.00           N
ATOM      7  CA  GLY B   1      24.030  15.200  26.875  1.00 20.00           C
ATOM      8  C   GLY B   1      22.662  15.200  26.146  1.00 20.00           C
ATOM      9  O   GLY B   1      22.662  15.200  24.921  1.00 20.00           O
TER
END
"""
    elif structure_type == "loop":
        # Structure with a loop region
        content = """ATOM      1  N   ALA A   1      20.154  11.200  19.962  1.00 20.00           N
ATOM      2  CA  ALA A   1      19.030  11.200  20.875  1.00 20.00           C
ATOM      3  C   ALA A   1      17.662  11.200  20.146  1.00 20.00           C
ATOM      4  O   ALA A   1      17.662  11.200  18.921  1.00 20.00           O
ATOM      5  CB  ALA A   1      19.030  12.424  21.788  1.00 20.00           C
ATOM      6  N   GLY A   2      16.662  11.200  20.846  1.00 20.00           N
ATOM      7  CA  GLY A   2      15.330  11.200  20.275  1.00 20.00           C
ATOM      8  C   GLY A   2      14.162  11.200  21.146  1.00 20.00           C
ATOM      9  O   GLY A   2      14.162  11.200  22.371  1.00 20.00           O
ATOM     10  N   LEU A   3      13.062  11.200  20.546  1.00 20.00           N
ATOM     11  CA  LEU A   3      11.730  11.200  21.175  1.00 20.00           C
ATOM     12  C   LEU A   3      10.562  11.200  20.346  1.00 20.00           C
ATOM     13  O   LEU A   3      10.562  11.200  19.121  1.00 20.00           O
ATOM     14  CB  LEU A   3      11.630  12.424  22.088  1.00 20.00           C
ATOM     15  CG  LEU A   3      11.630  13.824  21.488  1.00 20.00           C
ATOM     16  CD1 LEU A   3      11.630  13.824  19.988  1.00 20.00           C
ATOM     17  CD2 LEU A   3      11.630  15.024  22.188  1.00 20.00           C
TER
END
"""
    else:
        raise ValueError(f"Unknown structure_type: {structure_type}")

    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, 'w') as f:
        f.write(content)
    logger.info(f"Created demo PDB: {file_path}")