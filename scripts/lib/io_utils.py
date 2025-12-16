"""I/O utilities for MCP Rosetta scripts.

Common file operations and validation utilities.
"""

import os
import json
import logging
from pathlib import Path
from typing import Union, Dict, Any, Optional

logger = logging.getLogger(__name__)

def validate_pdb_file(file_path: Union[str, Path]) -> bool:
    """Validate that a PDB file exists and has basic PDB content.

    Args:
        file_path: Path to PDB file

    Returns:
        True if file is valid, False otherwise
    """
    file_path = Path(file_path)

    if not file_path.exists():
        logger.error(f"PDB file not found: {file_path}")
        return False

    if file_path.suffix.lower() not in ['.pdb', '.ent']:
        logger.warning(f"File does not have PDB extension: {file_path}")

    try:
        with open(file_path, 'r') as f:
            content = f.read(1000)  # Read first 1000 characters

        if 'ATOM' not in content and 'HETATM' not in content:
            logger.error(f"File does not contain ATOM records: {file_path}")
            return False

        return True
    except Exception as e:
        logger.error(f"Error reading PDB file {file_path}: {e}")
        return False

def load_config(config_path: Union[str, Path]) -> Dict[str, Any]:
    """Load configuration from JSON file.

    Args:
        config_path: Path to config file

    Returns:
        Configuration dictionary
    """
    config_path = Path(config_path)

    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")

    with open(config_path, 'r') as f:
        config = json.load(f)

    logger.info(f"Loaded configuration from: {config_path}")
    return config

def save_results(results: Dict[str, Any], output_path: Union[str, Path]) -> None:
    """Save results to JSON file.

    Args:
        results: Results dictionary
        output_path: Path to save results
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    logger.info(f"Results saved to: {output_path}")

def ensure_output_dir(output_file: Optional[Union[str, Path]]) -> Optional[Path]:
    """Ensure output directory exists.

    Args:
        output_file: Output file path

    Returns:
        Path object or None if no output file specified
    """
    if output_file is None:
        return None

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    return output_path