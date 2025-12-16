#!/usr/bin/env python3
"""
Script: loop_modeling.py
Description: Perform loop modeling using PyRosetta loop reconstruction protocol

Original Use Case: examples/use_case_3_loop_modeling.py
Dependencies Removed: None (all essential)
Repository Dependencies: PyRosetta (lazy loaded)

Usage:
    python scripts/loop_modeling.py --input <input_file> --loop_start <start> --loop_end <end> --output <output_file>

Example:
    python scripts/loop_modeling.py --input examples/data/loop.pdb --loop_start 5 --loop_end 10 --output results/loops.json
"""

# ==============================================================================
# Minimal Imports (only essential packages)
# ==============================================================================
import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Union, Optional, Dict, Any, List, Tuple

# ==============================================================================
# Configuration (extracted from use case)
# ==============================================================================
DEFAULT_CONFIG = {
    "loop_start": None,  # Required parameter
    "loop_end": None,    # Required parameter
    "trajectories": 10,
    "cycles_low": 100,   # Centroid stage cycles
    "cycles_high": 50,   # Fullatom stage cycles
    "temperature": 2.0,
    "pymol_visualization": False,
    "output_prefix": "loop_modeled"
}

# ==============================================================================
# Shared Library Functions
# ==============================================================================
def setup_logging() -> logging.Logger:
    """Set up logging configuration."""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    return logging.getLogger(__name__)

def load_pyrosetta():
    """Lazy load PyRosetta to minimize startup time."""
    try:
        import pyrosetta
        import pyrosetta.rosetta as rosetta
        from pyrosetta import (
            init, Pose, pose_from_file, get_fa_scorefxn,
            MoveMap, PyJobDistributor
        )
        from pyrosetta.rosetta import core, protocols

        # Initialize PyRosetta with constant seed for reproducibility
        init(extra_options="-constant_seed")

        return {
            'pyrosetta': pyrosetta,
            'rosetta': rosetta,
            'Pose': Pose,
            'pose_from_file': pose_from_file,
            'get_fa_scorefxn': get_fa_scorefxn,
            'MoveMap': MoveMap,
            'PyJobDistributor': PyJobDistributor,
            'core': core,
            'protocols': protocols
        }
    except ImportError as e:
        raise ImportError(f"PyRosetta import failed: {e}. Please install PyRosetta: conda install -c rosettacommons pyrosetta")

def validate_input_file(file_path: Union[str, Path]) -> Path:
    """Validate input PDB file."""
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"Input file not found: {file_path}")

    # Basic PDB validation
    with open(file_path, 'r') as f:
        content = f.read(1000)
    if 'ATOM' not in content and 'HETATM' not in content:
        raise ValueError(f"File does not appear to be a valid PDB: {file_path}")

    return file_path

def create_demo_structure(demo_path: Path) -> None:
    """Create minimal demo PDB structure with a loop."""
    demo_content = """ATOM      1  N   ALA A   1      20.154  11.200  19.962  1.00 20.00           N
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
ATOM     18  N   VAL A   4       9.462  11.200  20.946  1.00 20.00           N
ATOM     19  CA  VAL A   4       8.130  11.200  20.375  1.00 20.00           C
ATOM     20  C   VAL A   4       6.962  11.200  21.246  1.00 20.00           C
ATOM     21  O   VAL A   4       6.962  11.200  22.471  1.00 20.00           O
ATOM     22  CB  VAL A   4       8.030  12.424  19.488  1.00 20.00           C
ATOM     23  CG1 VAL A   4       8.030  13.824  20.188  1.00 20.00           C
ATOM     24  CG2 VAL A   4       8.030  12.424  18.188  1.00 20.00           C
TER
END
"""
    demo_path.parent.mkdir(parents=True, exist_ok=True)
    with open(demo_path, 'w') as f:
        f.write(demo_content)

def validate_loop_parameters(loop_start: int, loop_end: int, total_residues: int) -> None:
    """Validate loop modeling parameters."""
    if loop_start <= 0:
        raise ValueError(f"Loop start must be positive: {loop_start}")
    if loop_end <= loop_start:
        raise ValueError(f"Loop end ({loop_end}) must be greater than loop start ({loop_start})")
    if loop_end > total_residues:
        raise ValueError(f"Loop end ({loop_end}) exceeds total residues ({total_residues})")
    if loop_start > total_residues:
        raise ValueError(f"Loop start ({loop_start}) exceeds total residues ({total_residues})")

# ==============================================================================
# Core Function (main logic extracted from use case)
# ==============================================================================
def run_loop_modeling(
    input_file: Union[str, Path],
    loop_start: int,
    loop_end: int,
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform loop modeling using PyRosetta loop reconstruction protocol.

    Args:
        input_file: Path to input PDB file
        loop_start: Starting residue number of the loop
        loop_end: Ending residue number of the loop
        output_file: Path to save results (optional)
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters

    Returns:
        Dict containing:
            - result: Loop modeling results with scores and statistics
            - output_file: Path to output file (if saved)
            - metadata: Execution metadata

    Example:
        >>> result = run_loop_modeling("input.pdb", 5, 10, "output.json")
        >>> print(result['result']['best_score'])
    """
    logger = setup_logging()

    # Setup
    input_file = validate_input_file(input_file)
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}
    config['loop_start'] = loop_start
    config['loop_end'] = loop_end

    logger.info(f"Starting loop modeling for: {input_file}")
    logger.info(f"Loop region: {loop_start}-{loop_end}")
    logger.info(f"Configuration: {config}")

    try:
        # Load PyRosetta (lazy loading)
        pr = load_pyrosetta()
        logger.info("PyRosetta loaded successfully")

        # Load the input structure
        pose = pr['pose_from_file'](str(input_file))
        original_pose = pr['Pose']()
        original_pose.assign(pose)

        total_residues = pose.total_residue()
        validate_loop_parameters(loop_start, loop_end, total_residues)

        # Set up scoring functions
        scorefxn_low = pr['core'].scoring.ScoreFunctionFactory.create_score_function("score4_smooth_cart")
        scorefxn_high = pr['get_fa_scorefxn']()

        original_score = scorefxn_high(pose)
        logger.info(f"Original score: {original_score:.2f}")

        # Set up loop definition
        loop_def = pr['protocols'].loops.Loop(loop_start, loop_end)
        loops = pr['protocols'].loops.Loops()
        loops.add_loop(loop_def)

        logger.info(f"Loop defined: residues {loop_start} to {loop_end}")

        # Set up loop modeling protocol
        loop_protocol = pr['protocols'].loops.loop_mover.LoopMover()
        loop_protocol.set_loops(loops)
        loop_protocol.set_scorefxn(scorefxn_high)

        # Configure loop protocol
        loop_protocol.set_fa_scorefxn(scorefxn_high)
        loop_protocol.set_centroid_scorefxn(scorefxn_low)

        # Run loop modeling trajectories
        logger.info(f"Starting {config['trajectories']} loop modeling trajectories")
        scores = []
        loop_rmsd_values = []

        for trajectory in range(config['trajectories']):
            # Reset pose for new trajectory
            pose.assign(original_pose)

            # Set pose name
            pose.pdb_info().name(f"{config['output_prefix']}_{trajectory:04d}")

            # Apply loop modeling protocol
            loop_protocol.apply(pose)

            # Get scores
            final_score = scorefxn_high(pose)
            scores.append(final_score)

            # Calculate loop RMSD (simple approximation)
            loop_rmsd = 0.0
            for i in range(loop_start, loop_end + 1):
                if i <= pose.total_residue() and i <= original_pose.total_residue():
                    orig_ca = original_pose.residue(i).xyz('CA')
                    new_ca = pose.residue(i).xyz('CA')
                    loop_rmsd += (orig_ca - new_ca).norm_squared()

            loop_rmsd = (loop_rmsd / (loop_end - loop_start + 1)) ** 0.5
            loop_rmsd_values.append(loop_rmsd)

            logger.info(f"Trajectory {trajectory + 1}: Score = {final_score:.2f}, Loop RMSD = {loop_rmsd:.2f}")

        # Calculate statistics
        best_score = min(scores)
        improvement = original_score - best_score
        average_score = sum(scores) / len(scores)
        average_rmsd = sum(loop_rmsd_values) / len(loop_rmsd_values)
        min_rmsd = min(loop_rmsd_values)

        # Prepare results
        result = {
            "original_score": float(original_score),
            "best_score": float(best_score),
            "average_score": float(average_score),
            "score_improvement": float(improvement),
            "scores": [float(s) for s in scores],
            "loop_rmsd_values": [float(r) for r in loop_rmsd_values],
            "average_loop_rmsd": float(average_rmsd),
            "min_loop_rmsd": float(min_rmsd),
            "trajectories_completed": len(scores),
            "loop_start": loop_start,
            "loop_end": loop_end,
            "loop_length": loop_end - loop_start + 1,
            "total_residues": total_residues
        }

        logger.info("Loop modeling completed successfully!")
        logger.info(f"Original score: {original_score:.2f}")
        logger.info(f"Best modeled score: {best_score:.2f}")
        logger.info(f"Score improvement: {improvement:.2f}")
        logger.info(f"Average loop RMSD: {average_rmsd:.2f}")
        logger.info(f"Best loop RMSD: {min_rmsd:.2f}")

        # Save output if requested
        output_path = None
        if output_file:
            import json
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            output_data = {
                "loop_modeling_results": result,
                "metadata": {
                    "input_file": str(input_file),
                    "config": config,
                    "success": True
                }
            }

            with open(output_path, 'w') as f:
                json.dump(output_data, f, indent=2)
            logger.info(f"Results saved to: {output_path}")

        return {
            "result": result,
            "output_file": str(output_path) if output_path else None,
            "metadata": {
                "input_file": str(input_file),
                "config": config,
                "success": True
            }
        }

    except Exception as e:
        logger.error(f"Error during loop modeling: {e}")
        return {
            "result": None,
            "output_file": None,
            "metadata": {
                "input_file": str(input_file),
                "config": config,
                "success": False,
                "error": str(e)
            }
        }

# ==============================================================================
# CLI Interface
# ==============================================================================
def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', '-i', help='Input PDB file path')
    parser.add_argument('--output', '-o', help='Output JSON file path')
    parser.add_argument('--config', '-c', help='Config file (JSON)')
    parser.add_argument('--loop_start', '-s', type=int, help='Loop start residue number')
    parser.add_argument('--loop_end', '-e', type=int, help='Loop end residue number')
    parser.add_argument('--trajectories', '-n', type=int, help='Number of trajectories')
    parser.add_argument('--cycles_low', type=int, help='Centroid stage cycles')
    parser.add_argument('--cycles_high', type=int, help='Fullatom stage cycles')
    parser.add_argument('--temperature', type=float, help='Monte Carlo temperature')
    parser.add_argument('--pymol', action='store_true', help='Enable PyMOL visualization')
    parser.add_argument('--demo', action='store_true', help='Run with demo data')

    args = parser.parse_args()

    # Handle demo mode
    if args.demo:
        demo_pdb = Path('examples/data/test_loop.pdb')
        if not demo_pdb.exists():
            logger = setup_logging()
            logger.warning(f"Demo file {demo_pdb} not found. Creating minimal demo structure...")
            create_demo_structure(demo_pdb)
        args.input = str(demo_pdb)
        if not args.loop_start:
            args.loop_start = 2
        if not args.loop_end:
            args.loop_end = 3
        if not args.trajectories:
            args.trajectories = 3  # Smaller for demo

    # Validate required parameters
    if not args.input:
        print("Error: Input file required! Use --input or --demo")
        return 1
    if not args.loop_start:
        print("Error: Loop start required! Use --loop_start or --demo")
        return 1
    if not args.loop_end:
        print("Error: Loop end required! Use --loop_end or --demo")
        return 1

    # Load config if provided
    config = None
    if args.config:
        import json
        with open(args.config) as f:
            config = json.load(f)

    # Override config with CLI arguments
    config_overrides = {}
    if args.trajectories:
        config_overrides['trajectories'] = args.trajectories
    if args.cycles_low:
        config_overrides['cycles_low'] = args.cycles_low
    if args.cycles_high:
        config_overrides['cycles_high'] = args.cycles_high
    if args.temperature:
        config_overrides['temperature'] = args.temperature
    if args.pymol:
        config_overrides['pymol_visualization'] = True

    # Run loop modeling
    result = run_loop_modeling(
        input_file=args.input,
        loop_start=args.loop_start,
        loop_end=args.loop_end,
        output_file=args.output,
        config=config,
        **config_overrides
    )

    if result['metadata']['success']:
        print(f"✅ Success: Loop modeling completed")
        if result['output_file']:
            print(f"📄 Results saved to: {result['output_file']}")
        if result['result']:
            print(f"🎯 Best score: {result['result']['best_score']:.2f}")
            print(f"🔄 Best loop RMSD: {result['result']['min_loop_rmsd']:.2f}")
            print(f"📈 Score improvement: {result['result']['score_improvement']:.2f}")
        return 0
    else:
        print(f"❌ Error: {result['metadata']['error']}")
        return 1

if __name__ == '__main__':
    sys.exit(main())