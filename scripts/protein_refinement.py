#!/usr/bin/env python3
"""
Script: protein_refinement.py
Description: Perform protein structure refinement using PyRosetta Relax protocol

Original Use Case: examples/use_case_1_protein_refinement.py
Dependencies Removed: None (all essential)
Repository Dependencies: PyRosetta (lazy loaded)

Usage:
    python scripts/protein_refinement.py --input <input_file> --output <output_file>

Example:
    python scripts/protein_refinement.py --input examples/data/sample.pdb --output results/refined.json
"""

# ==============================================================================
# Minimal Imports (only essential packages)
# ==============================================================================
import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Union, Optional, Dict, Any

# ==============================================================================
# Configuration (extracted from use case)
# ==============================================================================
DEFAULT_CONFIG = {
    "trajectories": 5,
    "cycles": 100,
    "temperature": 2.0,
    "pymol_visualization": False,
    "output_prefix": "refined"
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
            init, Pose, pose_from_file, get_fa_scorefxn, standard_packer_task,
            SequenceMover, MonteCarlo, TrialMover, RepeatMover, MoveMap,
            PyMOLMover, PyJobDistributor
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
            'SequenceMover': SequenceMover,
            'MonteCarlo': MonteCarlo,
            'TrialMover': TrialMover,
            'RepeatMover': RepeatMover,
            'MoveMap': MoveMap,
            'PyMOLMover': PyMOLMover,
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
    """Create minimal demo PDB structure."""
    demo_content = """ATOM      1  N   ALA A   1      20.154  11.200  19.962  1.00 20.00           N
ATOM      2  CA  ALA A   1      19.030  11.200  20.875  1.00 20.00           C
ATOM      3  C   ALA A   1      17.662  11.200  20.146  1.00 20.00           C
ATOM      4  O   ALA A   1      17.662  11.200  18.921  1.00 20.00           O
ATOM      5  CB  ALA A   1      19.030  12.424  21.788  1.00 20.00           C
TER
END
"""
    demo_path.parent.mkdir(parents=True, exist_ok=True)
    with open(demo_path, 'w') as f:
        f.write(demo_content)

# ==============================================================================
# Core Function (main logic extracted from use case)
# ==============================================================================
def run_protein_refinement(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform protein structure refinement using PyRosetta Relax protocol.

    Args:
        input_file: Path to input PDB file
        output_file: Path to save results (optional)
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters

    Returns:
        Dict containing:
            - result: Refinement results with scores and statistics
            - output_file: Path to output file (if saved)
            - metadata: Execution metadata

    Example:
        >>> result = run_protein_refinement("input.pdb", "output.json")
        >>> print(result['result']['best_score'])
    """
    logger = setup_logging()

    # Setup
    input_file = validate_input_file(input_file)
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    logger.info(f"Starting protein refinement for: {input_file}")
    logger.info(f"Configuration: {config}")

    try:
        # Load PyRosetta (lazy loading)
        pr = load_pyrosetta()
        logger.info("PyRosetta loaded successfully")

        # Load the input structure
        pose = pr['pose_from_file'](str(input_file))
        original_pose = pr['Pose']()
        original_pose.assign(pose)

        # Set up scoring function
        score_function = pr['get_fa_scorefxn']()
        original_score = score_function(pose)
        logger.info(f"Original score: {original_score:.2f}")

        # Set up move map (allow all backbone torsions to move)
        movemap = pr['MoveMap']()
        movemap.set_bb(True)
        movemap.set_chi(True)

        # Set up movers
        small_mover = pr['protocols'].simple_moves.SmallMover(movemap, 0.8, 1)
        shear_mover = pr['protocols'].simple_moves.ShearMover(movemap, 0.8, 1)

        # Minimization
        min_mover = pr['protocols'].simple_moves.MinMover()
        min_mover.movemap(movemap)
        min_mover.score_function(score_function)

        # Sidechain packing
        task_factory = pr['core'].pack.task.TaskFactory()
        packer = pr['protocols'].simple_moves.PackRotamersMover(score_function)
        packer.task_factory(task_factory)

        # Set up PyMOL mover if requested
        pymol_mover = None
        if config['pymol_visualization']:
            try:
                pymol_mover = pr['PyMOLMover']()
                pymol_mover.keep_history(True)
                logger.info("PyMOL mover initialized")
            except Exception:
                logger.warning("PyMOL mover failed to initialize, continuing without visualization")

        # Create sequence of moves
        sequence_mover = pr['SequenceMover']()
        sequence_mover.add_mover(small_mover)
        sequence_mover.add_mover(shear_mover)
        sequence_mover.add_mover(min_mover)
        sequence_mover.add_mover(packer)
        if pymol_mover:
            sequence_mover.add_mover(pymol_mover)

        # Set up Monte Carlo
        mc = pr['MonteCarlo'](pose, score_function, config['temperature'])

        # Create trial mover
        trial_mover = pr['TrialMover'](sequence_mover, mc)

        # Repeat mover
        repeat_mover = pr['RepeatMover'](trial_mover, config['cycles'])

        # Run refinement trajectories
        logger.info(f"Starting {config['trajectories']} trajectories with {config['cycles']} cycles each")
        scores = []

        for trajectory in range(config['trajectories']):
            # Reset pose for new trajectory
            pose.assign(original_pose)
            mc.reset(pose)

            # Set pose name
            pose.pdb_info().name(f"{config['output_prefix']}_{trajectory:04d}")

            if pymol_mover:
                pymol_mover.apply(pose)

            # Run refinement protocol
            repeat_mover.apply(pose)

            # Get final score
            final_score = score_function(pose)
            scores.append(final_score)

            logger.info(f"Trajectory {trajectory + 1}: Score = {final_score:.2f}")

        # Calculate statistics
        best_score = min(scores)
        improvement = original_score - best_score
        average_score = sum(scores) / len(scores)

        # Prepare results
        result = {
            "original_score": float(original_score),
            "best_score": float(best_score),
            "average_score": float(average_score),
            "score_improvement": float(improvement),
            "scores": [float(s) for s in scores],
            "trajectories_completed": len(scores)
        }

        logger.info("Refinement completed successfully!")
        logger.info(f"Original score: {original_score:.2f}")
        logger.info(f"Best refined score: {best_score:.2f}")
        logger.info(f"Score improvement: {improvement:.2f}")
        logger.info(f"Average score: {average_score:.2f}")

        # Save output if requested
        output_path = None
        if output_file:
            import json
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            output_data = {
                "refinement_results": result,
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
        logger.error(f"Error during protein refinement: {e}")
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
    parser.add_argument('--trajectories', '-n', type=int, help='Number of trajectories')
    parser.add_argument('--cycles', type=int, help='Monte Carlo cycles per trajectory')
    parser.add_argument('--temperature', type=float, help='Monte Carlo temperature')
    parser.add_argument('--pymol', action='store_true', help='Enable PyMOL visualization')
    parser.add_argument('--demo', action='store_true', help='Run with demo data')

    args = parser.parse_args()

    # Handle demo mode
    if args.demo:
        demo_pdb = Path('examples/data/test_input.pdb')
        if not demo_pdb.exists():
            logger = setup_logging()
            logger.warning(f"Demo file {demo_pdb} not found. Creating minimal demo structure...")
            create_demo_structure(demo_pdb)
        args.input = str(demo_pdb)
        if not args.trajectories:
            args.trajectories = 2  # Smaller for demo
        if not args.cycles:
            args.cycles = 10

    # Validate input file
    if not args.input:
        print("Error: Input file required! Use --input or --demo")
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
    if args.cycles:
        config_overrides['cycles'] = args.cycles
    if args.temperature:
        config_overrides['temperature'] = args.temperature
    if args.pymol:
        config_overrides['pymol_visualization'] = True

    # Run refinement
    result = run_protein_refinement(
        input_file=args.input,
        output_file=args.output,
        config=config,
        **config_overrides
    )

    if result['metadata']['success']:
        print(f"✅ Success: Refinement completed")
        if result['output_file']:
            print(f"📄 Results saved to: {result['output_file']}")
        if result['result']:
            print(f"🎯 Best score: {result['result']['best_score']:.2f}")
            print(f"📈 Score improvement: {result['result']['score_improvement']:.2f}")
        return 0
    else:
        print(f"❌ Error: {result['metadata']['error']}")
        return 1

if __name__ == '__main__':
    sys.exit(main())