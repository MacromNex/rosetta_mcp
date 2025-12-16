#!/usr/bin/env python3
"""
Script: protein_docking.py
Description: Perform protein-protein docking using PyRosetta DockingProtocol

Original Use Case: examples/use_case_2_protein_docking.py
Dependencies Removed: None (all essential)
Repository Dependencies: PyRosetta (lazy loaded)

Usage:
    python scripts/protein_docking.py --input <input_file> --output <output_file>

Example:
    python scripts/protein_docking.py --input examples/data/complex.pdb --output results/docked.json
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
    "trajectories": 10,
    "chain_break": None,  # Auto-detect if None
    "use_centroid_stage": True,
    "use_fullatom_stage": True,
    "pymol_visualization": False,
    "output_prefix": "docked"
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
            FoldTree, MoveMap, PyJobDistributor
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
            'FoldTree': FoldTree,
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
    """Create minimal demo PDB structure with two chains."""
    demo_content = """ATOM      1  N   ALA A   1      20.154  11.200  19.962  1.00 20.00           N
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
    demo_path.parent.mkdir(parents=True, exist_ok=True)
    with open(demo_path, 'w') as f:
        f.write(demo_content)

def setup_docking_foldtree(pose, chain_break, pr):
    """Set up fold tree for docking."""
    ft = pr['FoldTree']()

    total_res = pose.total_residue()
    ft.add_edge(1, chain_break, -1)  # Chain A backbone
    ft.add_edge(chain_break + 1, total_res, -1)  # Chain B backbone
    ft.add_edge(chain_break, chain_break + 1, 1)  # Jump between chains

    pose.fold_tree(ft)
    return ft

# ==============================================================================
# Core Function (main logic extracted from use case)
# ==============================================================================
def run_protein_docking(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform protein-protein docking using PyRosetta DockingProtocol.

    Args:
        input_file: Path to input PDB file with both protein chains
        output_file: Path to save results (optional)
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters

    Returns:
        Dict containing:
            - result: Docking results with scores and interface analysis
            - output_file: Path to output file (if saved)
            - metadata: Execution metadata

    Example:
        >>> result = run_protein_docking("complex.pdb", "docked.json")
        >>> print(result['result']['best_score'])
    """
    logger = setup_logging()

    # Setup
    input_file = validate_input_file(input_file)
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    logger.info(f"Starting protein-protein docking for: {input_file}")
    logger.info(f"Configuration: {config}")

    try:
        # Load PyRosetta (lazy loading)
        pr = load_pyrosetta()
        logger.info("PyRosetta loaded successfully")

        # Load the input structure
        pose = pr['pose_from_file'](str(input_file))
        original_pose = pr['Pose']()
        original_pose.assign(pose)

        # Determine chain break if not provided
        chain_break = config['chain_break']
        if chain_break is None:
            # Simple heuristic: find middle or use chain information
            chain_break = pose.total_residue() // 2
            logger.info(f"Auto-detected chain break at residue {chain_break}")

        # Set up docking fold tree
        setup_docking_foldtree(pose, chain_break, pr)
        logger.info(f"Docking fold tree set up: Chain A (1-{chain_break}), Chain B ({chain_break + 1}-{pose.total_residue()})")

        # Set up scoring functions
        scorefxn_low = pr['core'].scoring.ScoreFunctionFactory.create_score_function("interchain_cen")
        scorefxn_high = pr['get_fa_scorefxn']()

        # Adjust interface weights
        scorefxn_high.set_weight(pr['core'].scoring.fa_rep, 0.1)
        scorefxn_high.set_weight(pr['core'].scoring.fa_atr, 0.8)

        original_score = scorefxn_high(pose)
        logger.info(f"Original score: {original_score:.2f}")

        # Set up docking protocol
        dock_protocol = pr['protocols'].docking.DockingProtocol()
        dock_protocol.set_scorefxn(scorefxn_high)
        dock_protocol.set_scorefxn_low(scorefxn_low)

        # Run docking trajectories
        logger.info(f"Starting {config['trajectories']} docking trajectories")
        scores = []
        interface_scores = []

        for trajectory in range(config['trajectories']):
            # Reset pose for new trajectory
            pose.assign(original_pose)

            # Set pose name
            pose.pdb_info().name(f"{config['output_prefix']}_{trajectory:04d}")

            # Apply random perturbation
            perturb = pr['protocols'].rigid.RigidBodyPerturbMover(1, 3.0, 8.0)
            perturb.apply(pose)

            # Run docking protocol
            dock_protocol.apply(pose)

            # Get scores
            final_score = scorefxn_high(pose)

            # Calculate interface score
            interface_analyzer = pr['protocols'].analysis.InterfaceAnalyzerMover()
            interface_analyzer.apply(pose)
            interface_score = pose.scores[pr['core'].scoring.interface_delta]

            scores.append(final_score)
            interface_scores.append(interface_score)

            logger.info(f"Trajectory {trajectory + 1}: Score = {final_score:.2f}, Interface = {interface_score:.2f}")

        # Calculate statistics
        best_score = min(scores)
        best_interface = min(interface_scores)
        improvement = original_score - best_score
        average_score = sum(scores) / len(scores)
        average_interface = sum(interface_scores) / len(interface_scores)

        # Prepare results
        result = {
            "original_score": float(original_score),
            "best_score": float(best_score),
            "best_interface_score": float(best_interface),
            "average_score": float(average_score),
            "average_interface_score": float(average_interface),
            "score_improvement": float(improvement),
            "scores": [float(s) for s in scores],
            "interface_scores": [float(s) for s in interface_scores],
            "trajectories_completed": len(scores),
            "chain_break": chain_break,
            "total_residues": pose.total_residue()
        }

        logger.info("Docking completed successfully!")
        logger.info(f"Original score: {original_score:.2f}")
        logger.info(f"Best docked score: {best_score:.2f}")
        logger.info(f"Best interface score: {best_interface:.2f}")
        logger.info(f"Score improvement: {improvement:.2f}")

        # Save output if requested
        output_path = None
        if output_file:
            import json
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            output_data = {
                "docking_results": result,
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
        logger.error(f"Error during protein docking: {e}")
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
    parser.add_argument('--chain_break', type=int, help='Chain break residue number')
    parser.add_argument('--pymol', action='store_true', help='Enable PyMOL visualization')
    parser.add_argument('--demo', action='store_true', help='Run with demo data')

    args = parser.parse_args()

    # Handle demo mode
    if args.demo:
        demo_pdb = Path('examples/data/test_docking.pdb')
        if not demo_pdb.exists():
            logger = setup_logging()
            logger.warning(f"Demo file {demo_pdb} not found. Creating minimal demo structure...")
            create_demo_structure(demo_pdb)
        args.input = str(demo_pdb)
        if not args.trajectories:
            args.trajectories = 3  # Smaller for demo
        if not args.chain_break:
            args.chain_break = 5  # Demo chain break

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
    if args.chain_break:
        config_overrides['chain_break'] = args.chain_break
    if args.pymol:
        config_overrides['pymol_visualization'] = True

    # Run docking
    result = run_protein_docking(
        input_file=args.input,
        output_file=args.output,
        config=config,
        **config_overrides
    )

    if result['metadata']['success']:
        print(f"✅ Success: Docking completed")
        if result['output_file']:
            print(f"📄 Results saved to: {result['output_file']}")
        if result['result']:
            print(f"🎯 Best score: {result['result']['best_score']:.2f}")
            print(f"🔗 Best interface: {result['result']['best_interface_score']:.2f}")
            print(f"📈 Score improvement: {result['result']['score_improvement']:.2f}")
        return 0
    else:
        print(f"❌ Error: {result['metadata']['error']}")
        return 1

if __name__ == '__main__':
    sys.exit(main())