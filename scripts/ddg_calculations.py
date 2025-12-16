#!/usr/bin/env python3
"""
Script: ddg_calculations.py
Description: Perform ddG (ΔΔG) calculations using PyRosetta for stability analysis

Original Use Case: examples/use_case_5_ddg_calculations.py
Dependencies Removed: None (all essential)
Repository Dependencies: PyRosetta (lazy loaded)

Usage:
    python scripts/ddg_calculations.py --input <input_file> --mutations <mutations> --output <output_file>

Example:
    python scripts/ddg_calculations.py --input examples/data/protein.pdb --mutations "A10G,B15L" --output results/ddg.json
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
    "mutations": [],  # List of mutations in format "A10G" (chain, position, new_aa)
    "trajectories": 5,
    "repack_radius": 8.0,  # Radius around mutation for repacking
    "minimize": True,
    "use_cartesian": False,
    "temperature": 2.0,
    "output_prefix": "ddg"
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
    """Create minimal demo PDB structure."""
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
TER
END
"""
    demo_path.parent.mkdir(parents=True, exist_ok=True)
    with open(demo_path, 'w') as f:
        f.write(demo_content)

def parse_mutations(mutation_string: str) -> List[Tuple[str, int, str]]:
    """Parse mutation string into list of (chain, position, new_aa) tuples.

    Args:
        mutation_string: String like "A10G,B15L" or "10G,15L" (assumes chain A)

    Returns:
        List of (chain, position, new_aa) tuples
    """
    if not mutation_string:
        return []

    mutations = []
    for mut in mutation_string.split(','):
        mut = mut.strip()
        if len(mut) < 2:
            continue

        # Check if chain is specified
        if mut[0].isalpha() and len(mut) >= 3:
            # Format: A10G (chain specified)
            chain = mut[0]
            position = int(mut[1:-1])
            new_aa = mut[-1]
        else:
            # Format: 10G (assumes chain A)
            chain = 'A'
            position = int(mut[:-1])
            new_aa = mut[-1]

        mutations.append((chain, position, new_aa))

    return mutations

def apply_mutation(pose, chain: str, position: int, new_aa: str, pr):
    """Apply a single point mutation to the pose."""
    # Find the residue number in the pose
    residue_num = None
    for i in range(1, pose.total_residue() + 1):
        if (pose.pdb_info().chain(i) == chain and
            pose.pdb_info().number(i) == position):
            residue_num = i
            break

    if residue_num is None:
        raise ValueError(f"Could not find residue {chain}{position}")

    # Create mutator
    mutator = pr['protocols'].simple_moves.MutateResidue()
    mutator.set_target(residue_num)
    mutator.set_res_name(new_aa)

    # Apply mutation
    mutator.apply(pose)

    return residue_num

# ==============================================================================
# Core Function (main logic extracted from use case)
# ==============================================================================
def run_ddg_calculations(
    input_file: Union[str, Path],
    mutations: Union[str, List[Tuple[str, int, str]]],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]:
    """
    Perform ddG (ΔΔG) calculations using PyRosetta for stability analysis.

    Args:
        input_file: Path to input PDB file
        mutations: Either mutation string ("A10G,B15L") or list of (chain, pos, aa) tuples
        output_file: Path to save results (optional)
        config: Configuration dict (uses DEFAULT_CONFIG if not provided)
        **kwargs: Override specific config parameters

    Returns:
        Dict containing:
            - result: ddG calculation results for each mutation
            - output_file: Path to output file (if saved)
            - metadata: Execution metadata

    Example:
        >>> result = run_ddg_calculations("input.pdb", "A10G", "ddg.json")
        >>> print(result['result']['A10G']['ddg'])
    """
    logger = setup_logging()

    # Setup
    input_file = validate_input_file(input_file)
    config = {**DEFAULT_CONFIG, **(config or {}), **kwargs}

    # Parse mutations
    if isinstance(mutations, str):
        mutation_list = parse_mutations(mutations)
    else:
        mutation_list = mutations

    if not mutation_list:
        raise ValueError("No valid mutations provided")

    logger.info(f"Starting ddG calculations for: {input_file}")
    logger.info(f"Mutations: {mutation_list}")
    logger.info(f"Configuration: {config}")

    try:
        # Load PyRosetta (lazy loading)
        pr = load_pyrosetta()
        logger.info("PyRosetta loaded successfully")

        # Load the input structure
        wild_type_pose = pr['pose_from_file'](str(input_file))

        # Set up scoring function
        scorefxn = pr['get_fa_scorefxn']()

        # Calculate wild-type score
        wild_type_score = scorefxn(wild_type_pose)
        logger.info(f"Wild-type score: {wild_type_score:.2f}")

        # Results dictionary
        mutation_results = {}

        # Process each mutation
        for chain, position, new_aa in mutation_list:
            mutation_name = f"{chain}{position}{new_aa}"
            logger.info(f"Processing mutation: {mutation_name}")

            # Create mutant pose
            mutant_pose = pr['Pose']()
            mutant_pose.assign(wild_type_pose)

            try:
                # Apply mutation
                residue_num = apply_mutation(mutant_pose, chain, position, new_aa, pr)

                # Set up repacking around mutation
                if config['minimize'] or config['repack_radius'] > 0:
                    # Create task factory for repacking
                    task_factory = pr['core'].pack.task.TaskFactory()

                    # Restrict to sphere around mutation
                    if config['repack_radius'] > 0:
                        sphere_selector = pr['core'].select.residue_selector.NeighborhoodResidueSelector()
                        sphere_selector.set_focus_selector(
                            pr['core'].select.residue_selector.ResidueIndexSelector(str(residue_num))
                        )
                        sphere_selector.set_distance(config['repack_radius'])

                        # Only repack selected residues
                        restrict_to_sphere = pr['core'].pack.task.operation.RestrictToRepackingRLT()
                        task_factory.push_back(restrict_to_sphere)

                    # Pack rotamers
                    packer = pr['protocols'].simple_moves.PackRotamersMover(scorefxn)
                    packer.task_factory(task_factory)
                    packer.apply(mutant_pose)

                # Minimize if requested
                if config['minimize']:
                    movemap = pr['MoveMap']()

                    # Allow movement around mutation site
                    for i in range(max(1, residue_num - 5), min(mutant_pose.total_residue() + 1, residue_num + 6)):
                        movemap.set_chi(i, True)

                    min_mover = pr['protocols'].simple_moves.MinMover()
                    min_mover.movemap(movemap)
                    min_mover.score_function(scorefxn)
                    min_mover.apply(mutant_pose)

                # Calculate mutant score
                mutant_score = scorefxn(mutant_pose)

                # Calculate ddG
                ddg = mutant_score - wild_type_score

                mutation_results[mutation_name] = {
                    "chain": chain,
                    "position": position,
                    "new_aa": new_aa,
                    "residue_number": residue_num,
                    "wild_type_score": float(wild_type_score),
                    "mutant_score": float(mutant_score),
                    "ddg": float(ddg),
                    "destabilizing": ddg > 0.0,
                    "success": True
                }

                logger.info(f"  {mutation_name}: ΔΔG = {ddg:.2f} kcal/mol ({'destabilizing' if ddg > 0 else 'stabilizing'})")

            except Exception as e:
                logger.error(f"  {mutation_name}: Failed - {e}")
                mutation_results[mutation_name] = {
                    "chain": chain,
                    "position": position,
                    "new_aa": new_aa,
                    "error": str(e),
                    "success": False
                }

        # Calculate summary statistics
        successful_mutations = [r for r in mutation_results.values() if r['success']]
        if successful_mutations:
            ddg_values = [r['ddg'] for r in successful_mutations]
            average_ddg = sum(ddg_values) / len(ddg_values)
            destabilizing_count = sum(1 for r in successful_mutations if r['destabilizing'])
        else:
            average_ddg = 0.0
            destabilizing_count = 0

        # Prepare results
        result = {
            "wild_type_score": float(wild_type_score),
            "mutations": mutation_results,
            "summary": {
                "total_mutations": len(mutation_list),
                "successful_mutations": len(successful_mutations),
                "failed_mutations": len(mutation_list) - len(successful_mutations),
                "average_ddg": float(average_ddg),
                "destabilizing_mutations": destabilizing_count,
                "stabilizing_mutations": len(successful_mutations) - destabilizing_count
            }
        }

        logger.info("ddG calculations completed successfully!")
        logger.info(f"Successful mutations: {len(successful_mutations)}/{len(mutation_list)}")
        if successful_mutations:
            logger.info(f"Average ΔΔG: {average_ddg:.2f} kcal/mol")
            logger.info(f"Destabilizing: {destabilizing_count}/{len(successful_mutations)}")

        # Save output if requested
        output_path = None
        if output_file:
            import json
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            output_data = {
                "ddg_results": result,
                "metadata": {
                    "input_file": str(input_file),
                    "mutations": mutation_list,
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
                "mutations": mutation_list,
                "config": config,
                "success": True
            }
        }

    except Exception as e:
        logger.error(f"Error during ddG calculations: {e}")
        return {
            "result": None,
            "output_file": None,
            "metadata": {
                "input_file": str(input_file),
                "mutations": mutation_list if 'mutation_list' in locals() else mutations,
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
    parser.add_argument('--mutations', '-m', help='Mutations (e.g., "A10G,B15L")')
    parser.add_argument('--trajectories', '-n', type=int, help='Number of trajectories')
    parser.add_argument('--repack_radius', '-r', type=float, help='Repacking radius (Angstroms)')
    parser.add_argument('--no-minimize', action='store_true', help='Skip minimization')
    parser.add_argument('--cartesian', action='store_true', help='Use Cartesian minimization')
    parser.add_argument('--demo', action='store_true', help='Run with demo data')

    args = parser.parse_args()

    # Handle demo mode
    if args.demo:
        demo_pdb = Path('examples/data/test_complex.pdb')
        if not demo_pdb.exists():
            logger = setup_logging()
            logger.warning(f"Demo file {demo_pdb} not found. Creating minimal demo structure...")
            create_demo_structure(demo_pdb)
        args.input = str(demo_pdb)
        if not args.mutations:
            args.mutations = "A1G"  # Demo mutation
        if not args.trajectories:
            args.trajectories = 2  # Smaller for demo

    # Validate required parameters
    if not args.input:
        print("Error: Input file required! Use --input or --demo")
        return 1
    if not args.mutations:
        print("Error: Mutations required! Use --mutations or --demo")
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
    if args.repack_radius:
        config_overrides['repack_radius'] = args.repack_radius
    if args.no_minimize:
        config_overrides['minimize'] = False
    if args.cartesian:
        config_overrides['use_cartesian'] = True

    # Run ddG calculations
    result = run_ddg_calculations(
        input_file=args.input,
        mutations=args.mutations,
        output_file=args.output,
        config=config,
        **config_overrides
    )

    if result['metadata']['success']:
        print(f"✅ Success: ddG calculations completed")
        if result['output_file']:
            print(f"📄 Results saved to: {result['output_file']}")
        if result['result']:
            summary = result['result']['summary']
            print(f"📊 Successful: {summary['successful_mutations']}/{summary['total_mutations']}")
            if summary['successful_mutations'] > 0:
                print(f"📈 Average ΔΔG: {summary['average_ddg']:.2f} kcal/mol")
                print(f"⚠️ Destabilizing: {summary['destabilizing_mutations']}")
                print(f"✅ Stabilizing: {summary['stabilizing_mutations']}")
        return 0
    else:
        print(f"❌ Error: {result['metadata']['error']}")
        return 1

if __name__ == '__main__':
    sys.exit(main())