#!/usr/bin/env python3

"""
PyRosetta ddG Calculations (Alanine Scanning) - Use Case 5

Description: Performs alanine scanning mutagenesis to calculate binding energy changes (ddG)
by iteratively mutating each interface residue to alanine and determining the change in
interaction energy. Can be used for protein-protein or protein-ligand interfaces.

Based on: D090_Ala_scan.py from PyRosetta demos
Category: ddG calculations, Structure quality analysis
Environment: ./env

Requirements:
- PyRosetta installed and initialized
- Input PDB file with protein complex or protein-ligand complex
"""

import os
import sys
import argparse
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def setup_pyrosetta():
    """Initialize PyRosetta with standard settings"""
    try:
        import pyrosetta
        import pyrosetta.rosetta as rosetta
        from pyrosetta import (
            init, Pose, pose_from_file, get_fa_scorefxn,
            MoveMap
        )
        from pyrosetta.rosetta import core, protocols

        # Initialize PyRosetta with constant seed for reproducibility
        init(extra_options="-constant_seed")
        logger.info("PyRosetta initialized successfully")
        return True, pyrosetta, rosetta
    except ImportError as e:
        logger.error(f"PyRosetta import failed: {e}")
        logger.error("Please install PyRosetta: conda install -c rosettacommons pyrosetta")
        return False, None, None

def find_interface_residues(pose, chain1, chain2, distance_cutoff=10.0):
    """
    Find residues at the interface between two chains

    Args:
        pose: PyRosetta pose
        chain1: First chain ID
        chain2: Second chain ID
        distance_cutoff: Maximum distance for interface definition

    Returns:
        list: Interface residues from chain1
    """
    interface_residues = []

    # Get residues from each chain
    chain1_residues = []
    chain2_residues = []

    for i in range(1, pose.total_residue() + 1):
        chain = pose.pdb_info().chain(i)
        if chain == chain1:
            chain1_residues.append(i)
        elif chain == chain2:
            chain2_residues.append(i)

    # Find interface residues from chain1
    for res1 in chain1_residues:
        for res2 in chain2_residues:
            # Calculate distance between CA atoms
            ca1 = pose.residue(res1).atom("CA")
            ca2 = pose.residue(res2).atom("CA")
            distance = ca1.distance(ca2)

            if distance <= distance_cutoff:
                if res1 not in interface_residues:
                    interface_residues.append(res1)
                break

    return sorted(interface_residues)

def mutate_residue(pose, position, target_aa="ALA"):
    """
    Mutate a residue to target amino acid

    Args:
        pose: PyRosetta pose
        position: Residue position to mutate
        target_aa: Target amino acid (3-letter code)

    Returns:
        bool: Success status
    """
    try:
        from pyrosetta.rosetta import core

        # Get the target amino acid type
        aa_type = core.chemical.aa_from_name(target_aa)

        # Create mutator
        mutator = protocols.simple_moves.MutateResidue(position, aa_type)
        mutator.apply(pose)

        return True
    except Exception as e:
        logger.error(f"Failed to mutate residue {position} to {target_aa}: {e}")
        return False

def pack_around_residue(pose, position, scorefxn, radius=8.0):
    """
    Pack rotamers around a specific residue

    Args:
        pose: PyRosetta pose
        position: Central residue position
        scorefxn: Scoring function
        radius: Packing radius in Angstroms
    """
    try:
        from pyrosetta.rosetta import core

        # Create task factory
        task_factory = core.pack.task.TaskFactory()

        # Create operation to restrict packing to residues near the mutated position
        center = pose.residue(position).atom("CA").xyz()

        for i in range(1, pose.total_residue() + 1):
            if i == position:
                continue

            ca_pos = pose.residue(i).atom("CA").xyz()
            distance = center.distance(ca_pos)

            if distance > radius:
                task_factory.push_back(core.pack.task.operation.PreventRepackingRLT([i]))

        # Create packer
        packer = protocols.simple_moves.PackRotamersMover(scorefxn)
        packer.task_factory(task_factory)

        # Apply packing
        packer.apply(pose)

    except Exception as e:
        logger.warning(f"Packing around residue {position} failed: {e}")

def calculate_binding_energy(pose, chain1, chain2, scorefxn):
    """
    Calculate binding energy by separation method

    Args:
        pose: PyRosetta pose
        chain1: First chain ID
        chain2: Second chain ID
        scorefxn: Scoring function

    Returns:
        float: Binding energy (complex_score - separated_scores)
    """
    # Score the complex
    complex_score = scorefxn(pose)

    # Create separated pose
    separated_pose = pose.clone()

    # Find jump between chains
    jump_num = 1  # Assume first jump separates the chains

    # Move chains apart by 1000 Angstroms
    trans_mover = protocols.rigid.RigidBodyTransMover(separated_pose, jump_num)
    trans_mover.step_size(1000.0)
    trans_mover.apply(separated_pose)

    # Score the separated complex
    separated_score = scorefxn(separated_pose)

    return complex_score - separated_score

def ddg_calculation(pdb_file, chain1="A", chain2="B", output_file="ddg_results.txt",
                   target_residue=None, distance_cutoff=10.0, mutation_aa="ALA"):
    """
    Perform ddG calculations using alanine scanning

    Args:
        pdb_file (str): Input PDB file path
        chain1 (str): First chain ID (will be scanned for mutations)
        chain2 (str): Second chain ID (partner)
        output_file (str): Output file for results
        target_residue (int): Specific residue to mutate (scan all interface if None)
        distance_cutoff (float): Interface definition cutoff
        mutation_aa (str): Target amino acid for mutation

    Returns:
        bool: Success status
    """
    success, pyrosetta, rosetta = setup_pyrosetta()
    if not success:
        return False

    try:
        from pyrosetta import Pose, pose_from_file, get_fa_scorefxn
        from pyrosetta.rosetta import core, protocols

        # 1. Load the input structure
        logger.info(f"Loading complex from {pdb_file}")
        original_pose = pose_from_file(pdb_file)

        # 2. Set up scoring function
        scorefxn = get_fa_scorefxn()
        original_score = scorefxn(original_pose)
        logger.info(f"Original complex score: {original_score:.2f}")

        # 3. Calculate original binding energy
        original_binding = calculate_binding_energy(original_pose, chain1, chain2, scorefxn)
        logger.info(f"Original binding energy: {original_binding:.2f}")

        # 4. Determine residues to mutate
        if target_residue:
            residues_to_scan = [target_residue]
            logger.info(f"Scanning single residue: {target_residue}")
        else:
            residues_to_scan = find_interface_residues(original_pose, chain1, chain2, distance_cutoff)
            logger.info(f"Found {len(residues_to_scan)} interface residues in chain {chain1}")
            logger.info(f"Interface residues: {residues_to_scan}")

        if not residues_to_scan:
            logger.error("No interface residues found!")
            return False

        # 5. Perform mutations and calculate ddG
        results = []

        with open(output_file, 'w') as f:
            f.write("# ddG Calculation Results\\n")
            f.write("# Position\\tOriginal_AA\\tMutant_AA\\tWild_Binding\\tMutant_Binding\\tddG\\n")

            for position in residues_to_scan:
                logger.info(f"\\nMutating position {position}...")

                try:
                    # Get original amino acid
                    original_aa = original_pose.residue(position).name3()

                    # Skip if already the target amino acid
                    if original_aa == mutation_aa:
                        logger.info(f"Position {position} is already {mutation_aa}, skipping")
                        continue

                    # Create mutant pose
                    mutant_pose = original_pose.clone()

                    # Mutate residue
                    if not mutate_residue(mutant_pose, position, mutation_aa):
                        logger.error(f"Failed to mutate position {position}")
                        continue

                    # Pack around mutation site
                    pack_around_residue(mutant_pose, position, scorefxn)

                    # Calculate mutant binding energy
                    mutant_binding = calculate_binding_energy(mutant_pose, chain1, chain2, scorefxn)

                    # Calculate ddG
                    ddg = mutant_binding - original_binding

                    # Store results
                    result = {
                        'position': position,
                        'original_aa': original_aa,
                        'mutant_aa': mutation_aa,
                        'wild_binding': original_binding,
                        'mutant_binding': mutant_binding,
                        'ddg': ddg
                    }
                    results.append(result)

                    logger.info(f"Position {position}: {original_aa} -> {mutation_aa}")
                    logger.info(f"Wild binding: {original_binding:.2f}")
                    logger.info(f"Mutant binding: {mutant_binding:.2f}")
                    logger.info(f"ddG: {ddg:.2f}")

                    # Write to output file
                    f.write(f"{position}\\t{original_aa}\\t{mutation_aa}\\t{original_binding:.2f}\\t{mutant_binding:.2f}\\t{ddg:.2f}\\n")
                    f.flush()

                except Exception as e:
                    logger.error(f"Error processing position {position}: {e}")
                    continue

        # 6. Summary
        if results:
            logger.info(f"\\nddG calculation completed for {len(results)} positions")
            logger.info(f"Results saved to {output_file}")

            # Find hotspots (ddG > 1.0 kcal/mol)
            hotspots = [r for r in results if r['ddg'] > 1.0]
            if hotspots:
                logger.info(f"\\nHotspots found (ddG > 1.0): {len(hotspots)} residues")
                for hs in hotspots:
                    logger.info(f"  Position {hs['position']} ({hs['original_aa']}): ddG = {hs['ddg']:.2f}")

            # Statistics
            ddg_values = [r['ddg'] for r in results]
            avg_ddg = sum(ddg_values) / len(ddg_values)
            max_ddg = max(ddg_values)
            min_ddg = min(ddg_values)

            logger.info(f"\\nStatistics:")
            logger.info(f"Average ddG: {avg_ddg:.2f}")
            logger.info(f"Maximum ddG: {max_ddg:.2f}")
            logger.info(f"Minimum ddG: {min_ddg:.2f}")

        return True

    except Exception as e:
        logger.error(f"Error during ddG calculation: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='ddG Calculations (Alanine Scanning) using PyRosetta')
    parser.add_argument('--input', '-i', help='Input PDB file with protein complex')
    parser.add_argument('--chain1', '-c1', default='A', help='First chain ID (will be mutated) (default: A)')
    parser.add_argument('--chain2', '-c2', default='B', help='Second chain ID (partner) (default: B)')
    parser.add_argument('--output', '-o', default='ddg_results.txt', help='Output file (default: ddg_results.txt)')
    parser.add_argument('--residue', '-r', type=int, help='Specific residue to mutate (scan all interface if not specified)')
    parser.add_argument('--distance', '-d', type=float, default=10.0, help='Interface distance cutoff (default: 10.0)')
    parser.add_argument('--mutation', '-m', default='ALA', help='Target amino acid for mutation (default: ALA)')
    parser.add_argument('--demo', action='store_true', help='Run with demo data')

    args = parser.parse_args()

    # Handle demo mode
    if args.demo:
        demo_pdb = os.path.join('examples', 'data', 'test_complex.pdb')
        if not os.path.exists(demo_pdb):
            logger.warning(f"Demo file {demo_pdb} not found. Creating minimal demo structure...")
            # Create a simple two-chain demo PDB
            demo_content = """ATOM      1  N   PHE A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  PHE A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   PHE A   1      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  O   PHE A   1      11.500  12.000  10.000  1.00 20.00           O
ATOM      5  CB  PHE A   1      11.500   9.000  10.000  1.00 20.00           C
TER
ATOM      6  N   TYR B   2      15.000  10.000  10.000  1.00 20.00           N
ATOM      7  CA  TYR B   2      16.000  10.000  10.000  1.00 20.00           C
ATOM      8  C   TYR B   2      16.500  11.000  10.000  1.00 20.00           C
ATOM      9  O   TYR B   2      16.500  12.000  10.000  1.00 20.00           O
ATOM     10  CB  TYR B   2      16.500   9.000  10.000  1.00 20.00           C
TER
END
"""
            os.makedirs(os.path.dirname(demo_pdb), exist_ok=True)
            with open(demo_pdb, 'w') as f:
                f.write(demo_content)
        args.input = demo_pdb

    # Validate input file
    if not args.input:
        logger.error("Input file required! Use --input or --demo")
        return 1
    if not os.path.exists(args.input):
        logger.error(f"Input file {args.input} not found!")
        return 1

    # Run ddG calculation
    logger.info("Starting ddG calculations...")
    success = ddg_calculation(
        pdb_file=args.input,
        chain1=args.chain1,
        chain2=args.chain2,
        output_file=args.output,
        target_residue=args.residue,
        distance_cutoff=args.distance,
        mutation_aa=args.mutation
    )

    if success:
        logger.info("ddG calculations completed successfully!")
        return 0
    else:
        logger.error("ddG calculations failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())