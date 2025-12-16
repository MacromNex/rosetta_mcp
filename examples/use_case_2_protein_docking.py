#!/usr/bin/env python3

"""
PyRosetta Protein-Protein Docking - Use Case 2

Description: Performs protein-protein docking to predict the conformations of two binding proteins.
Uses Rosetta's DockingProtocol with both low-resolution (centroid) and high-resolution (fullatom) stages.

Based on: D100_Docking.py from PyRosetta demos
Category: Protein-protein docking
Environment: ./env

Requirements:
- PyRosetta installed and initialized
- Input PDB file containing both protein chains
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
            FoldTree, MoveMap, PyJobDistributor
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

def setup_docking_foldtree(pose, chain_break):
    """
    Set up fold tree for docking

    Args:
        pose: PyRosetta pose
        chain_break: Residue number where chain A ends
    """
    # Create new fold tree for docking
    ft = FoldTree()

    # Chain A: residues 1 to chain_break
    # Chain B: residues chain_break+1 to total_residue
    # Jump between chains at the interface

    total_res = pose.total_residue()

    ft.add_edge(1, chain_break, -1)  # Chain A backbone
    ft.add_edge(chain_break + 1, total_res, -1)  # Chain B backbone
    ft.add_edge(chain_break, chain_break + 1, 1)  # Jump between chains

    pose.fold_tree(ft)
    logger.info(f"Docking fold tree set up: Chain A (1-{chain_break}), Chain B ({chain_break + 1}-{total_res})")

def protein_docking(pdb_file, output_prefix="docked", num_trajectories=10, chain_break=None, use_pymol=False):
    """
    Perform protein-protein docking using PyRosetta

    Args:
        pdb_file (str): Input PDB file path with both protein chains
        output_prefix (str): Prefix for output files
        num_trajectories (int): Number of docking trajectories
        chain_break (int): Residue number where chain A ends (auto-detect if None)
        use_pymol (bool): Enable PyMOL visualization

    Returns:
        bool: Success status
    """
    success, pyrosetta, rosetta = setup_pyrosetta()
    if not success:
        return False

    try:
        from pyrosetta import (
            Pose, pose_from_file, get_fa_scorefxn,
            FoldTree, MoveMap, PyJobDistributor
        )
        from pyrosetta.rosetta import core, protocols

        # 1. Load the input structure
        logger.info(f"Loading protein complex from {pdb_file}")
        pose = pose_from_file(pdb_file)

        # 2. Determine chain break if not provided
        if chain_break is None:
            # Simple heuristic: find the largest gap in residue numbering
            chain_break = pose.total_residue() // 2  # Default to middle
            logger.info(f"Auto-detected chain break at residue {chain_break}")

        # 3. Set up docking fold tree
        setup_docking_foldtree(pose, chain_break)

        # 4. Set up scoring functions
        scorefxn_low = core.scoring.ScoreFunctionFactory.create_score_function("interchain_cen")
        scorefxn_high = get_fa_scorefxn()

        # Add interface weight
        scorefxn_high.set_weight(core.scoring.fa_rep, 0.1)
        scorefxn_high.set_weight(core.scoring.fa_atr, 0.8)

        logger.info(f"Original score: {scorefxn_high(pose):.2f}")

        # 5. Set up movers for centroid/fullatom conversion
        to_centroid = protocols.simple_moves.SwitchResidueTypeSetMover("centroid")
        to_fullatom = protocols.simple_moves.SwitchResidueTypeSetMover("fa_standard")
        recover_sidechains = protocols.simple_moves.ReturnSidechainMover(pose)

        # 6. Set up rigid body perturbation movers
        randomize1 = protocols.rigid.RigidBodyRandomizeMover(pose, 1, protocols.rigid.partner_upstream)
        randomize2 = protocols.rigid.RigidBodyRandomizeMover(pose, 1, protocols.rigid.partner_downstream)

        # Perturbation parameters
        perturb = protocols.rigid.RigidBodyPerturbMover(1, 3.0, 8.0)  # 3Å translation, 8° rotation
        spin = protocols.rigid.RigidBodySpinMover(1)
        slide_into_contact = protocols.docking.DockingSlideIntoContactMover(1)

        # 7. Set up minimization
        movemap = MoveMap()
        movemap.set_jump(1, True)  # Allow the jump to move
        min_mover = protocols.simple_moves.MinMover()
        min_mover.movemap(movemap)
        min_mover.score_function(scorefxn_high)

        # 8. Set up docking protocol
        dock_protocol = protocols.docking.DockingProtocol()
        dock_protocol.set_scorefxn(scorefxn_high)
        dock_protocol.set_scorefxn_low(scorefxn_low)

        # 9. Set up job distributor
        job_distributor = PyJobDistributor(output_prefix, num_trajectories, scorefxn_high)

        # 10. Run docking trajectories
        logger.info(f"Starting {num_trajectories} docking trajectories")
        scores = []
        interface_scores = []

        while not job_distributor.job_complete:
            # Reset pose for new trajectory
            test_pose = Pose()
            test_pose.assign(pose)

            # Convert to centroid
            to_centroid.apply(test_pose)

            # Set pose name
            job_name = job_distributor.output_name()
            test_pose.pdb_info().name(job_name)

            # 1. Randomize starting position
            randomize1.apply(test_pose)
            randomize2.apply(test_pose)
            spin.apply(test_pose)

            # 2. Slide partners into contact
            slide_into_contact.apply(test_pose)

            # 3. Convert back to fullatom
            to_fullatom.apply(test_pose)
            recover_sidechains.apply(test_pose)

            # 4. Run docking protocol
            dock_protocol.apply(test_pose)

            # 5. Calculate scores
            final_score = scorefxn_high(test_pose)

            # Calculate interface score (separation energy)
            separated_pose = Pose()
            separated_pose.assign(test_pose)

            # Move chains far apart to calculate separation energy
            trans_mover = protocols.rigid.RigidBodyTransMover(separated_pose, 1)
            trans_mover.step_size(1000.0)  # Move 1000 Å apart
            trans_mover.apply(separated_pose)

            separated_score = scorefxn_high(separated_pose)
            interface_score = final_score - separated_score

            scores.append(final_score)
            interface_scores.append(interface_score)

            logger.info(f"Trajectory {len(scores)}: Score = {final_score:.2f}, Interface = {interface_score:.2f}")

            # Output structure
            job_distributor.output_decoy(test_pose)

        # 11. Summary
        best_score_idx = interface_scores.index(min(interface_scores))
        best_interface_score = interface_scores[best_score_idx]
        best_total_score = scores[best_score_idx]

        logger.info("Protein docking completed successfully!")
        logger.info(f"Best interface score: {best_interface_score:.2f}")
        logger.info(f"Best total score: {best_total_score:.2f}")
        logger.info(f"Average interface score: {sum(interface_scores)/len(interface_scores):.2f}")

        return True

    except Exception as e:
        logger.error(f"Error during protein docking: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Protein-Protein Docking using PyRosetta')
    parser.add_argument('--input', '-i', help='Input PDB file with both protein chains')
    parser.add_argument('--output', '-o', default='docked', help='Output prefix (default: docked)')
    parser.add_argument('--trajectories', '-n', type=int, default=10, help='Number of trajectories (default: 10)')
    parser.add_argument('--chain_break', '-c', type=int, help='Residue number where chain A ends (auto-detect if not specified)')
    parser.add_argument('--pymol', action='store_true', help='Enable PyMOL visualization')
    parser.add_argument('--demo', action='store_true', help='Run with demo data')

    args = parser.parse_args()

    # Handle demo mode
    if args.demo:
        demo_pdb = os.path.join('examples', 'data', 'test_docking.pdb')
        if not os.path.exists(demo_pdb):
            logger.warning(f"Demo file {demo_pdb} not found. Creating minimal demo structure...")
            # Create a simple two-chain demo PDB
            demo_content = """ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      11.500  12.000  10.000  1.00 20.00           O
TER
ATOM      5  N   ALA B   2      20.000  10.000  10.000  1.00 20.00           N
ATOM      6  CA  ALA B   2      21.000  10.000  10.000  1.00 20.00           C
ATOM      7  C   ALA B   2      21.500  11.000  10.000  1.00 20.00           C
ATOM      8  O   ALA B   2      21.500  12.000  10.000  1.00 20.00           O
TER
END
"""
            os.makedirs(os.path.dirname(demo_pdb), exist_ok=True)
            with open(demo_pdb, 'w') as f:
                f.write(demo_content)
            args.chain_break = 1  # Chain A has 1 residue
        args.input = demo_pdb
        args.trajectories = 3  # Smaller for demo

    # Validate input file
    if not args.input:
        logger.error("Input file required! Use --input or --demo")
        return 1
    if not os.path.exists(args.input):
        logger.error(f"Input file {args.input} not found!")
        return 1

    # Run docking
    logger.info("Starting protein-protein docking...")
    success = protein_docking(
        pdb_file=args.input,
        output_prefix=args.output,
        num_trajectories=args.trajectories,
        chain_break=args.chain_break,
        use_pymol=args.pymol
    )

    if success:
        logger.info("Protein docking completed successfully!")
        return 0
    else:
        logger.error("Protein docking failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())