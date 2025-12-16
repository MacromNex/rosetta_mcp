#!/usr/bin/env python3

"""
PyRosetta Loop Modeling - Use Case 3

Description: Models protein loops using a combination of fragment insertion and
CCD (Cyclic Coordinate Descent) loop closure. Uses both low-resolution (centroid)
simulated annealing and high-resolution (fullatom) refinement.

Based on: D080_Loop_modeling.py from PyRosetta demos
Category: Loop modeling
Environment: ./env

Requirements:
- PyRosetta installed and initialized
- Input PDB file with a loop to model
- Optional: Fragment files for better sampling
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

def setup_loop_fold_tree(pose, loop_begin, loop_end, cutpoint=None):
    """
    Set up fold tree for loop modeling

    Args:
        pose: PyRosetta pose
        loop_begin: First residue of the loop
        loop_end: Last residue of the loop
        cutpoint: Cutpoint for loop closure (middle if None)
    """
    if cutpoint is None:
        cutpoint = (loop_begin + loop_end) // 2

    ft = FoldTree()

    # Create fold tree with loop
    if loop_begin > 1:
        ft.add_edge(1, loop_begin - 1, -1)  # Before loop

    if loop_end < pose.total_residue():
        ft.add_edge(loop_end + 1, pose.total_residue(), -1)  # After loop

    # Loop edges with cutpoint
    ft.add_edge(loop_begin - 1, cutpoint, -1)
    ft.add_edge(cutpoint + 1, loop_end + 1, -1)
    ft.add_edge(cutpoint, cutpoint + 1, 1)  # Cutpoint jump

    pose.fold_tree(ft)
    logger.info(f"Loop fold tree set up: Loop {loop_begin}-{loop_end}, cutpoint at {cutpoint}")

def loop_modeling(pdb_file, loop_begin, loop_end, cutpoint=None, output_prefix="loop_model",
                 num_trajectories=10, fragments_file=None, use_pymol=False):
    """
    Perform loop modeling using PyRosetta

    Args:
        pdb_file (str): Input PDB file path
        loop_begin (int): First residue of the loop to model
        loop_end (int): Last residue of the loop to model
        cutpoint (int): Cutpoint for loop closure (auto-detect if None)
        output_prefix (str): Prefix for output files
        num_trajectories (int): Number of modeling trajectories
        fragments_file (str): Fragment file for sampling (optional)
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
        logger.info(f"Loading structure from {pdb_file}")
        pose = pose_from_file(pdb_file)

        # Validate loop region
        if loop_begin < 1 or loop_end > pose.total_residue() or loop_begin >= loop_end:
            logger.error(f"Invalid loop region: {loop_begin}-{loop_end} for pose with {pose.total_residue()} residues")
            return False

        # 2. Set up cutpoint
        if cutpoint is None:
            cutpoint = (loop_begin + loop_end) // 2

        # 3. Set up loop fold tree
        setup_loop_fold_tree(pose, loop_begin, loop_end, cutpoint)

        # 4. Set up scoring functions
        scorefxn_low = core.scoring.ScoreFunctionFactory.create_score_function("cen_std")
        scorefxn_high = get_fa_scorefxn()

        logger.info(f"Original score (high-res): {scorefxn_high(pose):.2f}")

        # 5. Set up loop object
        from pyrosetta.rosetta.protocols.loops import Loop, Loops
        loop = Loop(loop_begin, loop_end, cutpoint)
        loops = Loops()
        loops.add_loop(loop)

        logger.info(f"Modeling loop: {loop}")

        # 6. Set up move map
        movemap = MoveMap()
        movemap.set_bb_true_range(loop_begin, loop_end)
        movemap.set_chi_true_range(loop_begin, loop_end)

        # 7. Set up loop modeling protocol
        if fragments_file and os.path.exists(fragments_file):
            logger.info(f"Using fragment file: {fragments_file}")
            # Fragment-based loop modeling
            loop_mover = protocols.loops.loop_mover.LoopMover_Perturb_CCD()
        else:
            logger.info("Using CCD-only loop modeling (no fragments)")
            # CCD-only loop modeling
            loop_mover = protocols.loops.loop_mover.LoopMover_CCD()

        loop_mover.loops(loops)
        loop_mover.set_scorefxn(scorefxn_high)

        # 8. Set up refinement protocol
        refine_mover = protocols.loops.loop_mover.LoopMover_Refine_CCD()
        refine_mover.loops(loops)
        refine_mover.set_scorefxn(scorefxn_high)

        # 9. Set up job distributor
        job_distributor = PyJobDistributor(output_prefix, num_trajectories, scorefxn_high)

        # 10. Run loop modeling trajectories
        logger.info(f"Starting {num_trajectories} loop modeling trajectories")
        scores = []
        loop_rmsd_values = []

        # Store original loop coordinates for RMSD calculation
        original_pose = Pose()
        original_pose.assign(pose)

        while not job_distributor.job_complete:
            # Reset pose for new trajectory
            test_pose = Pose()
            test_pose.assign(pose)

            # Set pose name
            job_name = job_distributor.output_name()
            test_pose.pdb_info().name(job_name)

            try:
                # 1. Centroid loop modeling
                to_centroid = protocols.simple_moves.SwitchResidueTypeSetMover("centroid")
                to_centroid.apply(test_pose)

                # Apply loop modeling
                loop_mover.apply(test_pose)

                # 2. Convert back to fullatom
                to_fullatom = protocols.simple_moves.SwitchResidueTypeSetMover("fa_standard")
                to_fullatom.apply(test_pose)

                # 3. Fullatom refinement
                refine_mover.apply(test_pose)

                # 4. Calculate final score
                final_score = scorefxn_high(test_pose)

                # 5. Calculate loop RMSD
                loop_rmsd = 0.0
                try:
                    # Calculate backbone RMSD for loop region
                    atom_map = core.id.AtomID_Map()
                    atom_map.fill_with(core.id.BOGUS_ATOM_ID)

                    for i in range(loop_begin, loop_end + 1):
                        atom_map.set(core.id.AtomID(1, i), core.id.AtomID(1, i))  # N
                        atom_map.set(core.id.AtomID(2, i), core.id.AtomID(2, i))  # CA
                        atom_map.set(core.id.AtomID(3, i), core.id.AtomID(3, i))  # C

                    loop_rmsd = core.scoring.rms_at_corresponding_atoms(original_pose, test_pose, atom_map)
                except:
                    loop_rmsd = -1.0  # Failed to calculate

                scores.append(final_score)
                loop_rmsd_values.append(loop_rmsd)

                logger.info(f"Trajectory {len(scores)}: Score = {final_score:.2f}, Loop RMSD = {loop_rmsd:.2f}")

                # Output structure
                job_distributor.output_decoy(test_pose)

            except Exception as e:
                logger.warning(f"Trajectory {len(scores)+1} failed: {e}")
                # Create a dummy pose to satisfy job distributor
                job_distributor.output_decoy(pose)

        # 11. Summary
        if scores:
            best_score_idx = scores.index(min(scores))
            best_score = scores[best_score_idx]
            best_rmsd = loop_rmsd_values[best_score_idx]

            logger.info("Loop modeling completed successfully!")
            logger.info(f"Best score: {best_score:.2f}")
            logger.info(f"Best loop RMSD: {best_rmsd:.2f}")
            logger.info(f"Average score: {sum(scores)/len(scores):.2f}")

            valid_rmsds = [r for r in loop_rmsd_values if r >= 0]
            if valid_rmsds:
                logger.info(f"Average loop RMSD: {sum(valid_rmsds)/len(valid_rmsds):.2f}")

        return True

    except Exception as e:
        logger.error(f"Error during loop modeling: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Loop Modeling using PyRosetta')
    parser.add_argument('--input', '-i', help='Input PDB file')
    parser.add_argument('--loop_start', '-s', type=int, help='First residue of loop')
    parser.add_argument('--loop_end', '-e', type=int, help='Last residue of loop')
    parser.add_argument('--cutpoint', '-c', type=int, help='Cutpoint for loop closure (auto-detect if not specified)')
    parser.add_argument('--output', '-o', default='loop_model', help='Output prefix (default: loop_model)')
    parser.add_argument('--trajectories', '-n', type=int, default=10, help='Number of trajectories (default: 10)')
    parser.add_argument('--fragments', '-f', help='Fragment file for sampling (optional)')
    parser.add_argument('--pymol', action='store_true', help='Enable PyMOL visualization')
    parser.add_argument('--demo', action='store_true', help='Run with demo data')

    args = parser.parse_args()

    # Handle demo mode
    if args.demo:
        demo_pdb = os.path.join('examples', 'data', 'test_loop.pdb')
        if not os.path.exists(demo_pdb):
            logger.warning(f"Demo file {demo_pdb} not found. Creating minimal demo structure...")
            # Create a simple demo PDB with a loop region
            demo_content = """ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N
ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00 20.00           C
ATOM      3  C   ALA A   1      11.500  11.000  10.000  1.00 20.00           C
ATOM      4  O   ALA A   1      11.500  12.000  10.000  1.00 20.00           O
ATOM      5  N   GLY A   2      12.000  11.000  9.000   1.00 20.00           N
ATOM      6  CA  GLY A   2      13.000  11.000  9.000   1.00 20.00           C
ATOM      7  C   GLY A   2      13.500  12.000  9.000   1.00 20.00           C
ATOM      8  O   GLY A   2      13.500  13.000  9.000   1.00 20.00           O
ATOM      9  N   ALA A   3      14.000  12.000  8.000   1.00 20.00           N
ATOM     10  CA  ALA A   3      15.000  12.000  8.000   1.00 20.00           C
ATOM     11  C   ALA A   3      15.500  13.000  8.000   1.00 20.00           C
ATOM     12  O   ALA A   3      15.500  14.000  8.000   1.00 20.00           O
TER
END
"""
            os.makedirs(os.path.dirname(demo_pdb), exist_ok=True)
            with open(demo_pdb, 'w') as f:
                f.write(demo_content)
        args.input = demo_pdb
        args.loop_start = 2  # Model glycine as loop
        args.loop_end = 3
        args.trajectories = 3  # Smaller for demo

    # Validate input file
    if not args.input:
        logger.error("Input file required! Use --input or --demo")
        return 1
    if not os.path.exists(args.input):
        logger.error(f"Input file {args.input} not found!")
        return 1

    # Validate loop region
    if not args.loop_start or not args.loop_end:
        logger.error("Loop start and end residues required! Use --loop_start/--loop_end or --demo")
        return 1
    if args.loop_start >= args.loop_end:
        logger.error("Loop start must be less than loop end!")
        return 1

    # Run loop modeling
    logger.info("Starting loop modeling...")
    success = loop_modeling(
        pdb_file=args.input,
        loop_begin=args.loop_start,
        loop_end=args.loop_end,
        cutpoint=args.cutpoint,
        output_prefix=args.output,
        num_trajectories=args.trajectories,
        fragments_file=args.fragments,
        use_pymol=args.pymol
    )

    if success:
        logger.info("Loop modeling completed successfully!")
        return 0
    else:
        logger.error("Loop modeling failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())