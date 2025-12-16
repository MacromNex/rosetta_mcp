#!/usr/bin/env python3

"""
PyRosetta Protein Refinement (Relax) - Use Case 1

Description: Performs high-resolution (fullatom) refinement on an input protein structure.
The protocol uses small backbone torsion angle perturbations followed by backbone minimization
and sidechain packing to improve the structure quality.

Based on: D070_Refinement.py from PyRosetta demos
Category: Relax, Structure quality analysis
Environment: ./env

Requirements:
- PyRosetta installed and initialized
- Input PDB file
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
            init, Pose, pose_from_file, get_fa_scorefxn, standard_packer_task,
            SequenceMover, MonteCarlo, TrialMover, RepeatMover, MoveMap,
            PyMOLMover, PyJobDistributor
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

def protein_refinement(pdb_file, output_prefix="refined", num_trajectories=5, num_cycles=100, use_pymol=False):
    """
    Perform protein structure refinement using PyRosetta

    Args:
        pdb_file (str): Input PDB file path
        output_prefix (str): Prefix for output files
        num_trajectories (int): Number of refinement trajectories
        num_cycles (int): Number of Monte Carlo cycles per trajectory
        use_pymol (bool): Enable PyMOL visualization

    Returns:
        bool: Success status
    """
    success, pyrosetta, rosetta = setup_pyrosetta()
    if not success:
        return False

    try:
        from pyrosetta import (
            Pose, pose_from_file, get_fa_scorefxn, standard_packer_task,
            SequenceMover, MonteCarlo, TrialMover, RepeatMover, MoveMap,
            PyMOLMover, PyJobDistributor
        )
        from pyrosetta.rosetta import core, protocols

        # 1. Load the input structure
        logger.info(f"Loading structure from {pdb_file}")
        pose = pose_from_file(pdb_file)
        original_pose = Pose()
        original_pose.assign(pose)

        # 2. Set up scoring function
        score_function = get_fa_scorefxn()
        original_score = score_function(pose)
        logger.info(f"Original score: {original_score:.2f}")

        # 3. Set up move map (allow all backbone torsions to move)
        movemap = MoveMap()
        movemap.set_bb(True)
        movemap.set_chi(True)

        # 4. Set up movers
        # Small backbone moves
        small_mover = protocols.simple_moves.SmallMover(movemap, 0.8, 1)
        shear_mover = protocols.simple_moves.ShearMover(movemap, 0.8, 1)

        # Minimization
        min_mover = protocols.simple_moves.MinMover()
        min_mover.movemap(movemap)
        min_mover.score_function(score_function)

        # Sidechain packing
        task_factory = core.pack.task.TaskFactory()
        packer = protocols.simple_moves.PackRotamersMover(score_function)
        packer.task_factory(task_factory)

        # 5. Set up PyMOL mover if requested
        pymol_mover = None
        if use_pymol:
            try:
                pymol_mover = PyMOLMover()
                pymol_mover.keep_history(True)
                logger.info("PyMOL mover initialized")
            except:
                logger.warning("PyMOL mover failed to initialize, continuing without visualization")

        # 6. Create sequence of moves
        sequence_mover = SequenceMover()
        sequence_mover.add_mover(small_mover)
        sequence_mover.add_mover(shear_mover)
        sequence_mover.add_mover(min_mover)
        sequence_mover.add_mover(packer)
        if pymol_mover:
            sequence_mover.add_mover(pymol_mover)

        # 7. Set up Monte Carlo
        mc = MonteCarlo(pose, score_function, 2.0)  # kT = 2.0

        # 8. Create trial mover
        trial_mover = TrialMover(sequence_mover, mc)

        # 9. Repeat mover
        repeat_mover = RepeatMover(trial_mover, num_cycles)

        # 10. Set up job distributor
        job_distributor = PyJobDistributor(output_prefix, num_trajectories, score_function)

        # 11. Run refinement trajectories
        logger.info(f"Starting {num_trajectories} refinement trajectories with {num_cycles} cycles each")
        scores = []

        while not job_distributor.job_complete:
            # Reset pose for new trajectory
            pose.assign(original_pose)

            # Reset Monte Carlo
            mc.reset(pose)

            # Set pose name for PyMOL
            job_name = job_distributor.output_name()
            pose.pdb_info().name(job_name)

            if pymol_mover:
                pymol_mover.apply(pose)

            # Run refinement protocol
            repeat_mover.apply(pose)

            # Get final score
            final_score = score_function(pose)
            scores.append(final_score)

            logger.info(f"Trajectory {len(scores)}: Score = {final_score:.2f}")

            # Output structure
            job_distributor.output_decoy(pose)

        # 12. Summary
        best_score = min(scores)
        improvement = original_score - best_score

        logger.info("Refinement completed successfully!")
        logger.info(f"Original score: {original_score:.2f}")
        logger.info(f"Best refined score: {best_score:.2f}")
        logger.info(f"Score improvement: {improvement:.2f}")
        logger.info(f"Average score: {sum(scores)/len(scores):.2f}")

        return True

    except Exception as e:
        logger.error(f"Error during protein refinement: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description='Protein Structure Refinement using PyRosetta')
    parser.add_argument('--input', '-i', help='Input PDB file')
    parser.add_argument('--output', '-o', default='refined', help='Output prefix (default: refined)')
    parser.add_argument('--trajectories', '-n', type=int, default=5, help='Number of trajectories (default: 5)')
    parser.add_argument('--cycles', '-c', type=int, default=100, help='Monte Carlo cycles per trajectory (default: 100)')
    parser.add_argument('--pymol', action='store_true', help='Enable PyMOL visualization')
    parser.add_argument('--demo', action='store_true', help='Run with demo data')

    args = parser.parse_args()

    # Handle demo mode
    if args.demo:
        demo_pdb = os.path.join('examples', 'data', 'test_input.pdb')
        if not os.path.exists(demo_pdb):
            logger.warning(f"Demo file {demo_pdb} not found. Creating minimal demo structure...")
            # Create a simple demo PDB (just for testing structure)
            demo_content = """ATOM      1  N   ALA A   1      20.154  11.200  19.962  1.00 20.00           N
ATOM      2  CA  ALA A   1      19.030  11.200  20.875  1.00 20.00           C
ATOM      3  C   ALA A   1      17.662  11.200  20.146  1.00 20.00           C
ATOM      4  O   ALA A   1      17.662  11.200  18.921  1.00 20.00           O
ATOM      5  CB  ALA A   1      19.030  12.424  21.788  1.00 20.00           C
TER
END
"""
            os.makedirs(os.path.dirname(demo_pdb), exist_ok=True)
            with open(demo_pdb, 'w') as f:
                f.write(demo_content)
        args.input = demo_pdb
        args.trajectories = 2  # Smaller for demo
        args.cycles = 10

    # Validate input file
    if not args.input:
        logger.error("Input file required! Use --input or --demo")
        return 1
    if not os.path.exists(args.input):
        logger.error(f"Input file {args.input} not found!")
        return 1

    # Run refinement
    logger.info("Starting protein structure refinement...")
    success = protein_refinement(
        pdb_file=args.input,
        output_prefix=args.output,
        num_trajectories=args.trajectories,
        num_cycles=args.cycles,
        use_pymol=args.pymol
    )

    if success:
        logger.info("Protein refinement completed successfully!")
        return 0
    else:
        logger.error("Protein refinement failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())