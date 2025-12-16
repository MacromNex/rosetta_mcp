"""MCP Server for Rosetta

Provides both synchronous and asynchronous (submit) APIs for Rosetta molecular modeling tools.
"""

from fastmcp import FastMCP
from pathlib import Path
from typing import Optional, List
import sys

# Setup paths
SCRIPT_DIR = Path(__file__).parent.resolve()
MCP_ROOT = SCRIPT_DIR.parent
SCRIPTS_DIR = MCP_ROOT / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))
sys.path.insert(0, str(SCRIPTS_DIR))

from jobs.manager import job_manager
from loguru import logger

# Create MCP server
mcp = FastMCP("rosetta")

# ==============================================================================
# Job Management Tools (for async operations)
# ==============================================================================

@mcp.tool()
def get_job_status(job_id: str) -> dict:
    """
    Get the status of a submitted job.

    Args:
        job_id: The job ID returned from a submit_* function

    Returns:
        Dictionary with job status, timestamps, and any errors
    """
    return job_manager.get_job_status(job_id)

@mcp.tool()
def get_job_result(job_id: str) -> dict:
    """
    Get the results of a completed job.

    Args:
        job_id: The job ID of a completed job

    Returns:
        Dictionary with the job results or error if not completed
    """
    return job_manager.get_job_result(job_id)

@mcp.tool()
def get_job_log(job_id: str, tail: int = 50) -> dict:
    """
    Get log output from a running or completed job.

    Args:
        job_id: The job ID to get logs for
        tail: Number of lines from end (default: 50, use 0 for all)

    Returns:
        Dictionary with log lines and total line count
    """
    return job_manager.get_job_log(job_id, tail)

@mcp.tool()
def cancel_job(job_id: str) -> dict:
    """
    Cancel a running job.

    Args:
        job_id: The job ID to cancel

    Returns:
        Success or error message
    """
    return job_manager.cancel_job(job_id)

@mcp.tool()
def list_jobs(status: Optional[str] = None) -> dict:
    """
    List all submitted jobs.

    Args:
        status: Filter by status (pending, running, completed, failed, cancelled)

    Returns:
        List of jobs with their status
    """
    return job_manager.list_jobs(status)

# ==============================================================================
# Synchronous Tools (for fast operations < 10 min)
# ==============================================================================

@mcp.tool()
def refine_protein_structure(
    input_file: str,
    trajectories: int = 5,
    cycles: int = 100,
    temperature: float = 2.0,
    output_file: Optional[str] = None
) -> dict:
    """
    Refine protein structure using Rosetta Relax protocol (fast mode).

    Fast operation that completes in ~2-5 minutes. Use this for quick refinements
    with small trajectory counts.

    Args:
        input_file: Path to input PDB file
        trajectories: Number of refinement trajectories (default: 5)
        cycles: Number of Monte Carlo cycles per trajectory (default: 100)
        temperature: Monte Carlo temperature for acceptance (default: 2.0)
        output_file: Optional path to save results JSON

    Returns:
        Dictionary with refinement results and scores
    """
    try:
        from protein_refinement import run_protein_refinement

        result = run_protein_refinement(
            input_file=input_file,
            output_file=output_file,
            trajectories=trajectories,
            cycles=cycles,
            temperature=temperature
        )
        return {"status": "success", **result}
    except FileNotFoundError as e:
        return {"status": "error", "error": f"File not found: {e}"}
    except ImportError as e:
        return {"status": "error", "error": f"PyRosetta not available: {e}"}
    except Exception as e:
        logger.error(f"Protein refinement failed: {e}")
        return {"status": "error", "error": str(e)}

@mcp.tool()
def calculate_ddg(
    input_file: str,
    mutations: str,
    trajectories: int = 5,
    repack_radius: float = 8.0,
    minimize: bool = True,
    output_file: Optional[str] = None
) -> dict:
    """
    Calculate ΔΔG (stability change) for protein mutations (fast mode).

    Fast operation that completes in ~1-3 minutes for single mutations.
    Use this for quick stability analysis.

    Args:
        input_file: Path to input PDB file
        mutations: Mutations in format "A10G,B15L" or "10G,15L"
        trajectories: Number of trajectories (default: 5)
        repack_radius: Radius for sidechain repacking in Angstroms (default: 8.0)
        minimize: Whether to minimize after mutation (default: True)
        output_file: Optional path to save results JSON

    Returns:
        Dictionary with ΔΔG values and mutation analysis
    """
    try:
        from ddg_calculations import run_ddg_calculations

        result = run_ddg_calculations(
            input_file=input_file,
            mutations=mutations,
            output_file=output_file,
            trajectories=trajectories,
            repack_radius=repack_radius,
            minimize=minimize
        )
        return {"status": "success", **result}
    except FileNotFoundError as e:
        return {"status": "error", "error": f"File not found: {e}"}
    except ImportError as e:
        return {"status": "error", "error": f"PyRosetta not available: {e}"}
    except Exception as e:
        logger.error(f"ΔΔG calculation failed: {e}")
        return {"status": "error", "error": str(e)}

# ==============================================================================
# Submit Tools (for long-running operations > 10 min)
# ==============================================================================

@mcp.tool()
def submit_protein_docking(
    input_file: str,
    trajectories: int = 10,
    chain_break: Optional[int] = None,
    use_centroid_stage: bool = True,
    use_fullatom_stage: bool = True,
    output_dir: Optional[str] = None,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit protein-protein docking for background processing.

    This operation may take >10 minutes for thorough docking searches.
    Returns a job_id for tracking.

    Args:
        input_file: Path to input PDB file with both proteins
        trajectories: Number of docking trajectories (default: 10)
        chain_break: Residue where chain A ends (auto-detect if None)
        use_centroid_stage: Use low-resolution stage (default: True)
        use_fullatom_stage: Use high-resolution stage (default: True)
        output_dir: Directory for outputs
        job_name: Optional name for tracking

    Returns:
        Dictionary with job_id. Use:
        - get_job_status(job_id) to check progress
        - get_job_result(job_id) to get results
        - get_job_log(job_id) to see logs
    """
    script_path = str(SCRIPTS_DIR / "protein_docking.py")

    args = {
        "input": input_file,
        "trajectories": trajectories
    }

    if chain_break is not None:
        args["chain_break"] = chain_break

    # Note: use_centroid_stage, use_fullatom_stage, and output_dir are controlled by script's default config
    # These parameters are noted in the function documentation but handled internally by the script

    return job_manager.submit_job(
        script_path=script_path,
        args=args,
        job_name=job_name or f"docking_{Path(input_file).stem}"
    )

@mcp.tool()
def submit_loop_modeling(
    input_file: str,
    loop_start: int,
    loop_end: int,
    trajectories: int = 10,
    use_centroid_stage: bool = True,
    use_fullatom_stage: bool = True,
    output_dir: Optional[str] = None,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit loop modeling for background processing.

    This operation may take >10 minutes for complex loop regions.
    Returns a job_id for tracking.

    Args:
        input_file: Path to input PDB file
        loop_start: Starting residue number for loop region
        loop_end: Ending residue number for loop region
        trajectories: Number of modeling trajectories (default: 10)
        use_centroid_stage: Use low-resolution stage (default: True)
        use_fullatom_stage: Use high-resolution stage (default: True)
        output_dir: Directory for outputs
        job_name: Optional name for tracking

    Returns:
        Dictionary with job_id for tracking the loop modeling job
    """
    script_path = str(SCRIPTS_DIR / "loop_modeling.py")

    return job_manager.submit_job(
        script_path=script_path,
        args={
            "input": input_file,
            "loop_start": loop_start,
            "loop_end": loop_end,
            "trajectories": trajectories
        },
        job_name=job_name or f"loop_{loop_start}_{loop_end}_{Path(input_file).stem}"
    )

@mcp.tool()
def submit_ligand_docking(
    input_file: str,
    ligand_chain: Optional[str] = None,
    trajectories: int = 10,
    perturbation_cycles: int = 8,
    repack_sidechains: bool = True,
    minimize_final: bool = True,
    output_dir: Optional[str] = None,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit protein-ligand docking for background processing.

    This operation may take >10 minutes for complex ligand optimization.
    Returns a job_id for tracking.

    Args:
        input_file: Path to input PDB file with protein-ligand complex
        ligand_chain: Ligand chain ID (auto-detect if None)
        trajectories: Number of docking trajectories (default: 10)
        perturbation_cycles: Ligand perturbation cycles (default: 8)
        repack_sidechains: Repack sidechains around ligand (default: True)
        minimize_final: Minimize final complex (default: True)
        output_dir: Directory for outputs
        job_name: Optional name for tracking

    Returns:
        Dictionary with job_id for tracking the ligand docking job
    """
    script_path = str(SCRIPTS_DIR / "ligand_docking.py")

    args = {
        "input": input_file,
        "trajectories": trajectories
    }

    if ligand_chain:
        args["ligand_chain"] = ligand_chain

    # Convert boolean flags to script's expected format
    if not minimize_final:
        args["no-minimize"] = True
    if not repack_sidechains:
        args["no-repack"] = True

    # Note: perturbation_cycles and output_dir are controlled by script's default config

    return job_manager.submit_job(
        script_path=script_path,
        args=args,
        job_name=job_name or f"ligand_docking_{Path(input_file).stem}"
    )

@mcp.tool()
def submit_large_refinement(
    input_file: str,
    trajectories: int = 50,
    cycles: int = 500,
    temperature: float = 2.0,
    output_dir: Optional[str] = None,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit large-scale protein refinement for background processing.

    Use this for generating many trajectories (>20) which may take longer.

    Args:
        input_file: Path to input PDB file
        trajectories: Number of refinement trajectories (large number)
        cycles: Number of Monte Carlo cycles per trajectory
        temperature: Monte Carlo temperature for acceptance
        output_dir: Directory for outputs
        job_name: Optional name for tracking

    Returns:
        Dictionary with job_id for tracking the refinement job
    """
    script_path = str(SCRIPTS_DIR / "protein_refinement.py")

    return job_manager.submit_job(
        script_path=script_path,
        args={
            "input": input_file,
            "trajectories": trajectories,
            "cycles": cycles,
            "temperature": temperature
        },
        job_name=job_name or f"large_refinement_{Path(input_file).stem}"
    )

# ==============================================================================
# Batch Processing Tools
# ==============================================================================

@mcp.tool()
def submit_batch_refinement(
    input_files: List[str],
    trajectories: int = 5,
    cycles: int = 100,
    temperature: float = 2.0,
    output_dir: Optional[str] = None,
    job_name: Optional[str] = None
) -> dict:
    """
    Submit batch refinement for multiple PDB files.

    This operation may take >10 minutes for large batches. Returns a job_id for tracking.

    Args:
        input_files: List of PDB file paths to refine
        trajectories: Number of trajectories per structure
        cycles: Number of cycles per trajectory
        temperature: Monte Carlo temperature
        output_dir: Directory for outputs
        job_name: Optional name for tracking

    Returns:
        Dictionary with job_id. Use get_job_status(job_id) to check progress
    """
    # For batch processing, we'll create a wrapper script that handles multiple inputs
    # For now, we'll process them sequentially
    if len(input_files) == 1:
        # Call the submit_large_refinement function directly
        script_path = str(SCRIPTS_DIR / "protein_refinement.py")

        return job_manager.submit_job(
            script_path=script_path,
            args={
                "input": input_files[0],
                "trajectories": trajectories,
                "cycles": cycles,
                "temperature": temperature,
                "output_dir": output_dir
            },
            job_name=job_name or f"batch_refinement_{Path(input_files[0]).stem}"
        )

    # TODO: Implement true batch processing
    return {"status": "error", "error": "Batch processing not yet implemented for multiple files"}

# ==============================================================================
# Utility Tools
# ==============================================================================

@mcp.tool()
def validate_pdb_structure(input_file: str) -> dict:
    """
    Validate a PDB structure for Rosetta compatibility.

    Quick validation of PDB file format and basic structure checks.

    Args:
        input_file: Path to PDB file to validate

    Returns:
        Dictionary with validation results and structure info
    """
    try:
        from lib.io_utils import validate_pdb_file

        # Basic file existence check
        pdb_path = Path(input_file)
        if not pdb_path.exists():
            return {"status": "error", "error": f"File not found: {input_file}"}

        # Validate PDB format
        is_valid = validate_pdb_file(input_file)

        # Get basic structure info
        with open(input_file, 'r') as f:
            lines = f.readlines()

        atom_count = sum(1 for line in lines if line.startswith('ATOM'))
        hetatm_count = sum(1 for line in lines if line.startswith('HETATM'))
        chains = set(line[21] for line in lines if line.startswith('ATOM') and len(line) > 21)

        return {
            "status": "success",
            "valid": is_valid,
            "file_path": str(pdb_path.resolve()),
            "file_size": pdb_path.stat().st_size,
            "atom_count": atom_count,
            "hetatm_count": hetatm_count,
            "chain_count": len(chains),
            "chains": list(sorted(chains))
        }
    except Exception as e:
        return {"status": "error", "error": str(e)}

@mcp.tool()
def list_example_structures() -> dict:
    """
    List available example PDB structures for testing.

    Returns paths to example structures included with the tools.

    Returns:
        Dictionary with example structure paths and descriptions
    """
    examples_dir = MCP_ROOT / "examples" / "data"
    examples = []

    if examples_dir.exists():
        for pdb_file in examples_dir.glob("*.pdb"):
            examples.append({
                "file": str(pdb_file),
                "name": pdb_file.stem,
                "size": pdb_file.stat().st_size if pdb_file.exists() else 0
            })

    # Add demo structures that scripts can create
    demo_structures = [
        {
            "file": "demo_protein.pdb",
            "name": "Demo Protein",
            "description": "Generated by scripts in demo mode",
            "size": "~1KB"
        },
        {
            "file": "demo_complex.pdb",
            "name": "Demo Protein Complex",
            "description": "Two-chain complex for docking",
            "size": "~2KB"
        },
        {
            "file": "demo_ligand_complex.pdb",
            "name": "Demo Protein-Ligand Complex",
            "description": "Protein with small molecule ligand",
            "size": "~1.5KB"
        }
    ]

    return {
        "status": "success",
        "example_files": examples,
        "demo_structures": demo_structures,
        "total_examples": len(examples),
        "note": "Use --demo flag with any script to create demo structures automatically"
    }

# ==============================================================================
# Entry Point
# ==============================================================================

if __name__ == "__main__":
    mcp.run()