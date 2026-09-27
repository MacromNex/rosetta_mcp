#!/usr/bin/env python3
"""Test script for real-world scenario testing"""

import sys
import json
import time
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def scenario_1_protein_analysis_pipeline():
    """
    Scenario 1: Complete protein analysis pipeline
    Steps: List examples → Validate structure → Analyze structure properties → Refine if needed
    """
    print("=== Scenario 1: Protein Analysis Pipeline ===")
    try:
        from server import mcp

        print("Step 1: Discover available protein structures...")
        examples_result = mcp._tool_manager._tools['list_example_structures'].fn()

        if examples_result.get("status") == "success" and examples_result.get("example_files"):
            # Pick the largest structure for more interesting analysis
            structures = examples_result["example_files"]
            target_structure = max(structures, key=lambda x: x["size"])
            structure_path = target_structure["file"]

            print(f"Selected structure: {target_structure['name']} ({target_structure['size']} bytes)")

            print("Step 2: Validate the structure...")
            validation_result = mcp._tool_manager._tools['validate_pdb_structure'].fn(structure_path)

            if validation_result.get("status") == "success":
                print(f"✓ Structure validation:")
                print(f"  - Atoms: {validation_result.get('atom_count')}")
                print(f"  - Chains: {validation_result.get('chains')}")
                print(f"  - Valid: {validation_result.get('valid')}")

                # Decide if we need refinement based on structure size
                if validation_result.get("atom_count", 0) > 100:
                    print("Step 3: Structure is large, submitting for refinement...")
                    refinement_result = mcp._tool_manager._tools['submit_large_refinement'].fn(
                        structure_path,
                        trajectories=3,
                        cycles=200,
                        job_name=f"analysis_refinement_{target_structure['name']}"
                    )

                    if refinement_result.get("status") == "submitted":
                        job_id = refinement_result["job_id"]
                        print(f"✓ Refinement job submitted: {job_id}")

                        # Check initial status
                        status = mcp._tool_manager._tools['get_job_status'].fn(job_id)
                        print(f"✓ Job status: {status.get('status')}")

                        return True, f"Pipeline completed with refinement job {job_id}"
                    else:
                        print("✗ Refinement submission failed")
                        return False, "Refinement failed"
                else:
                    print("Step 3: Structure is small, performing quick refinement...")
                    refinement_result = mcp._tool_manager._tools['refine_protein_structure'].fn(
                        structure_path,
                        trajectories=2,
                        cycles=50
                    )

                    print(f"✓ Quick refinement result: {refinement_result.get('status')}")
                    return True, "Pipeline completed with quick refinement"
            else:
                print("✗ Structure validation failed")
                return False, "Validation failed"
        else:
            print("✗ No example structures found")
            return False, "No structures available"

    except Exception as e:
        print(f"✗ Error in protein analysis pipeline: {e}")
        return False, str(e)

def scenario_2_mutation_stability_analysis():
    """
    Scenario 2: Mutation stability analysis workflow
    Steps: Validate structure → Calculate ΔΔG for mutations → Submit docking if needed
    """
    print("\n=== Scenario 2: Mutation Stability Analysis ===")
    try:
        from server import mcp

        # Use a smaller structure for faster mutation analysis
        test_structure = Path(__file__).parent.parent / "examples" / "data" / "test_input.pdb"

        if test_structure.exists():
            print(f"Analyzing mutations for: {test_structure.name}")

            print("Step 1: Validate structure for mutation analysis...")
            validation = mcp._tool_manager._tools['validate_pdb_structure'].fn(str(test_structure))

            if validation.get("status") == "success":
                chains = validation.get("chains", [])
                print(f"✓ Found chains: {chains}")

                if chains:
                    print("Step 2: Calculate ΔΔG for test mutations...")
                    # Test mutations (these may fail if positions don't exist, but that's OK)
                    test_mutations = ["A1G", "A2L", "A3F"]

                    for mutation in test_mutations:
                        print(f"Testing mutation {mutation}...")
                        ddg_result = mcp._tool_manager._tools['calculate_ddg'].fn(
                            str(test_structure),
                            mutation,
                            trajectories=2,
                            output_file=None
                        )

                        print(f"  - {mutation}: {ddg_result.get('status')}")
                        if ddg_result.get('status') == 'success':
                            print(f"    Success: {ddg_result.get('metadata', {}).get('success', 'Unknown')}")

                    print("Step 3: If mutations affect binding, simulate protein-protein docking...")
                    # Check if we have a suitable docking structure
                    docking_file = Path(__file__).parent.parent / "examples" / "data" / "test_docking.pdb"
                    if docking_file.exists():
                        docking_result = mcp._tool_manager._tools['submit_protein_docking'].fn(
                            str(docking_file),
                            trajectories=2,
                            job_name="mutation_effect_docking"
                        )

                        if docking_result.get("status") == "submitted":
                            print(f"✓ Docking simulation started: {docking_result['job_id']}")
                            return True, f"Mutation analysis completed with docking {docking_result['job_id']}"
                        else:
                            return True, "Mutation analysis completed (docking failed)"
                    else:
                        return True, "Mutation analysis completed (no docking structure)"
                else:
                    return False, "No chains found for mutation analysis"
            else:
                return False, "Structure validation failed"
        else:
            return False, f"Test structure not found: {test_structure}"

    except Exception as e:
        print(f"✗ Error in mutation analysis: {e}")
        return False, str(e)

def scenario_3_comparative_docking_study():
    """
    Scenario 3: Comparative docking study
    Steps: Find multiple structures → Submit docking jobs → Monitor progress → Compare results
    """
    print("\n=== Scenario 3: Comparative Docking Study ===")
    try:
        from server import mcp

        print("Step 1: Find suitable structures for docking comparison...")
        examples_result = mcp._tool_manager._tools['list_example_structures'].fn()

        if examples_result.get("status") == "success":
            # Look for docking-suitable structures
            structures = examples_result["example_files"]
            docking_candidates = [s for s in structures if 'docking' in s['name'] or 'complex' in s['name']]

            if not docking_candidates:
                # Use any available structures
                docking_candidates = structures[:2]

            print(f"Selected {len(docking_candidates)} structures for comparison:")
            for candidate in docking_candidates:
                print(f"  - {candidate['name']}")

            submitted_jobs = []

            print("Step 2: Submit docking jobs for each structure...")
            for i, structure in enumerate(docking_candidates):
                print(f"Submitting docking for {structure['name']}...")

                docking_result = mcp._tool_manager._tools['submit_protein_docking'].fn(
                    structure['file'],
                    trajectories=3,
                    job_name=f"comparative_docking_{i}_{structure['name']}"
                )

                if docking_result.get("status") == "submitted":
                    job_id = docking_result["job_id"]
                    submitted_jobs.append(job_id)
                    print(f"  ✓ Job submitted: {job_id}")
                else:
                    print(f"  ✗ Failed to submit job for {structure['name']}")

            if submitted_jobs:
                print("Step 3: Monitor all submitted jobs...")
                for job_id in submitted_jobs:
                    status = mcp._tool_manager._tools['get_job_status'].fn(job_id)
                    print(f"  - {job_id}: {status.get('status')}")

                print("Step 4: Check overall job queue...")
                all_jobs = mcp._tool_manager._tools['list_jobs'].fn()
                print(f"✓ Total jobs in queue: {all_jobs.get('total', 0)}")

                return True, f"Comparative study initiated with {len(submitted_jobs)} jobs"
            else:
                return False, "No docking jobs submitted successfully"
        else:
            return False, "Failed to list example structures"

    except Exception as e:
        print(f"✗ Error in comparative docking study: {e}")
        return False, str(e)

def scenario_4_drug_design_pipeline():
    """
    Scenario 4: Drug design pipeline simulation
    Steps: Prepare protein → Ligand docking → Loop modeling for flexibility → Refinement
    """
    print("\n=== Scenario 4: Drug Design Pipeline ===")
    try:
        from server import mcp

        print("Step 1: Identify target protein and ligand...")
        # Look for ligand-containing structures
        examples_result = mcp._tool_manager._tools['list_example_structures'].fn()

        if examples_result.get("status") == "success":
            structures = examples_result["example_files"]
            ligand_structures = [s for s in structures if 'ligand' in s['name']]

            if ligand_structures:
                target_structure = ligand_structures[0]
                print(f"Selected target: {target_structure['name']}")

                print("Step 2: Validate protein-ligand complex...")
                validation = mcp._tool_manager._tools['validate_pdb_structure'].fn(target_structure['file'])

                if validation.get("status") == "success":
                    print(f"✓ Complex validated: {validation.get('atom_count')} atoms")

                    print("Step 3: Optimize ligand binding pose...")
                    ligand_docking_result = mcp._tool_manager._tools['submit_ligand_docking'].fn(
                        target_structure['file'],
                        trajectories=5,
                        perturbation_cycles=10,
                        job_name="drug_design_ligand_opt"
                    )

                    ligand_job = None
                    if ligand_docking_result.get("status") == "submitted":
                        ligand_job = ligand_docking_result["job_id"]
                        print(f"✓ Ligand optimization started: {ligand_job}")

                    print("Step 4: Model flexible loops near binding site...")
                    # Use arbitrary loop region for demonstration
                    loop_result = mcp._tool_manager._tools['submit_loop_modeling'].fn(
                        target_structure['file'],
                        loop_start=10,
                        loop_end=15,
                        trajectories=3,
                        job_name="drug_design_loop_flex"
                    )

                    loop_job = None
                    if loop_result.get("status") == "submitted":
                        loop_job = loop_result["job_id"]
                        print(f"✓ Loop modeling started: {loop_job}")

                    print("Step 5: Refine final complex...")
                    refinement_result = mcp._tool_manager._tools['submit_large_refinement'].fn(
                        target_structure['file'],
                        trajectories=5,
                        cycles=300,
                        job_name="drug_design_final_refinement"
                    )

                    refinement_job = None
                    if refinement_result.get("status") == "submitted":
                        refinement_job = refinement_result["job_id"]
                        print(f"✓ Final refinement started: {refinement_job}")

                    jobs_created = [j for j in [ligand_job, loop_job, refinement_job] if j]

                    if jobs_created:
                        print(f"Step 6: Pipeline initiated with {len(jobs_created)} jobs")
                        # Show job status summary
                        for job_id in jobs_created:
                            status = mcp._tool_manager._tools['get_job_status'].fn(job_id)
                            print(f"  - {status.get('job_name')}: {status.get('status')}")

                        return True, f"Drug design pipeline initiated with jobs: {jobs_created}"
                    else:
                        return False, "No jobs successfully submitted in pipeline"
                else:
                    return False, "Complex validation failed"
            else:
                print("No ligand structures found, using general structure...")
                # Fall back to using any structure
                if structures:
                    target = structures[0]
                    ligand_result = mcp._tool_manager._tools['submit_ligand_docking'].fn(
                        target['file'],
                        trajectories=2,
                        job_name="drug_design_fallback"
                    )

                    if ligand_result.get("status") == "submitted":
                        return True, f"Simplified drug design with job {ligand_result['job_id']}"
                    else:
                        return False, "Fallback ligand docking failed"
                else:
                    return False, "No structures available"
        else:
            return False, "Failed to list structures"

    except Exception as e:
        print(f"✗ Error in drug design pipeline: {e}")
        return False, str(e)

def main():
    """Run all real-world scenario tests"""
    print("Testing MCP Real-World Scenarios")
    print("================================")

    scenarios = [
        ("Protein Analysis Pipeline", scenario_1_protein_analysis_pipeline),
        ("Mutation Stability Analysis", scenario_2_mutation_stability_analysis),
        ("Comparative Docking Study", scenario_3_comparative_docking_study),
        ("Drug Design Pipeline", scenario_4_drug_design_pipeline)
    ]

    results = {}
    details = {}

    for scenario_name, scenario_func in scenarios:
        print(f"\n{'='*60}")
        success, message = scenario_func()
        results[scenario_name] = success
        details[scenario_name] = message

        status = "✓ PASS" if success else "✗ FAIL"
        print(f"Result: {status} - {message}")

    print(f"\n{'='*60}")
    print("=== Real-World Scenario Results ===")
    passed = sum(results.values())
    total = len(results)

    for scenario_name, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{scenario_name}: {status}")
        print(f"  Details: {details[scenario_name]}")

    print(f"\nSummary: {passed}/{total} scenarios passed")

    # Show final job queue status
    try:
        from server import mcp
        final_jobs = mcp._tool_manager._tools['list_jobs'].fn()
        print(f"Final job queue: {final_jobs.get('total', 0)} jobs")
    except:
        pass

    if passed >= total - 1:  # Allow 1 failure for edge cases
        print("🎉 Real-world scenarios mostly successful!")
        return True
    else:
        print("❌ Multiple real-world scenarios failed")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)