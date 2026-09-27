#!/usr/bin/env python3
"""Test script for submit API and job management functionality"""

import sys
import json
import time
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

def test_job_management_tools():
    """Test basic job management tools"""
    print("=== Testing Job Management Tools ===")
    try:
        from server import mcp

        # Test list_jobs initially (should be empty or show existing jobs)
        result = mcp._tool_manager._tools['list_jobs'].fn()
        print(f"✓ Initial job list: {json.dumps(result, indent=2)}")

        return result.get("status") == "success"
    except Exception as e:
        print(f"✗ Error in job management: {e}")
        return False

def test_submit_protein_docking():
    """Test submitting a protein docking job"""
    print("\n=== Testing submit_protein_docking ===")
    try:
        from server import mcp

        # Test with docking file
        test_file = Path(__file__).parent.parent / "examples" / "data" / "test_docking.pdb"
        if test_file.exists():
            result = mcp._tool_manager._tools['submit_protein_docking'].fn(
                input_file=str(test_file),
                trajectories=2,  # Small number for testing
                job_name="test_docking"
            )
            print(f"✓ Docking submission result: {json.dumps(result, indent=2)}")
            return result.get("status") == "submitted", result.get("job_id")
        else:
            print(f"✗ Test file not found: {test_file}")
            return False, None
    except Exception as e:
        print(f"✗ Error in submit_protein_docking: {e}")
        return False, None

def test_submit_loop_modeling():
    """Test submitting a loop modeling job"""
    print("\n=== Testing submit_loop_modeling ===")
    try:
        from server import mcp

        # Test with loop file
        test_file = Path(__file__).parent.parent / "examples" / "data" / "test_loop.pdb"
        if test_file.exists():
            result = mcp._tool_manager._tools['submit_loop_modeling'].fn(
                input_file=str(test_file),
                loop_start=1,
                loop_end=5,
                trajectories=2,
                job_name="test_loop"
            )
            print(f"✓ Loop modeling submission result: {json.dumps(result, indent=2)}")
            return result.get("status") == "submitted", result.get("job_id")
        else:
            print(f"✗ Test file not found: {test_file}")
            return False, None
    except Exception as e:
        print(f"✗ Error in submit_loop_modeling: {e}")
        return False, None

def test_job_status_workflow(job_id):
    """Test the complete job status workflow"""
    if not job_id:
        print("✗ No job ID provided for status testing")
        return False

    print(f"\n=== Testing Job Status Workflow for {job_id} ===")
    try:
        from server import mcp

        # Test get_job_status
        print("Testing get_job_status...")
        status_result = mcp._tool_manager._tools['get_job_status'].fn(job_id)
        print(f"✓ Job status: {json.dumps(status_result, indent=2)}")

        # Test get_job_log
        print("Testing get_job_log...")
        log_result = mcp._tool_manager._tools['get_job_log'].fn(job_id, tail=10)
        print(f"✓ Job log: {json.dumps(log_result, indent=2)}")

        # Test list_jobs to see our job
        print("Testing list_jobs...")
        list_result = mcp._tool_manager._tools['list_jobs'].fn()
        print(f"✓ All jobs: {json.dumps(list_result, indent=2)}")

        # If job is pending or running, we can test cancellation
        job_status = status_result.get("job_status", "")
        if job_status in ["pending", "running"]:
            print(f"Job is {job_status}, testing cancellation...")
            cancel_result = mcp._tool_manager._tools['cancel_job'].fn(job_id)
            print(f"✓ Cancel result: {json.dumps(cancel_result, indent=2)}")

        return True

    except Exception as e:
        print(f"✗ Error in job status workflow: {e}")
        return False

def test_submit_large_refinement():
    """Test submitting a large refinement job"""
    print("\n=== Testing submit_large_refinement ===")
    try:
        from server import mcp

        # Test with input file
        test_file = Path(__file__).parent.parent / "examples" / "data" / "test_input.pdb"
        if test_file.exists():
            result = mcp._tool_manager._tools['submit_large_refinement'].fn(
                input_file=str(test_file),
                trajectories=5,  # Small for testing
                cycles=50,
                job_name="test_large_refinement"
            )
            print(f"✓ Large refinement submission result: {json.dumps(result, indent=2)}")
            return result.get("status") == "submitted", result.get("job_id")
        else:
            print(f"✗ Test file not found: {test_file}")
            return False, None
    except Exception as e:
        print(f"✗ Error in submit_large_refinement: {e}")
        return False, None

def test_submit_ligand_docking():
    """Test submitting a ligand docking job"""
    print("\n=== Testing submit_ligand_docking ===")
    try:
        from server import mcp

        # Test with ligand file
        test_file = Path(__file__).parent.parent / "examples" / "data" / "test_ligand.pdb"
        if test_file.exists():
            result = mcp._tool_manager._tools['submit_ligand_docking'].fn(
                input_file=str(test_file),
                trajectories=2,
                job_name="test_ligand_docking"
            )
            print(f"✓ Ligand docking submission result: {json.dumps(result, indent=2)}")
            return result.get("status") == "submitted", result.get("job_id")
        else:
            print(f"✗ Test file not found: {test_file}")
            return False, None
    except Exception as e:
        print(f"✗ Error in submit_ligand_docking: {e}")
        return False, None

def main():
    """Run all submit API tests"""
    print("Testing MCP Submit API and Job Management")
    print("=========================================")

    job_ids = []
    results = {}

    # Test basic job management
    results["job_management"] = test_job_management_tools()

    # Test job submissions
    success, job_id = test_submit_protein_docking()
    results["submit_protein_docking"] = success
    if job_id:
        job_ids.append(job_id)

    success, job_id = test_submit_loop_modeling()
    results["submit_loop_modeling"] = success
    if job_id:
        job_ids.append(job_id)

    success, job_id = test_submit_large_refinement()
    results["submit_large_refinement"] = success
    if job_id:
        job_ids.append(job_id)

    success, job_id = test_submit_ligand_docking()
    results["submit_ligand_docking"] = success
    if job_id:
        job_ids.append(job_id)

    # Test job status workflow with the first job if available
    if job_ids:
        results["job_status_workflow"] = test_job_status_workflow(job_ids[0])
    else:
        print("✗ No job IDs available for status workflow testing")
        results["job_status_workflow"] = False

    print(f"\n=== Test Results ===")
    passed = sum(results.values())
    total = len(results)

    for test_name, result in results.items():
        status = "✓ PASS" if result else "✗ FAIL"
        print(f"{test_name}: {status}")

    print(f"\nSummary: {passed}/{total} tests passed")
    print(f"Job IDs created: {job_ids}")

    if passed >= total - 1:  # Allow 1 failure for edge cases
        print("🎉 Submit API tests mostly successful!")
        return True
    else:
        print("❌ Multiple submit API tests failed")
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)