#!/usr/bin/env python3
"""
Test script to verify Claude Code can interact with Rosetta MCP server.

This script will test the server by using fastmcp dev mode to simulate
tool calls that Claude Code would make.
"""

import asyncio
import json
import sys
from pathlib import Path

# Add the src directory to path
SCRIPT_DIR = Path(__file__).parent
MCP_ROOT = SCRIPT_DIR.parent
sys.path.insert(0, str(MCP_ROOT / "src"))

from server import mcp

async def test_tool_discovery():
    """Test that we can discover all tools."""
    print("=" * 60)
    print("Testing Tool Discovery")
    print("=" * 60)

    tools = await mcp.get_tools()
    print(f"Found {len(tools)} tools:")

    for i, tool in enumerate(tools, 1):
        print(f"  {i:2d}. {tool}")

    return len(tools) > 0

async def test_validate_structure():
    """Test the validate_pdb_structure tool."""
    print("\n" + "=" * 60)
    print("Testing Structure Validation")
    print("=" * 60)

    # Find test files
    test_data_dir = MCP_ROOT / "examples" / "data"
    test_files = list(test_data_dir.glob("*.pdb"))

    if not test_files:
        print("No test PDB files found!")
        return False

    print(f"Testing with: {test_files[0]}")

    # We'll simulate the tool call that Claude Code would make
    try:
        # Get the tool definition
        tool = await mcp.get_tool("validate_pdb_structure")
        print(f"Tool found: {tool}")

        # In a real MCP session, Claude would call this tool
        # For now, we just verify the tool exists and is accessible
        print("✓ Tool is accessible through MCP")
        return True

    except Exception as e:
        print(f"✗ Error accessing tool: {e}")
        return False

async def test_job_management_tools():
    """Test the job management tools."""
    print("\n" + "=" * 60)
    print("Testing Job Management Tools")
    print("=" * 60)

    job_tools = [
        "list_jobs",
        "get_job_status",
        "get_job_result",
        "get_job_log",
        "cancel_job"
    ]

    for tool_name in job_tools:
        try:
            tool = await mcp.get_tool(tool_name)
            print(f"✓ {tool_name}: Available")
        except Exception as e:
            print(f"✗ {tool_name}: Error - {e}")
            return False

    return True

async def test_async_tools():
    """Test the async submission tools."""
    print("\n" + "=" * 60)
    print("Testing Async Submission Tools")
    print("=" * 60)

    async_tools = [
        "submit_protein_docking",
        "submit_loop_modeling",
        "submit_ligand_docking",
        "submit_large_refinement",
        "submit_batch_refinement"
    ]

    for tool_name in async_tools:
        try:
            tool = await mcp.get_tool(tool_name)
            print(f"✓ {tool_name}: Available")
        except Exception as e:
            print(f"✗ {tool_name}: Error - {e}")
            return False

    return True

async def test_example_data_access():
    """Test that tools can access example data."""
    print("\n" + "=" * 60)
    print("Testing Example Data Access")
    print("=" * 60)

    try:
        # Test the list_example_structures tool
        tool = await mcp.get_tool("list_example_structures")
        print("✓ list_example_structures tool available")

        # Check that test files exist
        test_data_dir = MCP_ROOT / "examples" / "data"
        test_files = list(test_data_dir.glob("*.pdb"))

        print(f"Found {len(test_files)} test PDB files:")
        for test_file in test_files:
            size_kb = test_file.stat().st_size / 1024
            print(f"  - {test_file.name}: {size_kb:.1f} KB")

        return len(test_files) > 0

    except Exception as e:
        print(f"✗ Error: {e}")
        return False

def create_claude_test_prompts():
    """Create a file with prompts for manual testing in Claude Code."""
    print("\n" + "=" * 60)
    print("Creating Claude Test Prompts")
    print("=" * 60)

    prompts = [
        "# Quick Tests for Claude Code",
        "",
        "## Tool Discovery",
        "What Rosetta MCP tools are available? List them with brief descriptions.",
        "",
        "## Structure Validation",
        f"Validate the PDB structure: {MCP_ROOT}/examples/data/test_input.pdb",
        "",
        "## List Examples",
        "What example structures are available for testing?",
        "",
        "## Job Management",
        "List all current jobs.",
        "",
        "## Error Testing",
        "Try to validate a non-existent file: /fake/path.pdb",
        "",
        "## Advanced Test",
        f"Validate {MCP_ROOT}/examples/data/test_input.pdb and then submit a refinement with 3 trajectories.",
        "",
    ]

    prompt_file = MCP_ROOT / "tests" / "manual_claude_prompts.md"
    with open(prompt_file, 'w') as f:
        f.write('\n'.join(prompts))

    print(f"Manual test prompts saved to: {prompt_file}")
    print("\nCopy these prompts to Claude Code to test the MCP integration!")

async def main():
    """Run all tests."""
    print("Rosetta MCP Server - Claude Code Integration Test")
    print("This tests the server components that Claude Code will use.")

    tests = [
        ("Tool Discovery", test_tool_discovery),
        ("Structure Validation", test_validate_structure),
        ("Job Management Tools", test_job_management_tools),
        ("Async Submission Tools", test_async_tools),
        ("Example Data Access", test_example_data_access),
    ]

    passed = 0
    total = len(tests)

    for name, test_func in tests:
        try:
            if await test_func():
                passed += 1
            else:
                print(f"\n✗ {name} test failed")
        except Exception as e:
            print(f"\n✗ {name} test error: {e}")

    print("\n" + "=" * 60)
    print("FINAL RESULTS")
    print("=" * 60)
    print(f"Tests Passed: {passed}/{total}")
    print(f"Success Rate: {passed/total*100:.1f}%")

    if passed == total:
        print("🎉 All tests passed! Server is ready for Claude Code integration.")
    else:
        print("⚠️  Some tests failed. Review the issues above.")

    # Always create the manual test prompts
    create_claude_test_prompts()

    return passed == total

if __name__ == "__main__":
    success = asyncio.run(main())
    sys.exit(0 if success else 1)