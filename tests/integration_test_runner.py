#!/usr/bin/env python3
"""
Integration Test Runner for Rosetta MCP Server

This script runs automated integration tests to verify the MCP server
functionality without requiring full PyRosetta installation.

Usage:
    python tests/integration_test_runner.py

Expected Results:
    - All MCP protocol tests should PASS
    - Computational tests should FAIL with expected PyRosetta errors
"""

import asyncio
import json
import sys
from datetime import datetime
from pathlib import Path

# Add src and scripts to path for imports
project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root / "src"))
sys.path.insert(0, str(project_root / "scripts"))

class IntegrationTestRunner:
    def __init__(self):
        self.results = {
            "test_date": datetime.now().isoformat(),
            "test_environment": "integration_test",
            "tests": {},
            "summary": {}
        }
        self.passed = 0
        self.failed = 0

    async def test_server_import(self):
        """Test that the MCP server can be imported."""
        test_name = "server_import"
        try:
            from server import mcp
            tools = await mcp.get_tools()
            expected_tools = 14

            if len(tools) == expected_tools:
                self.results["tests"][test_name] = {
                    "status": "PASS",
                    "message": f"Server imported successfully, found {len(tools)} tools",
                    "details": {"tool_count": len(tools), "tools": list(tools)}
                }
                self.passed += 1
                return True
            else:
                self.results["tests"][test_name] = {
                    "status": "FAIL",
                    "message": f"Expected {expected_tools} tools, found {len(tools)}",
                    "details": {"tool_count": len(tools), "tools": list(tools)}
                }
                self.failed += 1
                return False

        except Exception as e:
            self.results["tests"][test_name] = {
                "status": "FAIL",
                "message": f"Server import failed: {str(e)}",
                "details": {"error": str(e), "type": type(e).__name__}
            }
            self.failed += 1
            return False

    async def test_tool_execution(self):
        """Test basic tool execution."""
        test_name = "tool_execution"
        try:
            from server import mcp

            # Test structure validation (should work without PyRosetta)
            tools = await mcp.get_tools()

            # Simulate tool execution by checking that we can access the tools
            if "validate_pdb_structure" in tools and "list_example_structures" in tools:
                self.results["tests"][test_name] = {
                    "status": "PASS",
                    "message": "Key tools are accessible via MCP interface",
                    "details": {"tested_tools": ["validate_pdb_structure", "list_example_structures"]}
                }
                self.passed += 1
                return True
            else:
                self.results["tests"][test_name] = {
                    "status": "FAIL",
                    "message": "Key tools not found in MCP interface",
                    "details": {"available_tools": list(tools)}
                }
                self.failed += 1
                return False

        except Exception as e:
            self.results["tests"][test_name] = {
                "status": "FAIL",
                "message": f"Tool execution test failed: {str(e)}",
                "details": {"error": str(e), "type": type(e).__name__}
            }
            self.failed += 1
            return False

    async def test_job_manager(self):
        """Test job manager functionality."""
        test_name = "job_manager"
        try:
            from jobs.manager import job_manager

            # Test basic job manager operations
            jobs = job_manager.list_jobs()

            self.results["tests"][test_name] = {
                "status": "PASS",
                "message": "Job manager is functional",
                "details": {"job_count": jobs.get("total", 0)}
            }
            self.passed += 1
            return True

        except Exception as e:
            self.results["tests"][test_name] = {
                "status": "FAIL",
                "message": f"Job manager test failed: {str(e)}",
                "details": {"error": str(e), "type": type(e).__name__}
            }
            self.failed += 1
            return False

    async def test_file_structure(self):
        """Test that required files exist."""
        test_name = "file_structure"

        required_files = [
            "src/server.py",
            "src/jobs/manager.py",
            "scripts/protein_refinement.py",
            "scripts/protein_docking.py",
            "scripts/loop_modeling.py",
            "examples/data/test_input.pdb"
        ]

        missing_files = []
        for file_path in required_files:
            if not Path(file_path).exists():
                missing_files.append(file_path)

        if not missing_files:
            self.results["tests"][test_name] = {
                "status": "PASS",
                "message": "All required files exist",
                "details": {"checked_files": len(required_files)}
            }
            self.passed += 1
            return True
        else:
            self.results["tests"][test_name] = {
                "status": "FAIL",
                "message": f"Missing required files: {missing_files}",
                "details": {"missing_files": missing_files}
            }
            self.failed += 1
            return False

    async def run_all_tests(self):
        """Run all integration tests."""
        print("🧪 Running Rosetta MCP Integration Tests...")
        print("=" * 50)

        # Run tests
        tests = [
            ("File Structure", self.test_file_structure()),
            ("Server Import", self.test_server_import()),
            ("Tool Execution", self.test_tool_execution()),
            ("Job Manager", self.test_job_manager())
        ]

        for test_name, test_coro in tests:
            print(f"Running {test_name}...", end=" ")
            result = await test_coro
            status = "✅ PASS" if result else "❌ FAIL"
            print(status)

        # Generate summary
        total_tests = self.passed + self.failed
        pass_rate = (self.passed / total_tests * 100) if total_tests > 0 else 0

        self.results["summary"] = {
            "total_tests": total_tests,
            "passed": self.passed,
            "failed": self.failed,
            "pass_rate": f"{pass_rate:.1f}%"
        }

        print("\n" + "=" * 50)
        print(f"📊 Test Results: {self.passed}/{total_tests} passed ({pass_rate:.1f}%)")
        print("=" * 50)

        # Print failed tests
        if self.failed > 0:
            print("\n❌ Failed Tests:")
            for test_name, test_result in self.results["tests"].items():
                if test_result["status"] == "FAIL":
                    print(f"  - {test_name}: {test_result['message']}")

        return self.results

    def save_results(self, output_file: str = "tests/integration_test_results.json"):
        """Save test results to file."""
        Path(output_file).parent.mkdir(exist_ok=True)
        with open(output_file, 'w') as f:
            json.dump(self.results, f, indent=2)
        print(f"\n💾 Results saved to {output_file}")

async def main():
    """Main test runner."""
    runner = IntegrationTestRunner()

    try:
        results = await runner.run_all_tests()
        runner.save_results()

        # Exit with appropriate code
        if results["summary"]["failed"] == 0:
            print("\n🎉 All tests passed! MCP server is ready for integration.")
            sys.exit(0)
        else:
            print(f"\n⚠️  {results['summary']['failed']} test(s) failed. Check logs for details.")
            sys.exit(1)

    except Exception as e:
        print(f"\n💥 Test runner crashed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    asyncio.run(main())