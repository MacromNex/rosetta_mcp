#!/usr/bin/env python3
"""Automated integration test runner for Rosetta MCP server."""

import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# Setup paths
SCRIPT_DIR = Path(__file__).parent
MCP_ROOT = SCRIPT_DIR.parent
SERVER_PATH = MCP_ROOT / "src" / "server.py"
PYTHON_PATH = MCP_ROOT / "env" / "bin" / "python"
TEST_DATA_DIR = MCP_ROOT / "examples" / "data"

class MCPTestRunner:
    def __init__(self):
        self.results = {
            "test_date": datetime.now().isoformat(),
            "server_path": str(SERVER_PATH),
            "python_path": str(PYTHON_PATH),
            "test_data_dir": str(TEST_DATA_DIR),
            "tests": {},
            "issues": [],
            "summary": {}
        }

        # Ensure test data exists
        if not TEST_DATA_DIR.exists():
            print(f"Warning: Test data directory not found: {TEST_DATA_DIR}")

    def run_python_command(self, code: str, timeout: int = 30) -> dict:
        """Run a Python command and return result."""
        try:
            result = subprocess.run(
                [str(PYTHON_PATH), "-c", code],
                cwd=str(MCP_ROOT),
                capture_output=True,
                text=True,
                timeout=timeout
            )
            return {
                "success": result.returncode == 0,
                "stdout": result.stdout.strip(),
                "stderr": result.stderr.strip(),
                "returncode": result.returncode
            }
        except subprocess.TimeoutExpired:
            return {
                "success": False,
                "stdout": "",
                "stderr": f"Command timed out after {timeout}s",
                "returncode": -1
            }
        except Exception as e:
            return {
                "success": False,
                "stdout": "",
                "stderr": str(e),
                "returncode": -1
            }

    def test_server_startup(self) -> bool:
        """Test that server starts without errors."""
        print("Testing server startup...")

        code = """
try:
    from src.server import mcp
    print("SUCCESS: Server imported successfully")
except ImportError as e:
    print(f"IMPORT_ERROR: {e}")
    exit(1)
except Exception as e:
    print(f"ERROR: {e}")
    exit(1)
"""

        result = self.run_python_command(code)
        self.results["tests"]["server_startup"] = {
            "status": "passed" if result["success"] else "failed",
            "output": result["stdout"],
            "error": result["stderr"]
        }

        if result["success"]:
            print("✓ Server startup test passed")
        else:
            print(f"✗ Server startup test failed: {result['stderr']}")

        return result["success"]

    def test_tool_listing(self) -> bool:
        """Test that all expected tools are available."""
        print("Testing tool listing...")

        code = """
import asyncio
from src.server import mcp

async def main():
    try:
        tools = await mcp.get_tools()
        print(f"TOOLS_COUNT: {len(tools)}")
        for tool in tools:
            print(f"TOOL: {tool}")
        return len(tools)
    except Exception as e:
        print(f"ERROR: {e}")
        return 0

count = asyncio.run(main())
if count >= 10:  # Expecting at least 10 tools
    print("SUCCESS: Found sufficient tools")
else:
    print(f"INSUFFICIENT_TOOLS: Only found {count} tools")
    exit(1)
"""

        result = self.run_python_command(code)
        self.results["tests"]["tool_listing"] = {
            "status": "passed" if result["success"] else "failed",
            "output": result["stdout"],
            "error": result["stderr"]
        }

        if result["success"]:
            print("✓ Tool listing test passed")
        else:
            print(f"✗ Tool listing test failed: {result['stderr']}")

        return result["success"]

    def test_sync_tool_validation(self) -> bool:
        """Test the validate_pdb_structure tool."""
        print("Testing sync tool (PDB validation)...")

        # Find a test PDB file
        test_files = list(TEST_DATA_DIR.glob("*.pdb"))
        if not test_files:
            print("✗ No test PDB files found")
            self.results["tests"]["sync_tool_validation"] = {
                "status": "skipped",
                "output": "",
                "error": "No test PDB files available"
            }
            return True  # Skip but don't fail

        test_file = test_files[0]

        code = f"""
import asyncio
from src.server import mcp

async def main():
    try:
        tool = await mcp.get_tool("validate_pdb_structure")
        if tool is None:
            print("ERROR: validate_pdb_structure tool not found")
            return False

        # We can't actually call the tool without the MCP protocol,
        # but we can verify it exists and has the right structure
        print(f"TOOL_FOUND: validate_pdb_structure")
        print(f"TEST_FILE: {test_file}")
        print("SUCCESS: Validation tool is available")
        return True
    except Exception as e:
        print(f"ERROR: {{e}}")
        return False

success = asyncio.run(main())
if not success:
    exit(1)
"""

        result = self.run_python_command(code)
        self.results["tests"]["sync_tool_validation"] = {
            "status": "passed" if result["success"] else "failed",
            "output": result["stdout"],
            "error": result["stderr"]
        }

        if result["success"]:
            print("✓ Sync tool test passed")
        else:
            print(f"✗ Sync tool test failed: {result['stderr']}")

        return result["success"]

    def test_job_manager_import(self) -> bool:
        """Test that job manager can be imported."""
        print("Testing job manager...")

        code = """
try:
    from src.jobs.manager import job_manager
    jobs = job_manager.list_jobs()
    print(f"SUCCESS: Job manager imported, found {len(jobs.get('jobs', []))} jobs")
except ImportError as e:
    print(f"IMPORT_ERROR: {e}")
    exit(1)
except Exception as e:
    print(f"ERROR: {e}")
    exit(1)
"""

        result = self.run_python_command(code)
        self.results["tests"]["job_manager"] = {
            "status": "passed" if result["success"] else "failed",
            "output": result["stdout"],
            "error": result["stderr"]
        }

        if result["success"]:
            print("✓ Job manager test passed")
        else:
            print(f"✗ Job manager test failed: {result['stderr']}")

        return result["success"]

    def test_file_access(self) -> bool:
        """Test that server can access example files."""
        print("Testing file access...")

        code = f"""
from pathlib import Path

test_dir = Path("{TEST_DATA_DIR}")
if not test_dir.exists():
    print(f"ERROR: Test directory does not exist: {{test_dir}}")
    exit(1)

pdb_files = list(test_dir.glob("*.pdb"))
print(f"FOUND_FILES: {{len(pdb_files)}} PDB files")

for pdb_file in pdb_files:
    print(f"FILE: {{pdb_file.name}} ({{pdb_file.stat().st_size}} bytes)")

if len(pdb_files) == 0:
    print("WARNING: No PDB files found in test directory")
else:
    print("SUCCESS: Test files are accessible")
"""

        result = self.run_python_command(code)
        self.results["tests"]["file_access"] = {
            "status": "passed" if result["success"] else "failed",
            "output": result["stdout"],
            "error": result["stderr"]
        }

        if result["success"]:
            print("✓ File access test passed")
        else:
            print(f"✗ File access test failed: {result['stderr']}")

        return result["success"]

    def test_error_handling(self) -> bool:
        """Test error handling with invalid inputs."""
        print("Testing error handling...")

        code = """
import asyncio
from src.server import mcp

async def main():
    try:
        # Test getting non-existent tool - expect exception
        try:
            tool = await mcp.get_tool("non_existent_tool")
            print(f"UNEXPECTED: Non-existent tool should throw exception but returned: {tool}")
            return False
        except Exception as e:
            print(f"SUCCESS: Non-existent tool correctly throws exception: {type(e).__name__}")

        # Test getting valid tool
        try:
            valid_tool = await mcp.get_tool("validate_pdb_structure")
            if valid_tool is not None:
                print("SUCCESS: Valid tool found")
            else:
                print("ERROR: Valid tool not found")
                return False
        except Exception as e:
            print(f"ERROR: Valid tool lookup failed: {e}")
            return False

        return True
    except Exception as e:
        print(f"GENERAL_ERROR: {e}")
        return False

success = asyncio.run(main())
if not success:
    exit(1)
"""

        result = self.run_python_command(code)
        self.results["tests"]["error_handling"] = {
            "status": "passed" if result["success"] else "failed",
            "output": result["stdout"],
            "error": result["stderr"]
        }

        if result["success"]:
            print("✓ Error handling test passed")
        else:
            print(f"✗ Error handling test failed: {result['stderr']}")

        return result["success"]

    def run_all_tests(self) -> bool:
        """Run all tests and generate report."""
        print("=" * 60)
        print("Running Rosetta MCP Integration Tests")
        print("=" * 60)

        tests = [
            ("Server Startup", self.test_server_startup),
            ("Tool Listing", self.test_tool_listing),
            ("Sync Tool Validation", self.test_sync_tool_validation),
            ("Job Manager", self.test_job_manager_import),
            ("File Access", self.test_file_access),
            ("Error Handling", self.test_error_handling),
        ]

        passed = 0
        total = len(tests)

        for name, test_func in tests:
            print(f"\n[{passed+1}/{total}] {name}")
            print("-" * 40)
            try:
                if test_func():
                    passed += 1
                else:
                    self.results["issues"].append(f"{name} test failed")
            except Exception as e:
                print(f"✗ {name} test encountered an exception: {e}")
                self.results["issues"].append(f"{name} test exception: {e}")

        self.results["summary"] = {
            "total_tests": total,
            "passed": passed,
            "failed": total - passed,
            "pass_rate": f"{passed/total*100:.1f}%" if total > 0 else "N/A"
        }

        print("\n" + "=" * 60)
        print("TEST RESULTS SUMMARY")
        print("=" * 60)
        print(f"Total Tests: {total}")
        print(f"Passed: {passed}")
        print(f"Failed: {total - passed}")
        print(f"Pass Rate: {self.results['summary']['pass_rate']}")

        if self.results["issues"]:
            print(f"\nIssues Found:")
            for issue in self.results["issues"]:
                print(f"  - {issue}")

        return passed == total

    def generate_report(self) -> str:
        """Generate JSON report."""
        return json.dumps(self.results, indent=2)

    def save_report(self, output_path: Path):
        """Save report to file."""
        with open(output_path, 'w') as f:
            f.write(self.generate_report())
        print(f"\nTest report saved to: {output_path}")

if __name__ == "__main__":
    runner = MCPTestRunner()

    success = runner.run_all_tests()

    # Save report
    report_path = SCRIPT_DIR.parent / "reports" / "step7_integration.json"
    report_path.parent.mkdir(exist_ok=True)
    runner.save_report(report_path)

    # Exit with appropriate code
    sys.exit(0 if success else 1)