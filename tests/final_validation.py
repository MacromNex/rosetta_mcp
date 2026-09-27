#!/usr/bin/env python3
"""
Final validation checklist for Rosetta MCP server integration.

This script validates all critical components required for production deployment.
"""

import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# Setup paths
SCRIPT_DIR = Path(__file__).parent
MCP_ROOT = SCRIPT_DIR.parent
SERVER_PATH = MCP_ROOT / "src" / "server.py"
PYTHON_PATH = MCP_ROOT / "env" / "bin" / "python"

class FinalValidator:
    def __init__(self):
        self.checklist = {
            "timestamp": datetime.now().isoformat(),
            "server_validation": {},
            "claude_code_integration": {},
            "tool_functionality": {},
            "job_management": {},
            "error_handling": {},
            "file_access": {},
            "gemini_cli_integration": {},
            "documentation": {},
            "production_readiness": {}
        }

    def run_command(self, cmd: list, timeout: int = 30) -> dict:
        """Run a shell command and return the result."""
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=str(MCP_ROOT)
            )
            return {
                "success": result.returncode == 0,
                "stdout": result.stdout.strip(),
                "stderr": result.stderr.strip(),
                "returncode": result.returncode
            }
        except subprocess.TimeoutExpired:
            return {"success": False, "error": "Command timed out"}
        except Exception as e:
            return {"success": False, "error": str(e)}

    def validate_server(self) -> bool:
        """Validate server startup and basic functionality."""
        print("🔍 Validating server startup...")

        # Test 1: Server imports without errors
        result = self.run_command([str(PYTHON_PATH), "-c", "from src.server import mcp; print('Server imports OK')"])
        self.checklist["server_validation"]["imports"] = result["success"]

        if not result["success"]:
            print(f"❌ Server import failed: {result.get('stderr', 'Unknown error')}")
            return False

        # Test 2: All tools are available
        result = self.run_command([
            str(PYTHON_PATH), "-c", """
import asyncio
from src.server import mcp
async def main():
    tools = await mcp.get_tools()
    print(f"TOOL_COUNT:{len(tools)}")
    return len(tools) >= 14
print("SUCCESS" if asyncio.run(main()) else "FAIL")
"""
        ])

        tool_count_ok = result["success"] and "SUCCESS" in result["stdout"]
        self.checklist["server_validation"]["tool_count"] = tool_count_ok

        if not tool_count_ok:
            print("❌ Tool count validation failed")
            return False

        # Test 3: Dev server can start
        result = self.run_command(["timeout", "5s", "./env/bin/fastmcp", "dev", "src/server.py"])
        # Dev server timeout is expected, we just want to see it starts
        output = result.get("stdout", "") + " " + result.get("stderr", "") + " " + result.get("error", "")
        dev_server_ok = ("Session token:" in output or
                        "MCP Inspector" in output or
                        "Proxy server listening" in output)
        self.checklist["server_validation"]["dev_server"] = dev_server_ok

        if dev_server_ok:
            print("✅ Server validation passed")
        else:
            print("❌ Dev server validation failed")

        return tool_count_ok and dev_server_ok

    def validate_claude_code_integration(self) -> bool:
        """Validate Claude Code MCP integration."""
        print("🔍 Validating Claude Code integration...")

        # Test 1: Claude CLI is available
        result = self.run_command(["claude", "--version"])
        claude_available = result["success"]
        self.checklist["claude_code_integration"]["cli_available"] = claude_available

        if not claude_available:
            print("❌ Claude CLI not available")
            return False

        # Test 2: Server is registered
        result = self.run_command(["claude", "mcp", "list"])
        server_registered = result["success"] and "rosetta:" in result["stdout"]
        self.checklist["claude_code_integration"]["server_registered"] = server_registered

        if not server_registered:
            print("❌ Rosetta server not registered with Claude Code")
            return False

        # Test 3: Server is connected
        server_connected = "✓ Connected" in result["stdout"]
        self.checklist["claude_code_integration"]["server_connected"] = server_connected

        if server_connected:
            print("✅ Claude Code integration passed")
        else:
            print("❌ Rosetta server registered but not connected")

        return server_connected

    def validate_tool_functionality(self) -> bool:
        """Validate tool functionality and parameters."""
        print("🔍 Validating tool functionality...")

        # Test individual tool access
        tools_to_test = [
            "validate_pdb_structure",
            "refine_protein_structure",
            "calculate_ddg",
            "submit_protein_docking",
            "list_jobs"
        ]

        all_tools_ok = True

        for tool_name in tools_to_test:
            result = self.run_command([
                str(PYTHON_PATH), "-c", f"""
import asyncio
from src.server import mcp
async def main():
    try:
        tool = await mcp.get_tool('{tool_name}')
        return tool is not None
    except:
        return False
print("SUCCESS" if asyncio.run(main()) else "FAIL")
"""
            ])

            tool_ok = result["success"] and "SUCCESS" in result["stdout"]
            self.checklist["tool_functionality"][tool_name] = tool_ok

            if not tool_ok:
                print(f"❌ Tool {tool_name} not accessible")
                all_tools_ok = False

        if all_tools_ok:
            print("✅ Tool functionality validation passed")

        return all_tools_ok

    def validate_job_management(self) -> bool:
        """Validate job management system."""
        print("🔍 Validating job management...")

        # Test 1: Job manager can be imported
        result = self.run_command([
            str(PYTHON_PATH), "-c", """
from src.jobs.manager import job_manager
jobs = job_manager.list_jobs()
print(f"JOB_MANAGER_OK:{len(jobs.get('jobs', []))}")
"""
        ])

        job_manager_ok = result["success"] and "JOB_MANAGER_OK:" in result["stdout"]
        self.checklist["job_management"]["manager_import"] = job_manager_ok

        # Test 2: Jobs directory exists or can be created
        jobs_dir = MCP_ROOT / "jobs"
        jobs_dir_ok = jobs_dir.exists() or self._create_jobs_dir()
        self.checklist["job_management"]["jobs_directory"] = jobs_dir_ok

        success = job_manager_ok and jobs_dir_ok

        if success:
            print("✅ Job management validation passed")
        else:
            print("❌ Job management validation failed")

        return success

    def _create_jobs_dir(self) -> bool:
        """Create jobs directory if it doesn't exist."""
        try:
            jobs_dir = MCP_ROOT / "jobs"
            jobs_dir.mkdir(exist_ok=True)
            return True
        except Exception:
            return False

    def validate_error_handling(self) -> bool:
        """Validate error handling capabilities."""
        print("🔍 Validating error handling...")

        # Test 1: Invalid tool name handling
        result = self.run_command([
            str(PYTHON_PATH), "-c", """
import asyncio
from src.server import mcp
async def main():
    try:
        tool = await mcp.get_tool('nonexistent_tool')
        return False  # Should not reach here
    except Exception:
        return True  # Exception is expected
print("SUCCESS" if asyncio.run(main()) else "FAIL")
"""
        ])

        error_handling_ok = result["success"] and "SUCCESS" in result["stdout"]
        self.checklist["error_handling"]["invalid_tool"] = error_handling_ok

        if error_handling_ok:
            print("✅ Error handling validation passed")
        else:
            print("❌ Error handling validation failed")

        return error_handling_ok

    def validate_file_access(self) -> bool:
        """Validate file access and example data."""
        print("🔍 Validating file access...")

        # Test 1: Example data directory exists
        examples_dir = MCP_ROOT / "examples" / "data"
        examples_exist = examples_dir.exists()
        self.checklist["file_access"]["examples_directory"] = examples_exist

        # Test 2: PDB files are accessible
        pdb_files_count = 0
        if examples_exist:
            pdb_files = list(examples_dir.glob("*.pdb"))
            pdb_files_count = len(pdb_files)

        pdb_files_ok = pdb_files_count > 0
        self.checklist["file_access"]["pdb_files_count"] = pdb_files_count
        self.checklist["file_access"]["pdb_files_accessible"] = pdb_files_ok

        success = examples_exist and pdb_files_ok

        if success:
            print(f"✅ File access validation passed ({pdb_files_count} PDB files)")
        else:
            print("❌ File access validation failed")

        return success

    def validate_gemini_cli_integration(self) -> bool:
        """Validate Gemini CLI integration."""
        print("🔍 Validating Gemini CLI integration...")

        # Test 1: Gemini CLI is available
        result = self.run_command(["gemini", "--version"])
        gemini_available = result["success"]
        self.checklist["gemini_cli_integration"]["cli_available"] = gemini_available

        # Test 2: Configuration exists
        gemini_config = Path.home() / ".gemini" / "settings.json"
        config_exists = gemini_config.exists()
        self.checklist["gemini_cli_integration"]["config_exists"] = config_exists

        # Test 3: Rosetta server configured
        rosetta_configured = False
        if config_exists:
            try:
                with open(gemini_config, 'r') as f:
                    config = json.load(f)
                rosetta_configured = "mcpServers" in config and "rosetta" in config.get("mcpServers", {})
            except Exception:
                pass

        self.checklist["gemini_cli_integration"]["rosetta_configured"] = rosetta_configured

        success = gemini_available and config_exists and rosetta_configured

        if success:
            print("✅ Gemini CLI integration passed")
        elif not gemini_available:
            print("⚠️ Gemini CLI not available (optional)")
        else:
            print("❌ Gemini CLI integration failed")

        # Gemini CLI is optional, so we don't fail the overall validation
        return True

    def validate_documentation(self) -> bool:
        """Validate documentation completeness."""
        print("🔍 Validating documentation...")

        required_files = [
            "README.md",
            "reports/step7_integration.md",
            "tests/test_prompts.md",
            "tests/manual_claude_prompts.md",
            "tests/gemini_test_prompts.md",
            "tests/quick_fixes.md"
        ]

        all_docs_exist = True

        for file_path in required_files:
            full_path = MCP_ROOT / file_path
            exists = full_path.exists()
            self.checklist["documentation"][file_path] = exists

            if not exists:
                print(f"❌ Missing documentation: {file_path}")
                all_docs_exist = False

        if all_docs_exist:
            print("✅ Documentation validation passed")

        return all_docs_exist

    def validate_production_readiness(self) -> bool:
        """Final production readiness check."""
        print("🔍 Validating production readiness...")

        # Check all validation results
        validations = [
            ("Server", self.checklist["server_validation"]),
            ("Claude Code", self.checklist["claude_code_integration"]),
            ("Tools", self.checklist["tool_functionality"]),
            ("Jobs", self.checklist["job_management"]),
            ("Error Handling", self.checklist["error_handling"]),
            ("File Access", self.checklist["file_access"]),
            ("Documentation", self.checklist["documentation"])
        ]

        all_passed = True
        summary = {}

        for name, validation_data in validations:
            if isinstance(validation_data, dict):
                # Count successful validations
                total_checks = len(validation_data)
                passed_checks = sum(1 for v in validation_data.values()
                                  if isinstance(v, bool) and v)

                if isinstance(list(validation_data.values())[0], bool):
                    category_passed = all(v for v in validation_data.values())
                else:
                    # Handle mixed types (like file counts)
                    category_passed = passed_checks == total_checks or passed_checks >= (total_checks * 0.8)

                summary[name] = {
                    "passed": category_passed,
                    "checks": f"{passed_checks}/{total_checks}"
                }

                if not category_passed:
                    all_passed = False
            else:
                # Simple boolean check
                summary[name] = {
                    "passed": validation_data,
                    "checks": "1/1" if validation_data else "0/1"
                }

                if not validation_data:
                    all_passed = False

        self.checklist["production_readiness"]["overall_status"] = all_passed
        self.checklist["production_readiness"]["summary"] = summary

        return all_passed

    def run_final_validation(self) -> bool:
        """Run complete final validation."""
        print("🔧 Running Final Production Validation")
        print("=" * 60)

        validators = [
            ("Server Validation", self.validate_server),
            ("Claude Code Integration", self.validate_claude_code_integration),
            ("Tool Functionality", self.validate_tool_functionality),
            ("Job Management", self.validate_job_management),
            ("Error Handling", self.validate_error_handling),
            ("File Access", self.validate_file_access),
            ("Gemini CLI Integration", self.validate_gemini_cli_integration),
            ("Documentation", self.validate_documentation),
            ("Production Readiness", self.validate_production_readiness)
        ]

        all_passed = True

        for name, validator in validators:
            print(f"\n📋 {name}")
            print("-" * 40)

            try:
                if not validator():
                    all_passed = False
            except Exception as e:
                print(f"💥 {name} validation error: {e}")
                all_passed = False

        # Final summary
        print("\n" + "=" * 60)
        print("🎯 FINAL VALIDATION RESULTS")
        print("=" * 60)

        if all_passed:
            print("🎉 ALL VALIDATIONS PASSED")
            print("✅ Rosetta MCP Server is READY FOR PRODUCTION")
        else:
            print("❌ Some validations failed")
            print("⚠️ Review failed checks before production deployment")

        # Print summary
        if "production_readiness" in self.checklist and "summary" in self.checklist["production_readiness"]:
            print(f"\nValidation Summary:")
            for category, result in self.checklist["production_readiness"]["summary"].items():
                status = "✅" if result["passed"] else "❌"
                print(f"  {status} {category}: {result['checks']}")

        return all_passed

    def save_results(self, output_path: Path):
        """Save validation results to file."""
        with open(output_path, 'w') as f:
            json.dump(self.checklist, f, indent=2)
        print(f"\n💾 Validation results saved to: {output_path}")

def main():
    """Run final validation from command line."""
    import argparse

    parser = argparse.ArgumentParser(description="Run final validation for Rosetta MCP server")
    parser.add_argument("--save", help="Save results to file",
                       default="reports/final_validation.json")

    args = parser.parse_args()

    validator = FinalValidator()
    success = validator.run_final_validation()

    # Save results
    save_path = MCP_ROOT / args.save
    save_path.parent.mkdir(exist_ok=True)
    validator.save_results(save_path)

    sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()