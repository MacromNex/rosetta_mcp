#!/usr/bin/env python3
"""
Debugging utilities for Rosetta MCP server.

This module provides tools for troubleshooting common issues with the MCP server,
job management, and tool execution.
"""

import json
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

# Setup paths
SCRIPT_DIR = Path(__file__).parent
MCP_ROOT = SCRIPT_DIR.parent
SERVER_PATH = MCP_ROOT / "src" / "server.py"
PYTHON_PATH = MCP_ROOT / "env" / "bin" / "python"
JOBS_DIR = MCP_ROOT / "jobs"

class MCPDebugger:
    """Debugging utilities for MCP server issues."""

    def __init__(self):
        self.results = {}

    def check_server_startup(self) -> Dict:
        """Check if server can start up correctly."""
        print("🔍 Checking server startup...")

        result = {"status": "unknown", "issues": [], "output": ""}

        try:
            # Test import
            cmd_result = subprocess.run(
                [str(PYTHON_PATH), "-c", "from src.server import mcp; print('Server imports OK')"],
                cwd=str(MCP_ROOT),
                capture_output=True,
                text=True,
                timeout=30
            )

            if cmd_result.returncode != 0:
                result["status"] = "failed"
                result["issues"].append(f"Import failed: {cmd_result.stderr}")
                result["output"] = cmd_result.stderr
            else:
                result["status"] = "passed"
                result["output"] = cmd_result.stdout
                print("✅ Server startup OK")

        except subprocess.TimeoutExpired:
            result["status"] = "failed"
            result["issues"].append("Server startup timed out")
        except Exception as e:
            result["status"] = "failed"
            result["issues"].append(f"Startup check failed: {e}")

        return result

    def check_tool_availability(self) -> Dict:
        """Check that all expected tools are available."""
        print("🔍 Checking tool availability...")

        result = {"status": "unknown", "issues": [], "tools": [], "expected_count": 14}

        try:
            cmd_result = subprocess.run([
                str(PYTHON_PATH), "-c", """
import asyncio
from src.server import mcp

async def main():
    tools = await mcp.get_tools()
    print(f"TOOL_COUNT:{len(tools)}")
    for tool in tools:
        print(f"TOOL:{tool}")

asyncio.run(main())
"""
            ], cwd=str(MCP_ROOT), capture_output=True, text=True, timeout=30)

            if cmd_result.returncode == 0:
                lines = cmd_result.stdout.strip().split('\n')
                tool_count = 0
                tools = []

                for line in lines:
                    if line.startswith("TOOL_COUNT:"):
                        tool_count = int(line.split(":")[1])
                    elif line.startswith("TOOL:"):
                        tools.append(line.split(":", 1)[1])

                result["tools"] = tools

                if tool_count >= result["expected_count"]:
                    result["status"] = "passed"
                    print(f"✅ Found {tool_count} tools")
                else:
                    result["status"] = "failed"
                    result["issues"].append(f"Expected {result['expected_count']} tools, found {tool_count}")
            else:
                result["status"] = "failed"
                result["issues"].append(f"Tool check failed: {cmd_result.stderr}")

        except Exception as e:
            result["status"] = "failed"
            result["issues"].append(f"Tool availability check failed: {e}")

        return result

    def check_file_access(self) -> Dict:
        """Check access to example files and directories."""
        print("🔍 Checking file access...")

        result = {"status": "unknown", "issues": [], "accessible_files": []}

        # Check critical directories
        critical_dirs = [
            MCP_ROOT / "src",
            MCP_ROOT / "scripts",
            MCP_ROOT / "examples" / "data",
            MCP_ROOT / "jobs"
        ]

        for dir_path in critical_dirs:
            if not dir_path.exists():
                result["issues"].append(f"Missing directory: {dir_path}")
            elif not dir_path.is_dir():
                result["issues"].append(f"Not a directory: {dir_path}")

        # Check example files
        examples_dir = MCP_ROOT / "examples" / "data"
        if examples_dir.exists():
            pdb_files = list(examples_dir.glob("*.pdb"))
            result["accessible_files"] = [str(f) for f in pdb_files]

            if len(pdb_files) == 0:
                result["issues"].append("No PDB files found in examples/data")
            else:
                print(f"✅ Found {len(pdb_files)} example PDB files")
        else:
            result["issues"].append("Examples directory does not exist")

        # Check jobs directory
        if not JOBS_DIR.exists():
            try:
                JOBS_DIR.mkdir(parents=True)
                print("✅ Created jobs directory")
            except Exception as e:
                result["issues"].append(f"Cannot create jobs directory: {e}")

        result["status"] = "passed" if not result["issues"] else "failed"
        return result

    def check_job_manager(self) -> Dict:
        """Check job manager functionality."""
        print("🔍 Checking job manager...")

        result = {"status": "unknown", "issues": [], "job_count": 0}

        try:
            cmd_result = subprocess.run([
                str(PYTHON_PATH), "-c", """
from src.jobs.manager import job_manager
jobs = job_manager.list_jobs()
print(f"JOB_COUNT:{len(jobs.get('jobs', []))}")
print("Job manager working")
"""
            ], cwd=str(MCP_ROOT), capture_output=True, text=True, timeout=30)

            if cmd_result.returncode == 0:
                output_lines = cmd_result.stdout.strip().split('\n')
                for line in output_lines:
                    if line.startswith("JOB_COUNT:"):
                        result["job_count"] = int(line.split(":")[1])

                result["status"] = "passed"
                print(f"✅ Job manager OK, {result['job_count']} jobs")
            else:
                result["status"] = "failed"
                result["issues"].append(f"Job manager check failed: {cmd_result.stderr}")

        except Exception as e:
            result["status"] = "failed"
            result["issues"].append(f"Job manager check failed: {e}")

        return result

    def check_python_environment(self) -> Dict:
        """Check Python environment and dependencies."""
        print("🔍 Checking Python environment...")

        result = {"status": "unknown", "issues": [], "python_version": "", "packages": {}}

        try:
            # Check Python version
            cmd_result = subprocess.run([
                str(PYTHON_PATH), "--version"
            ], capture_output=True, text=True)

            if cmd_result.returncode == 0:
                result["python_version"] = cmd_result.stdout.strip()
                print(f"✅ Python: {result['python_version']}")
            else:
                result["issues"].append("Cannot get Python version")

            # Check required packages
            required_packages = ["fastmcp", "loguru"]
            for package in required_packages:
                try:
                    pkg_result = subprocess.run([
                        str(PYTHON_PATH), "-c", f"import {package}; print('{package}: OK')"
                    ], capture_output=True, text=True, timeout=10)

                    if pkg_result.returncode == 0:
                        result["packages"][package] = "OK"
                        print(f"✅ Package {package}: OK")
                    else:
                        result["packages"][package] = f"Error: {pkg_result.stderr}"
                        result["issues"].append(f"Package {package} not available")
                except Exception as e:
                    result["packages"][package] = f"Exception: {e}"
                    result["issues"].append(f"Cannot check package {package}: {e}")

            result["status"] = "passed" if not result["issues"] else "failed"

        except Exception as e:
            result["status"] = "failed"
            result["issues"].append(f"Environment check failed: {e}")

        return result

    def check_mcp_integration(self) -> Dict:
        """Check MCP client integration (Claude Code)."""
        print("🔍 Checking MCP integration...")

        result = {"status": "unknown", "issues": [], "integrations": {}}

        # Check Claude Code integration
        try:
            claude_result = subprocess.run([
                "claude", "mcp", "list"
            ], capture_output=True, text=True, timeout=30)

            if claude_result.returncode == 0:
                if "rosetta:" in claude_result.stdout and "✓ Connected" in claude_result.stdout:
                    result["integrations"]["claude_code"] = "Connected"
                    print("✅ Claude Code: Connected")
                elif "rosetta:" in claude_result.stdout:
                    result["integrations"]["claude_code"] = "Registered but connection issue"
                    result["issues"].append("Rosetta server registered but not connected in Claude Code")
                else:
                    result["integrations"]["claude_code"] = "Not registered"
                    result["issues"].append("Rosetta server not registered in Claude Code")
            else:
                result["integrations"]["claude_code"] = f"Error: {claude_result.stderr}"
                result["issues"].append("Claude Code MCP check failed")

        except FileNotFoundError:
            result["integrations"]["claude_code"] = "Claude Code not installed"
            result["issues"].append("Claude Code CLI not found")
        except Exception as e:
            result["integrations"]["claude_code"] = f"Exception: {e}"
            result["issues"].append(f"Claude Code check failed: {e}")

        # Check Gemini CLI integration
        try:
            gemini_config = Path.home() / ".gemini" / "settings.json"
            if gemini_config.exists():
                with open(gemini_config, 'r') as f:
                    config = json.load(f)

                if "mcpServers" in config and "rosetta" in config["mcpServers"]:
                    result["integrations"]["gemini_cli"] = "Configured"
                    print("✅ Gemini CLI: Configured")
                else:
                    result["integrations"]["gemini_cli"] = "Not configured"
                    result["issues"].append("Rosetta server not configured in Gemini CLI")
            else:
                result["integrations"]["gemini_cli"] = "Gemini CLI not configured"
                result["issues"].append("Gemini CLI settings not found")

        except Exception as e:
            result["integrations"]["gemini_cli"] = f"Exception: {e}"
            result["issues"].append(f"Gemini CLI check failed: {e}")

        result["status"] = "passed" if not result["issues"] else "warning"
        return result

    def diagnose_job_issues(self, job_id: Optional[str] = None) -> Dict:
        """Diagnose issues with job execution."""
        print(f"🔍 Diagnosing job issues{' for ' + job_id if job_id else ''}...")

        result = {"status": "unknown", "issues": [], "jobs": []}

        try:
            if not JOBS_DIR.exists():
                result["issues"].append("Jobs directory does not exist")
                result["status"] = "failed"
                return result

            # List all jobs
            job_dirs = [d for d in JOBS_DIR.iterdir() if d.is_dir()]

            if not job_dirs:
                result["status"] = "passed"
                result["issues"].append("No jobs found")
                print("ℹ️ No jobs to diagnose")
                return result

            for job_dir in job_dirs:
                if job_id and job_dir.name != job_id:
                    continue

                job_info = {
                    "job_id": job_dir.name,
                    "status": "unknown",
                    "issues": [],
                    "files": {}
                }

                # Check job files
                metadata_file = job_dir / "metadata.json"
                log_file = job_dir / "job.log"
                result_file = job_dir / "result.json"

                # Check metadata
                if metadata_file.exists():
                    try:
                        with open(metadata_file, 'r') as f:
                            metadata = json.load(f)
                        job_info["files"]["metadata"] = "OK"
                        job_info["status"] = metadata.get("status", "unknown")
                    except Exception as e:
                        job_info["files"]["metadata"] = f"Error: {e}"
                        job_info["issues"].append("Cannot read metadata file")
                else:
                    job_info["files"]["metadata"] = "Missing"
                    job_info["issues"].append("Metadata file missing")

                # Check log file
                if log_file.exists():
                    try:
                        log_size = log_file.stat().st_size
                        job_info["files"]["log"] = f"OK ({log_size} bytes)"

                        # Check for errors in log
                        with open(log_file, 'r') as f:
                            log_content = f.read()

                        if "ERROR" in log_content or "Exception" in log_content:
                            job_info["issues"].append("Errors found in log file")

                    except Exception as e:
                        job_info["files"]["log"] = f"Error: {e}"
                        job_info["issues"].append("Cannot read log file")
                else:
                    job_info["files"]["log"] = "Missing"
                    job_info["issues"].append("Log file missing")

                # Check result file
                if result_file.exists():
                    try:
                        with open(result_file, 'r') as f:
                            result_data = json.load(f)
                        job_info["files"]["result"] = "OK"
                    except Exception as e:
                        job_info["files"]["result"] = f"Error: {e}"
                        job_info["issues"].append("Cannot read result file")
                else:
                    if job_info["status"] in ["completed", "success"]:
                        job_info["files"]["result"] = "Missing (should exist)"
                        job_info["issues"].append("Result file missing for completed job")
                    else:
                        job_info["files"]["result"] = "Not created yet"

                result["jobs"].append(job_info)

            # Overall status
            if any(job["issues"] for job in result["jobs"]):
                result["status"] = "issues_found"
            else:
                result["status"] = "passed"
                print(f"✅ All {len(result['jobs'])} jobs look healthy")

        except Exception as e:
            result["status"] = "failed"
            result["issues"].append(f"Job diagnosis failed: {e}")

        return result

    def run_full_diagnosis(self) -> Dict:
        """Run complete diagnostic check."""
        print("🔧 Running Full MCP Server Diagnosis")
        print("=" * 50)

        diagnosis_result = {
            "timestamp": datetime.now().isoformat(),
            "overall_status": "unknown",
            "checks": {},
            "summary": {},
            "recommendations": []
        }

        # Define all checks
        checks = [
            ("Server Startup", self.check_server_startup),
            ("Python Environment", self.check_python_environment),
            ("Tool Availability", self.check_tool_availability),
            ("File Access", self.check_file_access),
            ("Job Manager", self.check_job_manager),
            ("MCP Integration", self.check_mcp_integration),
            ("Job Issues", lambda: self.diagnose_job_issues())
        ]

        passed = 0
        failed = 0
        warnings = 0

        for check_name, check_func in checks:
            print(f"\n📋 {check_name}")
            print("-" * 30)

            try:
                check_result = check_func()
                diagnosis_result["checks"][check_name] = check_result

                if check_result["status"] == "passed":
                    passed += 1
                elif check_result["status"] == "failed":
                    failed += 1
                    print(f"❌ {check_name} FAILED")
                    for issue in check_result.get("issues", []):
                        print(f"   - {issue}")
                elif check_result["status"] == "warning":
                    warnings += 1
                    print(f"⚠️ {check_name} has warnings")
                    for issue in check_result.get("issues", []):
                        print(f"   - {issue}")

            except Exception as e:
                failed += 1
                diagnosis_result["checks"][check_name] = {
                    "status": "error",
                    "issues": [str(e)]
                }
                print(f"💥 {check_name} ERROR: {e}")

        # Summary
        total_checks = len(checks)
        diagnosis_result["summary"] = {
            "total_checks": total_checks,
            "passed": passed,
            "failed": failed,
            "warnings": warnings
        }

        # Overall status
        if failed == 0 and warnings == 0:
            diagnosis_result["overall_status"] = "healthy"
        elif failed == 0:
            diagnosis_result["overall_status"] = "warning"
        else:
            diagnosis_result["overall_status"] = "issues"

        # Generate recommendations
        self._generate_recommendations(diagnosis_result)

        # Print summary
        print("\n" + "=" * 50)
        print("📊 DIAGNOSIS SUMMARY")
        print("=" * 50)
        print(f"Overall Status: {diagnosis_result['overall_status'].upper()}")
        print(f"Checks Passed: {passed}/{total_checks}")
        print(f"Checks Failed: {failed}/{total_checks}")
        print(f"Checks with Warnings: {warnings}/{total_checks}")

        if diagnosis_result["recommendations"]:
            print(f"\n📝 RECOMMENDATIONS:")
            for i, rec in enumerate(diagnosis_result["recommendations"], 1):
                print(f"{i}. {rec}")

        return diagnosis_result

    def _generate_recommendations(self, diagnosis_result: Dict):
        """Generate recommendations based on diagnosis results."""
        recommendations = []

        checks = diagnosis_result["checks"]

        # Server startup issues
        if checks.get("Server Startup", {}).get("status") == "failed":
            recommendations.append("Fix server startup issues before proceeding with other tests")

        # Environment issues
        if checks.get("Python Environment", {}).get("status") == "failed":
            recommendations.append("Install missing Python packages: pip install fastmcp loguru")

        # File access issues
        if checks.get("File Access", {}).get("status") == "failed":
            recommendations.append("Check file permissions and create missing directories")

        # Integration issues
        mcp_check = checks.get("MCP Integration", {})
        if mcp_check.get("status") in ["failed", "warning"]:
            claude_status = mcp_check.get("integrations", {}).get("claude_code")
            if "Not registered" in str(claude_status):
                recommendations.append("Register with Claude Code: claude mcp add rosetta -- $(pwd)/env/bin/python $(pwd)/src/server.py")
            elif "connection issue" in str(claude_status):
                recommendations.append("Restart Claude Code to reconnect to MCP server")

        # Job issues
        job_check = checks.get("Job Issues", {})
        if job_check.get("status") == "issues_found":
            recommendations.append("Check job logs for specific error messages: cat jobs/*/job.log")

        diagnosis_result["recommendations"] = recommendations

def main():
    """Run diagnosis from command line."""
    import argparse

    parser = argparse.ArgumentParser(description="Debug Rosetta MCP server")
    parser.add_argument("--job-id", help="Diagnose specific job")
    parser.add_argument("--quick", action="store_true", help="Run quick checks only")
    parser.add_argument("--save", help="Save results to file")

    args = parser.parse_args()

    debugger = MCPDebugger()

    if args.job_id:
        result = debugger.diagnose_job_issues(args.job_id)
    elif args.quick:
        # Quick checks only
        result = {}
        result["Server Startup"] = debugger.check_server_startup()
        result["Tool Availability"] = debugger.check_tool_availability()
    else:
        # Full diagnosis
        result = debugger.run_full_diagnosis()

    if args.save:
        save_path = Path(args.save)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        with open(save_path, 'w') as f:
            json.dump(result, f, indent=2)
        print(f"\n💾 Results saved to: {save_path}")

if __name__ == "__main__":
    main()