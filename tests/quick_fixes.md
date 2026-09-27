# Quick Fixes for Common Rosetta MCP Issues

This document provides quick solutions for common problems encountered with the Rosetta MCP server.

## 🚨 Emergency Debugging Commands

### Quick Health Check
```bash
# Run basic server health check
cd /home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp
./env/bin/python tests/debug_utils.py --quick
```

### Full Diagnosis
```bash
# Run complete diagnostic
./env/bin/python tests/debug_utils.py --save reports/diagnosis.json
```

### Check Specific Job
```bash
# Diagnose issues with specific job
./env/bin/python tests/debug_utils.py --job-id JOB_ID_HERE
```

## 🔧 Common Issues and Fixes

### Issue 1: Server Won't Import
```bash
# Symptoms: ImportError when starting server
# Fix 1: Check Python environment
./env/bin/python -c "import sys; print(sys.path)"

# Fix 2: Verify FastMCP installation
./env/bin/pip list | grep fastmcp

# Fix 3: Reinstall if needed
./env/bin/pip install --upgrade fastmcp loguru
```

### Issue 2: Claude Code Can't Connect
```bash
# Symptoms: "Failed to connect" in claude mcp list
# Fix 1: Check registration
claude mcp list | grep rosetta

# Fix 2: Re-register server
claude mcp remove rosetta
claude mcp add rosetta -- $(pwd)/env/bin/python $(pwd)/src/server.py

# Fix 3: Verify server starts manually
./env/bin/python -c "from src.server import mcp; print('Server OK')"
```

### Issue 3: Tools Not Found
```bash
# Symptoms: Tool not available errors
# Fix 1: List available tools
./env/bin/python -c "
import asyncio
from src.server import mcp
async def main():
    tools = await mcp.get_tools()
    for tool in tools: print(tool)
asyncio.run(main())
"

# Fix 2: Check server file
grep "@mcp.tool" src/server.py
```

### Issue 4: File Access Problems
```bash
# Symptoms: FileNotFoundError for example files
# Fix 1: Check example data exists
ls -la examples/data/*.pdb

# Fix 2: Create test file if missing
echo "ATOM      1  N   ALA A   1      20.154  16.967  10.000  1.00 20.00           N" > examples/data/test.pdb

# Fix 3: Check permissions
chmod 644 examples/data/*.pdb
```

### Issue 5: Job Manager Issues
```bash
# Symptoms: Job submission fails
# Fix 1: Check jobs directory
ls -la jobs/

# Fix 2: Create jobs directory if missing
mkdir -p jobs

# Fix 3: Check job manager
./env/bin/python -c "from src.jobs.manager import job_manager; print(job_manager.list_jobs())"
```

### Issue 6: Path Resolution Problems
```bash
# Symptoms: Scripts not found
# Fix 1: Check scripts directory
ls -la scripts/

# Fix 2: Verify PYTHONPATH
export PYTHONPATH=$(pwd):$PYTHONPATH

# Fix 3: Use absolute paths
./env/bin/python -c "
import sys
from pathlib import Path
sys.path.insert(0, str(Path.cwd()))
print('Path fixed')
"
```

## 🔍 Advanced Debugging

### Debug Server Startup Step by Step
```bash
# Step 1: Test Python environment
./env/bin/python --version

# Step 2: Test imports
./env/bin/python -c "import fastmcp; print('FastMCP OK')"
./env/bin/python -c "import loguru; print('Loguru OK')"

# Step 3: Test server creation
./env/bin/python -c "from fastmcp import FastMCP; mcp = FastMCP('test'); print('MCP creation OK')"

# Step 4: Test server import
./env/bin/python -c "from src.server import mcp; print('Server import OK')"

# Step 5: Test tool listing
./env/bin/python -c "
import asyncio
from src.server import mcp
async def main():
    tools = await mcp.get_tools()
    print(f'Found {len(tools)} tools')
asyncio.run(main())
"
```

### Debug Tool Execution
```bash
# Test individual tool functions (not through MCP)
./env/bin/python -c "
import sys
sys.path.append('src')
from pathlib import Path

# Test validation function exists
try:
    from server import validate_pdb_structure
    print('Validation tool function found')
except Exception as e:
    print(f'Tool function error: {e}')
"
```

### Debug Job Execution
```bash
# Check job execution environment
./env/bin/python -c "
from pathlib import Path
scripts_dir = Path('scripts')
if scripts_dir.exists():
    scripts = list(scripts_dir.glob('*.py'))
    print(f'Found {len(scripts)} script files')
    for script in scripts:
        print(f'  - {script.name}')
else:
    print('Scripts directory not found')
"
```

### Debug MCP Protocol
```bash
# Test FastMCP dev server
timeout 5s ./env/bin/fastmcp dev src/server.py || echo "Dev server test complete"

# Test with MCP inspector (if available)
./env/bin/fastmcp dev src/server.py &
sleep 2
curl -s http://localhost:6274/ && echo "MCP inspector accessible"
```

## 🚨 Emergency Recovery

### Reset Everything
```bash
# 1. Remove from Claude Code
claude mcp remove rosetta

# 2. Clean job directory
rm -rf jobs/*

# 3. Restart with fresh registration
claude mcp add rosetta -- $(pwd)/env/bin/python $(pwd)/src/server.py

# 4. Verify connection
claude mcp list | grep rosetta
```

### Reinstall Dependencies
```bash
# Clean reinstall of MCP dependencies
./env/bin/pip uninstall -y fastmcp
./env/bin/pip install fastmcp loguru

# Verify installation
./env/bin/pip show fastmcp
```

### Check System Resources
```bash
# Check disk space
df -h .

# Check memory
free -h

# Check running processes
ps aux | grep python | grep server.py
```

## 📝 Logging and Monitoring

### Enable Detailed Logging
```bash
# Create logs directory
mkdir -p logs

# Run server with debug logging
export LOGURU_LEVEL=DEBUG
./env/bin/python src/server.py
```

### Monitor Job Execution
```bash
# Watch job logs in real-time
tail -f jobs/*/job.log

# Monitor job directory
watch -n 2 'ls -la jobs/'

# Check job status periodically
watch -n 5 './env/bin/python -c "from src.jobs.manager import job_manager; print(job_manager.list_jobs())"'
```

### Performance Monitoring
```bash
# Monitor server resource usage
top -p $(pgrep -f "src/server.py")

# Check file handles
lsof | grep server.py

# Monitor network connections
netstat -tlnp | grep python
```

## 📞 When All Else Fails

### Collect Debug Information
```bash
# Create comprehensive debug report
./env/bin/python tests/debug_utils.py --save reports/emergency_debug.json

# Collect system information
cat /etc/os-release > reports/system_info.txt
./env/bin/python --version >> reports/system_info.txt
./env/bin/pip freeze >> reports/system_info.txt

# Collect log files
cp -r logs/ reports/logs_backup/ 2>/dev/null || true
cp -r jobs/ reports/jobs_backup/ 2>/dev/null || true
```

### Manual Server Test
```bash
# Test server completely manually
./env/bin/python -c "
print('=== Manual Server Test ===')
try:
    from src.server import mcp
    print('✓ Server imports OK')

    import asyncio
    async def test():
        tools = await mcp.get_tools()
        print(f'✓ Found {len(tools)} tools')

        tool = await mcp.get_tool('validate_pdb_structure')
        print('✓ Can access validation tool')

        return True

    success = asyncio.run(test())
    if success:
        print('✓ All manual tests passed')
    else:
        print('✗ Manual tests failed')

except Exception as e:
    print(f'✗ Manual test failed: {e}')
    import traceback
    traceback.print_exc()
"
```

Remember: Most issues can be resolved by following the diagnostic output and applying the suggested fixes systematically.