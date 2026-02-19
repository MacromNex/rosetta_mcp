# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Is

An MCP (Model Context Protocol) server exposing Rosetta molecular modeling tools via FastMCP. Provides 14 tools for protein structure refinement, docking, loop modeling, ligand docking, and mutation stability analysis (ddG). PyRosetta is optional — the server gracefully degrades to demo mode without it.

## Common Commands

```bash
# Setup environment
bash quick_setup.sh

# Run the MCP server
./env/bin/python src/server.py

# Run tests
./env/bin/python test_server.py

# Run a script directly (CLI mode)
./env/bin/python scripts/protein_refinement.py --input examples/data/test_input.pdb --trajectories 2
./env/bin/python scripts/ddg_calculations.py --demo

# Docker build (without PyRosetta)
docker build -t rosetta-mcp .

# Docker build (with PyRosetta, academic license)
docker build --build-arg PYROSETTA_CHANNEL=rosettacommons -t rosetta-mcp .

# Register with Claude Code
claude mcp add rosetta -- ./env/bin/python src/server.py
```

## Architecture

### Three-Layer Design

1. **MCP Server** (`src/server.py`) — FastMCP decorators exposing 14 tools. Adds `SCRIPTS_DIR` and `SCRIPT_DIR` to `sys.path` so scripts are importable.

2. **Job Manager** (`src/jobs/manager.py`) — Runs long tasks as subprocesses via `mamba run -p ./env python <script>`. Each job gets a UUID-based directory under `jobs/` with `metadata.json`, `job.log`, and `output.json`.

3. **Domain Scripts** (`scripts/`) — Five standalone Python scripts, each with identical structure: `DEFAULT_CONFIG` dict, `load_pyrosetta()` lazy loader, `run_*()` main function, and `main()` CLI entry point.

### Tool Categories

| Category | Tools | Execution |
|---|---|---|
| **Sync** (< 10 min) | `refine_protein_structure`, `calculate_ddg`, `validate_pdb_structure`, `list_example_structures` | Direct function import from `scripts/` |
| **Submit** (> 10 min) | `submit_protein_docking`, `submit_loop_modeling`, `submit_ligand_docking`, `submit_large_refinement`, `submit_batch_refinement` | Subprocess via `JobManager.submit_job()` |
| **Job Mgmt** | `get_job_status`, `get_job_result`, `get_job_log`, `cancel_job`, `list_jobs` | Read from `jobs/{id}/metadata.json` |

### Key Patterns

- **PyRosetta is always lazy-loaded** inside functions, never at module level. This allows the server and tests to start without PyRosetta installed.
- **Config merging**: `DEFAULT_CONFIG < file config < CLI args < kwargs`. Scripts merge with `{**DEFAULT_CONFIG, **(config or {}), **kwargs}`.
- **All tools return** `{"status": "success"|"error", ...}` dicts consistently.
- **Demo mode**: Every script supports `--demo` which creates synthetic PDB files and uses reduced parameters.
- **Job subprocess command**: `["mamba", "run", "-p", "./env", "python", script_path, "--input", ...]` with `cwd` set to the repo root.

### Adding a New Tool

1. Create `scripts/new_tool.py` following the existing template (DEFAULT_CONFIG, load_pyrosetta, run_*, main).
2. For sync: import `run_*` in `server.py` and add `@mcp.tool()` wrapper with try/except ImportError.
3. For async: add `@mcp.tool()` that calls `job_manager.submit_job(script_path, args)`.
4. Add config template in `configs/`.

## Environment

- Python 3.10 via conda/mamba in `./env`
- Dependencies: fastmcp, loguru, click, pandas, numpy, tqdm, biopython, matplotlib
- PyRosetta: optional, requires license from RosettaCommons (`conda install -c rosettacommons pyrosetta`)
- The `repo/` directory contains the full Rosetta C++ source (gitignored, not needed at runtime)
