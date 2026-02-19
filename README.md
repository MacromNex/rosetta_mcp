# Rosetta MCP

> Model Control Protocol (MCP) server for Rosetta molecular modeling and protein design suite

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Local Usage (Scripts)](#local-usage-scripts)
- [MCP Server Installation](#mcp-server-installation)
- [Using with Claude Code](#using-with-claude-code)
- [Using with Gemini CLI](#using-with-gemini-cli)
- [Available Tools](#available-tools)
- [Examples](#examples)
- [Demo Data](#demo-data)
- [Configuration Files](#configuration-files)
- [Troubleshooting](#troubleshooting)
- [Development](#development)

## Overview

The Rosetta MCP server provides easy access to Rosetta's powerful molecular modeling and protein design tools through a modern Model Control Protocol interface. It enables protein structure refinement, protein-protein docking, loop modeling, ligand docking, and stability analysis (ΔΔG calculations) via both synchronous and asynchronous APIs.

### Features
- **Dual API Design**: Sync tools for quick analysis (<10 min), Submit API for long-running tasks (>10 min)
- **Comprehensive Protein Modeling**: Structure refinement, docking, loop modeling, ligand binding, stability analysis
- **Robust Error Handling**: Graceful PyRosetta dependency management with clear error messages
- **Job Management**: Full lifecycle tracking for long-running computational jobs
- **Batch Processing**: Support for processing multiple structures simultaneously
- **Demo Mode**: Built-in test data generation for immediate evaluation
- **Claude Code & Gemini CLI**: Ready integration with modern AI development environments

### Directory Structure
```
./
├── README.md               # This file
├── env/                    # Conda environment
├── src/
│   └── server.py           # MCP server (14 tools)
├── scripts/
│   ├── protein_refinement.py      # Structure refinement
│   ├── protein_docking.py         # Protein-protein docking
│   ├── loop_modeling.py           # Loop reconstruction
│   ├── ligand_docking.py          # Protein-ligand docking
│   ├── ddg_calculations.py        # ΔΔG stability analysis
│   └── lib/                       # Shared utilities
├── examples/
│   └── data/               # Demo data (7 PDB files, ~313 KB)
├── configs/                # Configuration files
├── jobs/                   # Job execution directory
└── reports/                # Documentation and test results
```

---

## Installation

### Quick Setup (Recommended)

Run the automated setup script:

```bash
cd rosetta_mcp
bash quick_setup.sh
```

The script will create the conda environment, install all dependencies, and display the Claude Code configuration. See `quick_setup.sh --help` for options like `--skip-env`.

**Note:** PyRosetta requires a separate license. All scripts work in demo mode without PyRosetta.

### Prerequisites
- Conda or Mamba (mamba recommended for faster installation)
- Python 3.10+
- PyRosetta license (academic or commercial - optional for demo mode)

### Manual Installation (Alternative)

If you prefer manual installation or need to customize the setup, follow `reports/step3_environment.md`:

```bash
# Navigate to the MCP directory
cd /home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp

# Create conda environment (use mamba if available)
mamba create -p ./env python=3.10 -y
# or: conda create -p ./env python=3.10 -y

# Activate environment
mamba activate ./env
# or: conda activate ./env

# Install Dependencies
mamba run -p ./env pip install numpy pandas loguru click tqdm
mamba run -p ./env pip install --force-reinstall --no-cache-dir fastmcp
mamba run -p ./env pip install biopython matplotlib

# Install MCP dependencies (if not already installed)
pip install fastmcp loguru --ignore-installed
```

### PyRosetta Installation (Optional)

PyRosetta requires a license from RosettaCommons:

```bash
# Academic license
pip install pyrosetta --find-links https://west.rosettacommons.org/pyrosetta/quarterly/release
```

**Note:** All scripts work in demo mode without PyRosetta and provide clear installation instructions when needed.

---

## Local Usage (Scripts)

You can use the scripts directly without MCP for local processing.

### Available Scripts

| Script | Description | Example |
|--------|-------------|---------|
| `scripts/protein_refinement.py` | Structure refinement using Relax protocol | See below |
| `scripts/protein_docking.py` | Protein-protein docking | See below |
| `scripts/loop_modeling.py` | Loop reconstruction | See below |
| `scripts/ligand_docking.py` | Protein-ligand docking | See below |
| `scripts/ddg_calculations.py` | ΔΔG stability analysis | See below |

### Script Examples

#### Protein Structure Refinement

```bash
# Activate environment
mamba activate ./env

# Run refinement with custom parameters
python scripts/protein_refinement.py \
  --input examples/data/test_in.pdb \
  --output results/refined_structure.json \
  --trajectories 5 \
  --cycles 100

# Run demo mode (no PyRosetta required)
python scripts/protein_refinement.py --demo --trajectories 2 --cycles 10
```

**Parameters:**
- `--input, -i`: Input PDB file path (required)
- `--output, -o`: Output JSON file path (default: auto-generated)
- `--trajectories, -n`: Number of refinement trajectories (default: 5)
- `--cycles, -c`: Monte Carlo cycles per trajectory (default: 100)
- `--temperature, -t`: Monte Carlo temperature (default: 2.0)
- `--config`: Configuration file (optional)

#### Protein-Protein Docking

```bash
python scripts/protein_docking.py \
  --input examples/data/test_complex.pdb \
  --trajectories 10 \
  --chain_break 150

# Demo mode with auto-generated complex
python scripts/protein_docking.py --demo --trajectories 5
```

**Parameters:**
- `--input, -i`: PDB file with both protein chains (required)
- `--trajectories, -n`: Number of docking trajectories (default: 10)
- `--chain_break, -c`: Residue number where chain A ends (auto-detect)
- `--output, -o`: Output file path (optional)

#### Loop Modeling

```bash
python scripts/loop_modeling.py \
  --input examples/data/test_loop.pdb \
  --loop_start 45 \
  --loop_end 52 \
  --trajectories 10

# Demo mode
python scripts/loop_modeling.py --demo --loop_start 10 --loop_end 15
```

**Parameters:**
- `--input, -i`: Input PDB file (required)
- `--loop_start, -s`: First residue of loop region (required)
- `--loop_end, -e`: Last residue of loop region (required)
- `--trajectories, -n`: Number of modeling trajectories (default: 10)

#### Ligand Docking

```bash
python scripts/ligand_docking.py \
  --input examples/data/test_ligand.pdb \
  --ligand_chain L \
  --trajectories 10

# Demo mode
python scripts/ligand_docking.py --demo --trajectories 5
```

**Parameters:**
- `--input, -i`: PDB file with protein-ligand complex (required)
- `--ligand_chain, -l`: Ligand chain ID (auto-detect if not specified)
- `--trajectories, -n`: Number of docking trajectories (default: 10)

#### ΔΔG Calculations

```bash
python scripts/ddg_calculations.py \
  --input examples/data/test_in.pdb \
  --mutations "A10G,B15L" \
  --trajectories 5

# Demo mode
python scripts/ddg_calculations.py --demo --mutations "10G,15A"
```

**Parameters:**
- `--input, -i`: Input PDB file (required)
- `--mutations, -m`: Mutations in format "A10G,B15L" or "10G,15L" (required)
- `--trajectories, -n`: Number of trajectories (default: 5)
- `--repack_radius, -r`: Sidechain repacking radius in Å (default: 8.0)

---

## MCP Server Installation

### Option 1: Using fastmcp (Recommended)

```bash
# Install MCP server for Claude Code
fastmcp install src/server.py --name rosetta
```

### Option 2: Manual Installation for Claude Code

```bash
# Add MCP server to Claude Code
claude mcp add rosetta -- $(pwd)/env/bin/python $(pwd)/src/server.py

# Verify installation
claude mcp list
```

### Option 3: Configure in settings.json

Add to `~/.claude/settings.json`:

```json
{
  "mcpServers": {
    "rosetta": {
      "command": "/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/env/bin/python",
      "args": ["/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/src/server.py"]
    }
  }
}
```

---

## Using with Claude Code

After installing the MCP server, you can use it directly in Claude Code.

### Quick Start

```bash
# Start Claude Code
claude
```

### Example Prompts

#### Tool Discovery
```
What tools are available from rosetta?
```

#### Basic Structure Refinement
```
Use refine_protein_structure with input file @examples/data/test_in.pdb
```

#### Mutation Analysis
```
Use calculate_ddg with input file @examples/data/test_in.pdb and mutations "A10G,B15L"
```

#### Long-Running Tasks (Submit API)
```
Submit protein docking for @examples/data/test_complex.pdb with trajectories 20
Then check the job status
```

#### Batch Processing
```
Submit batch refinement for these files:
- @examples/data/test_in.pdb
- @examples/data/design_in.pdb
```

#### Job Management
```
List all running jobs
Get job status for job "abc123"
Get the latest logs for job "abc123"
```

### Using @ References

In Claude Code, use `@` to reference files and directories:

| Reference | Description |
|-----------|-------------|
| `@examples/data/test_in.pdb` | Reference a specific structure file |
| `@configs/default_config.json` | Reference a config file |
| `@results/` | Reference output directory |

---

## Using with Gemini CLI

### Configuration

Add to `~/.gemini/settings.json`:

```json
{
  "mcpServers": {
    "rosetta": {
      "command": "/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/env/bin/python",
      "args": ["/home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/src/server.py"]
    }
  }
}
```

### Example Prompts

```bash
# Start Gemini CLI
gemini

# Example prompts (same as Claude Code)
> What tools are available?
> Use refine_protein_structure with file examples/data/test_in.pdb
> Submit loop modeling for examples/data/test_loop.pdb with loop_start 45 and loop_end 52
```

---

## Available Tools

### Quick Operations (Sync API)

These tools return results immediately (< 10 minutes):

| Tool | Description | Parameters |
|------|-------------|------------|
| `refine_protein_structure` | Fast protein refinement | `input_file`, `trajectories`, `cycles`, `temperature` |
| `calculate_ddg` | Mutation stability analysis | `input_file`, `mutations`, `trajectories`, `repack_radius` |
| `validate_pdb_structure` | PDB file validation | `input_file` |
| `list_example_structures` | List available demo files | None |

### Long-Running Tasks (Submit API)

These tools return a job_id for tracking (> 10 minutes):

| Tool | Description | Parameters |
|------|-------------|------------|
| `submit_protein_docking` | Protein-protein docking | `input_file`, `trajectories`, `chain_break` |
| `submit_loop_modeling` | Loop reconstruction | `input_file`, `loop_start`, `loop_end`, `trajectories` |
| `submit_ligand_docking` | Protein-ligand docking | `input_file`, `ligand_chain`, `trajectories` |
| `submit_large_refinement` | Large-scale refinement | `input_file`, `trajectories`, `cycles` |
| `submit_batch_refinement` | Batch structure refinement | `input_files`, `trajectories`, `cycles` |

### Job Management Tools

| Tool | Description |
|------|-------------|
| `get_job_status` | Check job progress |
| `get_job_result` | Get results when completed |
| `get_job_log` | View execution logs |
| `cancel_job` | Cancel running job |
| `list_jobs` | List all jobs |

---

## Examples

### Example 1: Quick Structure Analysis

**Goal:** Quickly refine and analyze a protein structure

**Using Script:**
```bash
python scripts/protein_refinement.py \
  --input examples/data/test_in.pdb \
  --trajectories 3 \
  --cycles 50
```

**Using MCP (in Claude Code):**
```
Validate the structure at @examples/data/test_in.pdb then refine it with 3 trajectories and 50 cycles
```

**Expected Output:**
- Validation results with structure statistics
- Refined structure with energy scores
- Runtime: ~2-3 minutes

### Example 2: Mutation Stability Analysis

**Goal:** Analyze how mutations affect protein stability

**Using Script:**
```bash
python scripts/ddg_calculations.py \
  --input examples/data/test_in.pdb \
  --mutations "A10G,A10L,A10F" \
  --trajectories 5
```

**Using MCP (in Claude Code):**
```
Calculate ΔΔG for mutations "A10G,A10L,A10F" in @examples/data/test_in.pdb with 5 trajectories
```

**Expected Output:**
- ΔΔG values for each mutation
- Stabilizing vs destabilizing classification
- Per-mutation success/failure status

### Example 3: Protein-Protein Docking

**Goal:** Predict binding conformation of two proteins

**Using Script:**
```bash
python scripts/protein_docking.py \
  --input examples/data/test_complex.pdb \
  --trajectories 20 \
  --chain_break 150
```

**Using MCP (in Claude Code):**
```
Submit protein docking for @examples/data/test_complex.pdb with chain break at residue 150 and 20 trajectories
Check the job status every few minutes
When complete, get the results
```

**Expected Output:**
- Job ID for tracking
- Docking scores and interface analysis
- Runtime: ~15-30 minutes

### Example 4: Loop Reconstruction

**Goal:** Model a missing loop region

**Using Script:**
```bash
python scripts/loop_modeling.py \
  --input examples/data/test_loop.pdb \
  --loop_start 45 \
  --loop_end 52 \
  --trajectories 15
```

**Using MCP (in Claude Code):**
```
Submit loop modeling for @examples/data/test_loop.pdb for residues 45 to 52 with 15 trajectories
Monitor the job progress with logs
```

**Expected Output:**
- Loop models with RMSD scores
- Best conformations ranked by energy
- Runtime: ~10-20 minutes

### Example 5: Batch Processing

**Goal:** Process multiple structures at once

**Using Script:**
```bash
for f in examples/data/*.pdb; do
  python scripts/protein_refinement.py --input "$f" --trajectories 5
done
```

**Using MCP (in Claude Code):**
```
Submit batch refinement for all PDB files in @examples/data/ with 5 trajectories each
Track the progress and get results when complete
```

---

## Demo Data

The `examples/data/` directory contains sample data for testing:

| File | Description | Size | Use With |
|------|-------------|------|----------|
| `test_input.pdb` | Basic test structure | 0.4 KB | All tools |
| `test_loop.pdb` | Loop modeling test | 0.9 KB | Loop modeling |
| `test_ligand.pdb` | Protein-ligand complex | 0.6 KB | Ligand docking |
| `test_docking.pdb` | Protein docking test | 0.6 KB | Protein docking |
| `test_complex.pdb` | Multi-chain complex | 0.8 KB | Docking, ddG |
| `test_in.pdb` | Large test structure | 153.4 KB | Refinement, all tools |
| `design_in.pdb` | Design test structure | 148.3 KB | Refinement, design |
| `ZN1.params` | Zinc parameter file | 0.5 KB | Ligand docking |

---

## Configuration Files

The `configs/` directory contains configuration templates:

| Config | Description | Key Parameters |
|--------|-------------|---------------|
| `default_config.json` | Global defaults template | `trajectories`, `cycles`, `temperature` |
| `protein_refinement_config.json` | Refinement parameters | `trajectories`, `cycles`, `temperature` |
| `protein_docking_config.json` | Docking parameters | `trajectories`, `chain_break`, `use_centroid_stage` |
| `loop_modeling_config.json` | Loop modeling parameters | `loop_start`, `loop_end`, `trajectories` |
| `ligand_docking_config.json` | Ligand docking parameters | `trajectories`, `ligand_chain`, `perturbation_cycles` |
| `ddg_calculations_config.json` | ΔΔG parameters | `mutations`, `repack_radius`, `minimize` |

### Config Example

```json
{
  "trajectories": 5,
  "cycles": 100,
  "temperature": 2.0,
  "output_prefix": "rosetta_output"
}
```

---

## Troubleshooting

### Environment Issues

**Problem:** Environment not found
```bash
# Recreate environment
mamba create -p ./env python=3.10 -y
mamba activate ./env
pip install fastmcp loguru click pandas numpy tqdm biopython matplotlib
```

**Problem:** Import errors
```bash
# Verify installation
python -c "import fastmcp, loguru, numpy, pandas; print('All imports successful')"
```

### MCP Issues

**Problem:** Server not found in Claude Code
```bash
# Check MCP registration
claude mcp list

# Re-add if needed
claude mcp remove rosetta
claude mcp add rosetta -- $(pwd)/env/bin/python $(pwd)/src/server.py
```

**Problem:** Tools not working
```bash
# Test server directly
python -c "
import sys
sys.path.append('src')
from server import mcp
print('Available tools:', list(mcp.list_tools().keys()))
"
```

**Problem:** PyRosetta not available
```bash
# Check PyRosetta installation
python -c "import pyrosetta; print('PyRosetta available')"

# If not available, scripts work in demo mode with helpful error messages
```

**Problem:** Connection issues
```bash
# Test server startup
python src/server.py
# Should start without errors and show "MCP Server ready"
```

### Job Issues

**Problem:** Job stuck in pending
```bash
# Check job directory
ls -la jobs/

# View job metadata
cat jobs/<job_id>/metadata.json
```

**Problem:** Job failed
```
Use get_job_log with job_id "<job_id>" and tail 100 to see error details
```

**Problem:** Long queue times
```
Use list_jobs to see all pending jobs and cancel old ones if needed
```

### File Access Issues

**Problem:** File not found errors
```bash
# Check file paths are absolute
ls -la examples/data/
# Use full paths: /home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp/examples/data/
```

**Problem:** Permission errors
```bash
# Fix permissions
chmod +x scripts/*.py
chmod -R 755 examples/data/
```

---

## Development

### Running Tests

```bash
# Activate environment
mamba activate ./env

# Run test suite
python test_server.py

# Test individual scripts
python scripts/protein_refinement.py --demo
```

### Starting Dev Server

```bash
# Run MCP server in development mode
fastmcp dev src/server.py

# The server will be available at localhost:6274
```

### Adding New Tools

To add new tools to the MCP server:

1. Create a new script in `scripts/`
2. Add corresponding config in `configs/`
3. Add MCP tool wrapper in `src/server.py`
4. Test with demo mode first
5. Update this README with new tool documentation

---

## License

This project uses the Rosetta software suite. Please comply with Rosetta licensing terms:
- Academic use: Free with registration
- Commercial use: Requires license from RosettaCommons

## Credits

Based on [Rosetta Commons](https://www.rosettacommons.org/) software suite.
Original PyRosetta demos and protocols adapted for MCP integration.

---

**Note:** This MCP server enables AI assistants to perform sophisticated molecular modeling tasks. All tools handle PyRosetta licensing gracefully and provide demo modes for testing without the full Rosetta installation.