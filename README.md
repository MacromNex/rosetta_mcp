# Rosetta MCP

> Molecular modeling and protein design tools powered by PyRosetta and integrated with Claude Code through the Model Context Protocol (MCP)

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

The Rosetta MCP provides comprehensive molecular modeling capabilities through PyRosetta, offering both standalone scripts and MCP integration for use with Claude Code and other AI assistants. This toolkit enables protein structure refinement, protein-protein docking, loop modeling, ligand docking, and stability analysis.

### Features
- **Protein Structure Refinement**: High-resolution structure optimization using Rosetta Relax protocol
- **Protein-Protein Docking**: Rigid-body docking for predicting protein complexes
- **Loop Modeling**: Reconstruct missing or flexible loop regions
- **Ligand Docking**: Protein-ligand complex optimization with interface scoring
- **Stability Analysis**: ΔΔG calculations for mutation effects
- **Batch Processing**: Submit multiple structures for concurrent processing
- **Job Management**: Asynchronous execution with status tracking and log monitoring

### Directory Structure
```
./
├── README.md               # This file
├── env/                    # Conda environment
├── src/
│   └── server.py           # MCP server (12 tools)
├── scripts/
│   ├── protein_refinement.py      # Structure refinement
│   ├── protein_docking.py         # Protein-protein docking
│   ├── loop_modeling.py           # Loop reconstruction
│   ├── ligand_docking.py          # Protein-ligand docking
│   ├── ddg_calculations.py        # ΔΔG stability analysis
│   └── lib/                       # Shared utilities
├── examples/
│   └── data/                      # Demo structures (8 files)
├── configs/                       # Configuration files (6 configs)
├── jobs/                          # Job execution directory
└── repo/                          # Original Rosetta repository
```

---

## Installation

### Prerequisites
- Conda or Mamba (mamba recommended for faster installation)
- Python 3.10+
- PyRosetta license (academic or commercial - optional for demo mode)

### Step 1: Create Environment

```bash
# Navigate to the MCP directory
cd /home/xux/Desktop/ProteinMCP/ProteinMCP/tool-mcps/rosetta_mcp

# Create conda environment (use mamba if available)
mamba create -p ./env python=3.10 -y
# or: conda create -p ./env python=3.10 -y

# Activate environment
mamba activate ./env
# or: conda activate ./env
```

### Step 2: Install Dependencies

```bash
# Install Python dependencies
pip install fastmcp loguru click pandas numpy tqdm biopython matplotlib

# Install MCP dependencies (if not already installed)
pip install fastmcp loguru
```

### Step 3: Install PyRosetta (Optional)

PyRosetta requires a license from RosettaCommons:

```bash
# Academic license
conda install -c rosettacommons pyrosetta

# Commercial license
conda install -c rosettacommons/label/commercial pyrosetta
```

**Note**: All scripts work in demo mode without PyRosetta and provide clear installation instructions.

### Step 4: Verify Installation

```bash
# Test imports without PyRosetta
python -c "import sys; sys.path.insert(0, 'src'); from server import mcp; print(f'Found {len(mcp.get_name_to_tool_map())} tools')"

# Test script functionality
python scripts/protein_refinement.py --help
```

---

## Local Usage (Scripts)

You can use the scripts directly without MCP for local processing.

### Available Scripts

| Script | Description | Example |
|--------|-------------|---------|
| `scripts/protein_refinement.py` | Structure refinement using Relax protocol | See below |
| `scripts/protein_docking.py` | Protein-protein docking | See below |
| `scripts/loop_modeling.py` | Loop reconstruction with CCD closure | See below |
| `scripts/ligand_docking.py` | Protein-ligand docking optimization | See below |
| `scripts/ddg_calculations.py` | ΔΔG stability calculations | See below |

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
- `--input, -i`: Input PDB file (required)
- `--output, -o`: Results output file (default: auto-generated)
- `--trajectories, -n`: Number of refinement trajectories (default: 5)
- `--cycles, -c`: Monte Carlo cycles per trajectory (default: 100)
- `--temperature, -t`: MC temperature (default: 2.0)

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
- `--input, -i`: Protein complex PDB file (required)
- `--trajectories, -n`: Number of docking trajectories (default: 10)
- `--chain_break, -c`: Residue where chain A ends (auto-detect if not provided)

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
- `--input, -i`: Protein-ligand complex PDB file (required)
- `--ligand_chain, -l`: Ligand chain ID (auto-detect if not provided)
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
- `--repack_radius, -r`: Repacking radius in Angstroms (default: 8.0)

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

### Quick Operations (Synchronous API)

These tools return results immediately (< 10 minutes):

| Tool | Description | Parameters |
|------|-------------|------------|
| `refine_protein_structure` | Fast structure refinement | `input_file`, `trajectories`, `cycles`, `temperature`, `output_file` |
| `calculate_ddg` | ΔΔG stability analysis | `input_file`, `mutations`, `trajectories`, `repack_radius`, `minimize`, `output_file` |

### Long-Running Tasks (Submit API)

These tools return a job_id for tracking (> 10 minutes):

| Tool | Description | Parameters |
|------|-------------|------------|
| `submit_protein_docking` | Protein-protein docking | `input_file`, `trajectories`, `chain_break`, `use_centroid_stage`, `use_fullatom_stage`, `output_dir`, `job_name` |
| `submit_loop_modeling` | Loop reconstruction | `input_file`, `loop_start`, `loop_end`, `trajectories`, `use_centroid_stage`, `use_fullatom_stage`, `output_dir`, `job_name` |
| `submit_ligand_docking` | Protein-ligand docking | `input_file`, `ligand_chain`, `trajectories`, `perturbation_cycles`, `repack_sidechains`, `minimize_final`, `output_dir`, `job_name` |
| `submit_large_refinement` | Large-scale refinement | `input_file`, `trajectories`, `cycles`, `temperature`, `output_dir`, `job_name` |
| `submit_batch_refinement` | Batch processing | `input_files`, `trajectories`, `cycles`, `temperature`, `output_dir`, `job_name` |

### Job Management Tools

| Tool | Description |
|------|-------------|
| `get_job_status` | Check job progress and status |
| `get_job_result` | Get completed job results |
| `get_job_log` | View execution logs with tail option |
| `cancel_job` | Cancel running job |
| `list_jobs` | List all jobs with optional status filter |

### Utility Tools

| Tool | Description |
|------|-------------|
| `validate_pdb_structure` | Validate PDB format and compatibility |
| `list_example_structures` | List available demo structures |

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

| File | Size | Description | Use With |
|------|------|-------------|----------|
| `test_in.pdb` | 157KB | General test protein structure | All tools |
| `design_in.pdb` | 152KB | Protein design test structure | Refinement, ΔΔG |
| `test_complex.pdb` | 802B | Two-chain complex for docking | Protein docking |
| `test_docking.pdb` | 644B | Minimal docking test case | Protein docking |
| `test_ligand.pdb` | 565B | Protein-ligand complex | Ligand docking |
| `test_loop.pdb` | 956B | Structure with loop region | Loop modeling |
| `test_input.pdb` | 403B | Minimal test structure | All tools |
| `ZN1.params` | 528B | Zinc parameter file | Ligand docking |

---

## Configuration Files

The `configs/` directory contains configuration templates:

| Config | Description | Key Parameters |
|--------|-------------|----------------|
| `default_config.json` | Global template and defaults | `trajectories`, `cycles`, `temperature`, `repack_radius` |
| `protein_refinement_config.json` | Refinement-specific settings | `trajectories`, `cycles`, `temperature` |
| `protein_docking_config.json` | Docking protocol parameters | `trajectories`, `chain_break`, `centroid_stage` |
| `loop_modeling_config.json` | Loop modeling settings | `loop_start`, `loop_end`, `trajectories` |
| `ligand_docking_config.json` | Ligand optimization parameters | `trajectories`, `perturbation_cycles`, `minimize_final` |
| `ddg_calculations_config.json` | ΔΔG calculation settings | `mutations`, `repack_radius`, `minimize` |

### Example Config Usage

```bash
# Use custom config file
python scripts/protein_refinement.py \
  --input examples/data/test_in.pdb \
  --config configs/protein_refinement_config.json

# Override config parameters via CLI
python scripts/protein_refinement.py \
  --input examples/data/test_in.pdb \
  --config configs/protein_refinement_config.json \
  --trajectories 10  # Overrides config value
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
python -c "import sys; sys.path.insert(0, 'src'); from server import mcp; print('Server loaded successfully')"
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
python src/server.py --help

# Check tool availability
python -c "
import sys; sys.path.insert(0, 'src')
from server import mcp
tools = mcp.get_name_to_tool_map()
print(f'Available tools: {list(tools.keys())}')
"
```

### PyRosetta Issues

**Problem:** PyRosetta not available
- All scripts work in `--demo` mode without PyRosetta
- Use demo mode to test functionality: `python scripts/protein_refinement.py --demo`
- Install PyRosetta with academic or commercial license

**Problem:** PyRosetta license errors
- Ensure valid license from RosettaCommons
- Check license file location and permissions
- Contact RosettaCommons support for license issues

### Job Issues

**Problem:** Job stuck in pending
```bash
# Check job directory
ls -la jobs/

# View job logs
python -c "
import sys; sys.path.insert(0, 'src')
from jobs.manager import job_manager
jobs = job_manager.list_jobs()
print(jobs)
"
```

**Problem:** Job failed
```
Use get_job_log with job_id "<job_id>" and tail 100 to see detailed error information
```

**Problem:** Out of disk space
```bash
# Clean old job directories
rm -rf jobs/*/  # Be careful - removes all job results
```

---

## Development

### Running Tests

```bash
# Activate environment
mamba activate ./env

# Test server functionality
python test_server.py

# Test individual scripts
python scripts/protein_refinement.py --demo --trajectories 1
python scripts/ddg_calculations.py --demo --mutations "10G"
```

### Starting Dev Server

```bash
# Run MCP server in dev mode
fastmcp dev src/server.py

# Test with Claude Code dev setup
export PYTHONPATH="$PWD/src:$PWD/scripts"
python src/server.py
```

### Adding New Tools

1. Create script in `scripts/` following existing patterns
2. Add tool wrapper in `src/server.py`
3. Create config file in `configs/`
4. Add tests and documentation
5. Update this README

### Performance Monitoring

```bash
# Monitor job execution
watch -n 5 'ls -la jobs/'

# Check system resources
top -p $(pgrep -f "python.*server.py")

# Monitor log files
tail -f jobs/*/job.log
```

---

## License

This project is based on the Rosetta molecular modeling suite and requires appropriate licensing for PyRosetta. The MCP integration code is provided under the MIT License.

## Credits

- Based on [Rosetta](https://www.rosettacommons.org/) molecular modeling suite
- PyRosetta Python interface
- Uses [FastMCP](https://github.com/jlowin/fastmcp) for Model Context Protocol integration
- Demo structures from Rosetta benchmark datasets