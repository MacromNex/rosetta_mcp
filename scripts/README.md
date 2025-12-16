# MCP Rosetta Scripts

Clean, self-contained scripts extracted from use cases for MCP tool wrapping.

## Design Principles

1. **Minimal Dependencies**: Only essential packages imported (standard library only)
2. **Self-Contained**: PyRosetta functions are lazy loaded, utility functions inlined where possible
3. **Configurable**: Parameters in config files, not hardcoded in scripts
4. **MCP-Ready**: Each script has a main function ready for MCP wrapping
5. **Graceful Error Handling**: Clear PyRosetta dependency messages

## Scripts Overview

| Script | Description | Main Function | Config |
|--------|-------------|---------------|--------|
| `protein_refinement.py` | Structure refinement (Relax protocol) | `run_protein_refinement()` | `configs/protein_refinement_config.json` |
| `protein_docking.py` | Protein-protein docking | `run_protein_docking()` | `configs/protein_docking_config.json` |
| `loop_modeling.py` | Loop reconstruction | `run_loop_modeling()` | `configs/loop_modeling_config.json` |
| `ligand_docking.py` | Protein-ligand docking | `run_ligand_docking()` | `configs/ligand_docking_config.json` |
| `ddg_calculations.py` | ΔΔG stability analysis | `run_ddg_calculations()` | `configs/ddg_calculations_config.json` |

## Dependencies

### Essential Dependencies (Standard Library)
- `argparse` - CLI argument parsing
- `logging` - Structured logging
- `os`, `sys` - System operations
- `pathlib` - Path handling
- `typing` - Type hints
- `json` - Configuration loading

### Repository Dependencies
- **PyRosetta** (lazy loaded) - All scripts require PyRosetta but handle its absence gracefully

## Usage

### Basic Usage

```bash
# Activate environment (prefer mamba over conda)
mamba activate ./env  # or: conda activate ./env

# Run a script with demo data
python scripts/protein_refinement.py --demo

# Run with real data
python scripts/protein_refinement.py --input examples/data/protein.pdb --output results/refined.json

# Use custom configuration
python scripts/protein_refinement.py --input FILE --config configs/custom.json
```

### Demo Mode

All scripts support `--demo` mode with automatically generated test data:

```bash
python scripts/protein_refinement.py --demo --trajectories 2 --cycles 10
python scripts/protein_docking.py --demo --trajectories 3
python scripts/loop_modeling.py --demo --trajectories 3
python scripts/ligand_docking.py --demo --trajectories 3
python scripts/ddg_calculations.py --demo
```

### Configuration Files

Each script can use JSON configuration files:

```bash
python scripts/protein_refinement.py --input FILE --config configs/protein_refinement_config.json
```

Configuration files support:
- Parameter defaults
- Usage examples
- Documentation for each parameter
- Override capability via CLI arguments

### CLI Arguments Override

CLI arguments always override configuration file settings:

```bash
# Config says 10 trajectories, but CLI overrides to 5
python scripts/protein_refinement.py --config myconfig.json --trajectories 5
```

## Script Details

### 1. protein_refinement.py

**Purpose**: High-resolution structure refinement using Rosetta's Relax protocol

**Main Function**:
```python
run_protein_refinement(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]
```

**Key Parameters**:
- `trajectories`: Number of refinement trajectories (default: 5)
- `cycles`: Monte Carlo cycles per trajectory (default: 100)
- `temperature`: MC temperature in kT (default: 2.0)

**Output**: JSON with refinement scores and statistics

**Example**:
```bash
python scripts/protein_refinement.py \\
    --input examples/data/protein.pdb \\
    --output results/refined.json \\
    --trajectories 10 \\
    --cycles 200
```

### 2. protein_docking.py

**Purpose**: Protein-protein docking using Rosetta's DockingProtocol

**Main Function**:
```python
run_protein_docking(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]
```

**Key Parameters**:
- `trajectories`: Number of docking trajectories (default: 10)
- `chain_break`: Where chain A ends (auto-detect if None)

**Output**: JSON with docking scores and interface analysis

**Example**:
```bash
python scripts/protein_docking.py \\
    --input examples/data/complex.pdb \\
    --output results/docked.json \\
    --trajectories 20 \\
    --chain_break 150
```

### 3. loop_modeling.py

**Purpose**: Loop reconstruction using Rosetta's loop modeling protocol

**Main Function**:
```python
run_loop_modeling(
    input_file: Union[str, Path],
    loop_start: int,
    loop_end: int,
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]
```

**Key Parameters**:
- `loop_start`: Starting residue number (required)
- `loop_end`: Ending residue number (required)
- `trajectories`: Number of modeling trajectories (default: 10)

**Output**: JSON with loop scores and RMSD analysis

**Example**:
```bash
python scripts/loop_modeling.py \\
    --input examples/data/protein.pdb \\
    --loop_start 10 \\
    --loop_end 15 \\
    --output results/loops.json \\
    --trajectories 50
```

### 4. ligand_docking.py

**Purpose**: Protein-ligand docking using Rosetta's ligand docking protocol

**Main Function**:
```python
run_ligand_docking(
    input_file: Union[str, Path],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]
```

**Key Parameters**:
- `trajectories`: Number of docking trajectories (default: 10)
- `ligand_chain`: Chain containing ligand (auto-detect if None)
- `minimize`: Whether to minimize final complex (default: True)

**Output**: JSON with binding scores and ligand analysis

**Example**:
```bash
python scripts/ligand_docking.py \\
    --input examples/data/complex.pdb \\
    --output results/ligand_docked.json \\
    --ligand_chain L \\
    --trajectories 30
```

### 5. ddg_calculations.py

**Purpose**: ΔΔG calculations for protein stability analysis

**Main Function**:
```python
run_ddg_calculations(
    input_file: Union[str, Path],
    mutations: Union[str, List[Tuple[str, int, str]]],
    output_file: Optional[Union[str, Path]] = None,
    config: Optional[Dict[str, Any]] = None,
    **kwargs
) -> Dict[str, Any]
```

**Key Parameters**:
- `mutations`: Mutations like "A10G,B15L" (required)
- `repack_radius`: Radius for sidechain repacking (default: 8.0 Å)
- `minimize`: Whether to minimize after mutation (default: True)

**Output**: JSON with ΔΔG values for each mutation

**Example**:
```bash
python scripts/ddg_calculations.py \\
    --input examples/data/protein.pdb \\
    --mutations "A10G,A15L,B20V" \\
    --output results/ddg.json \\
    --repack_radius 10.0
```

## Shared Library

Common functions are in `scripts/lib/`:

### pyrosetta_utils.py
- `setup_pyrosetta()`: Lazy PyRosetta initialization with graceful error handling
- `create_demo_pdb()`: Generate minimal demo structures

### io_utils.py
- `validate_pdb_file()`: Basic PDB file validation
- `load_config()`: JSON configuration loading
- `save_results()`: Standardized results output

## Error Handling

All scripts provide graceful error handling:

1. **PyRosetta Not Available**: Clear installation instructions
2. **Invalid Input Files**: File validation with helpful messages
3. **Configuration Errors**: JSON parsing with line numbers
4. **Parameter Validation**: Range checking with suggestions

## Output Format

All scripts return standardized JSON output:

```json
{
  "result": {
    "original_score": -123.45,
    "best_score": -145.67,
    "score_improvement": 22.22,
    "trajectories_completed": 10,
    // ... protocol-specific results
  },
  "output_file": "path/to/results.json",
  "metadata": {
    "input_file": "path/to/input.pdb",
    "config": { /* full configuration used */ },
    "success": true
  }
}
```

On error:

```json
{
  "result": null,
  "output_file": null,
  "metadata": {
    "input_file": "path/to/input.pdb",
    "config": { /* configuration attempted */ },
    "success": false,
    "error": "Detailed error message"
  }
}
```

## For MCP Wrapping (Step 6)

Each script exports a main function that can be wrapped as an MCP tool:

```python
# Example MCP tool wrapper
from scripts.protein_refinement import run_protein_refinement

@mcp.tool()
def refine_protein_structure(input_file: str, output_file: str = None, trajectories: int = 5):
    """Refine protein structure using Rosetta Relax protocol."""
    return run_protein_refinement(
        input_file=input_file,
        output_file=output_file,
        trajectories=trajectories
    )
```

## Testing

All scripts have been tested with:
- Demo data generation ✓
- CLI argument parsing ✓
- Configuration file loading ✓
- Error handling ✓
- Help message generation ✓

## PyRosetta Dependency

**Important**: All scripts require PyRosetta which:
- Is not freely available via conda/pip
- Requires commercial license from RosettaCommons
- Must be installed manually: `conda install -c rosettacommons pyrosetta`

Scripts handle PyRosetta absence gracefully and provide clear installation instructions.

## Performance Notes

- **Startup Time**: Lazy loading minimizes import overhead when PyRosetta is not available
- **Memory Usage**: Demo structures are minimal to reduce memory footprint
- **Configuration**: JSON configs are cached to avoid repeated parsing
- **Error Recovery**: Scripts fail fast with clear messages to minimize debugging time

## File Structure

```
scripts/
├── lib/                    # Shared utilities
│   ├── __init__.py
│   ├── io_utils.py        # I/O and validation functions
│   └── pyrosetta_utils.py # PyRosetta interface functions
├── protein_refinement.py  # Structure refinement
├── protein_docking.py     # Protein-protein docking
├── loop_modeling.py       # Loop reconstruction
├── ligand_docking.py      # Protein-ligand docking
├── ddg_calculations.py    # ΔΔG stability analysis
└── README.md              # This file
```