# MCP Rosetta Configuration Files

Configuration files for MCP Rosetta scripts providing parameter defaults, documentation, and usage examples.

## Overview

| Config File | Script | Description |
|-------------|--------|-------------|
| `default_config.json` | All | Default template with global settings |
| `protein_refinement_config.json` | `protein_refinement.py` | Structure refinement parameters |
| `protein_docking_config.json` | `protein_docking.py` | Protein-protein docking parameters |
| `loop_modeling_config.json` | `loop_modeling.py` | Loop reconstruction parameters |
| `ligand_docking_config.json` | `ligand_docking.py` | Protein-ligand docking parameters |
| `ddg_calculations_config.json` | `ddg_calculations.py` | ΔΔG calculation parameters |

## Usage

### Basic Usage

```bash
# Use script-specific config
python scripts/protein_refinement.py --input FILE --config configs/protein_refinement_config.json

# Override specific parameters
python scripts/protein_refinement.py --config myconfig.json --trajectories 20
```

### Creating Custom Configs

```bash
# Copy default config
cp configs/protein_refinement_config.json myconfig.json

# Edit parameters
{
  "trajectories": 20,
  "cycles": 500,
  "temperature": 1.5
}
```

## Config File Structure

All config files follow this structure:

```json
{
  "_description": "Human-readable description of the configuration",
  "_source": "Original source file this config was extracted from",
  "_script": "Target script file that uses this config",

  "parameter1": "value1",
  "parameter2": 123,
  "parameter3": true,

  "_parameter_descriptions": {
    "parameter1": "Description of what parameter1 does",
    "parameter2": "Description of parameter2 with units/range",
    "parameter3": "Boolean parameter description"
  },

  "_usage_examples": {
    "example_name": {
      "parameter1": "example_value",
      "parameter2": 456,
      "description": "What this example configuration achieves"
    }
  }
}
```

## Configuration Parameters

### Global Parameters (All Scripts)

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `pymol_visualization` | bool | false | Enable PyMOL visualization during protocols |
| `output_prefix` | str | varies | Prefix for output files |
| `temperature` | float | 2.0 | Monte Carlo temperature (kT) |

### Script-Specific Parameters

#### protein_refinement.py

| Parameter | Type | Default | Description | Range |
|-----------|------|---------|-------------|-------|
| `trajectories` | int | 5 | Number of refinement trajectories | 1-100 |
| `cycles` | int | 100 | Monte Carlo cycles per trajectory | 10-1000 |

**Usage Examples**:
- `quick_refinement`: Fast testing (2 trajectories, 50 cycles)
- `thorough_refinement`: Publication quality (20 trajectories, 500 cycles)
- `high_temperature`: Escape local minima (temperature=5.0)

#### protein_docking.py

| Parameter | Type | Default | Description | Range |
|-----------|------|---------|-------------|-------|
| `trajectories` | int | 10 | Number of docking trajectories | 1-100 |
| `chain_break` | int/null | null | Where chain A ends (auto-detect) | 1-N_residues |
| `use_centroid_stage` | bool | true | Use low-resolution docking stage | - |
| `use_fullatom_stage` | bool | true | Use high-resolution docking stage | - |

**Usage Examples**:
- `quick_docking`: Fast mode (5 trajectories, fullatom only)
- `thorough_docking`: Many trajectories (50)
- `two_stage_docking`: Full protocol (centroid + fullatom)
- `manual_chain_break`: Specify chain boundary

#### loop_modeling.py

| Parameter | Type | Default | Description | Range |
|-----------|------|---------|-------------|-------|
| `loop_start` | int | null | Loop start residue (required) | 1-N_residues |
| `loop_end` | int | null | Loop end residue (required) | loop_start+1 to N_residues |
| `trajectories` | int | 10 | Number of modeling trajectories | 1-200 |
| `cycles_low` | int | 100 | Centroid stage cycles | 10-1000 |
| `cycles_high` | int | 50 | Fullatom stage cycles | 10-500 |

**Usage Examples**:
- `short_loop`: 6 residues (10-15), 20 trajectories
- `long_loop`: 16 residues (50-65), 50 trajectories, more cycles
- `quick_modeling`: Fast testing (5 trajectories, reduced cycles)
- `thorough_modeling`: Difficult cases (100 trajectories, many cycles)

#### ligand_docking.py

| Parameter | Type | Default | Description | Range |
|-----------|------|---------|-------------|-------|
| `trajectories` | int | 10 | Number of docking trajectories | 1-100 |
| `ligand_chain` | str/null | null | Chain ID with ligand (auto-detect) | A-Z |
| `high_res_docking` | bool | true | Use high-resolution docking | - |
| `repack_sidechains` | bool | true | Repack protein sidechains | - |
| `minimize` | bool | true | Minimize final complex | - |

**Usage Examples**:
- `quick_docking`: Fast mode (5 trajectories, no optimization)
- `thorough_docking`: Full optimization (50 trajectories)
- `manual_ligand_chain`: Specify ligand chain ID
- `rigid_docking`: Rigid body only (no high-res stage)

#### ddg_calculations.py

| Parameter | Type | Default | Description | Range |
|-----------|------|---------|-------------|-------|
| `mutations` | list | [] | List of mutations (CLI input preferred) | - |
| `repack_radius` | float | 8.0 | Repacking radius around mutations (Å) | 4.0-15.0 |
| `minimize` | bool | true | Minimize after mutation | - |
| `use_cartesian` | bool | false | Use Cartesian minimization | - |

**Usage Examples**:
- `single_mutation`: One mutation for testing
- `multiple_mutations`: Several mutations for analysis
- `large_repack`: Bigger radius for buried mutations
- `no_minimize`: Skip minimization for speed
- `cartesian_minimize`: Better accuracy with Cartesian
- `interface_mutations`: Mutations at protein interfaces

## Environment-Specific Configs

### Development/Testing

```json
{
  "trajectories": 2,
  "cycles": 10,
  "temperature": 2.0,
  "_description": "Fast parameters for development and testing"
}
```

### Production/Research

```json
{
  "trajectories": 50,
  "cycles": 500,
  "temperature": 1.0,
  "_description": "Thorough parameters for production analysis"
}
```

### High-Throughput

```json
{
  "trajectories": 1,
  "cycles": 50,
  "minimize": false,
  "_description": "Minimal parameters for high-throughput screening"
}
```

## Parameter Validation

Scripts validate parameters and provide helpful error messages:

```bash
# Invalid loop range
python scripts/loop_modeling.py --loop_start 10 --loop_end 8
# Error: Loop end (8) must be greater than loop start (10)

# Temperature out of range
python scripts/protein_refinement.py --temperature -1.0
# Warning: Negative temperature may cause issues
```

## Config File Loading

1. **Default values** loaded from script's DEFAULT_CONFIG
2. **Config file values** override defaults if specified
3. **CLI arguments** override config file values
4. **Parameter validation** ensures values are in acceptable ranges

Priority: `CLI args > Config file > Script defaults`

## Best Practices

### 1. Use Descriptive Names

```json
{
  "_description": "High-accuracy refinement for NMR structures"
}
```

### 2. Document Parameter Choices

```json
{
  "trajectories": 100,
  "_comment": "Increased from 50 based on convergence analysis"
}
```

### 3. Include Usage Context

```json
{
  "_usage_context": "For membrane proteins with low initial quality",
  "temperature": 3.0,
  "cycles": 1000
}
```

### 4. Version Your Configs

```json
{
  "_version": "1.2.0",
  "_date": "2025-12-15",
  "_author": "research_team"
}
```

## Common Config Patterns

### Quick Testing

```json
{
  "trajectories": 2,
  "cycles": 10,
  "minimize": false,
  "_description": "Minimal parameters for fast testing"
}
```

### Production Analysis

```json
{
  "trajectories": 50,
  "cycles": 500,
  "temperature": 1.0,
  "_description": "Thorough parameters for publication"
}
```

### High Throughput

```json
{
  "trajectories": 1,
  "cycles": 25,
  "repack_sidechains": false,
  "_description": "Fast screening parameters"
}
```

### Difficult Cases

```json
{
  "trajectories": 200,
  "cycles": 2000,
  "temperature": 4.0,
  "_description": "Extensive sampling for challenging targets"
}
```

## Error Handling

### Invalid JSON

```bash
# Malformed JSON file
python scripts/protein_refinement.py --config bad.json
# Error: Invalid JSON in config file: Expecting ',' delimiter: line 5 column 2
```

### Missing Required Parameters

```bash
# Missing loop_start in config
python scripts/loop_modeling.py --config incomplete.json
# Error: Loop start required! Use --loop_start or --demo
```

### Parameter Type Errors

```bash
# String instead of number
python scripts/protein_refinement.py --trajectories "many"
# Error: argument --trajectories: invalid int value: 'many'
```

## Integration with MCP Tools

Config files are designed for easy MCP tool integration:

```python
# MCP tool with config support
@mcp.tool()
def refine_protein(input_file: str, config_name: str = "default") -> dict:
    config_path = f"configs/{config_name}_config.json"
    return run_protein_refinement(input_file=input_file, config=load_config(config_path))
```

## File Structure

```
configs/
├── default_config.json              # Global defaults and template
├── protein_refinement_config.json   # Structure refinement
├── protein_docking_config.json      # Protein-protein docking
├── loop_modeling_config.json        # Loop reconstruction
├── ligand_docking_config.json       # Protein-ligand docking
├── ddg_calculations_config.json     # ΔΔG calculations
└── README.md                        # This file
```