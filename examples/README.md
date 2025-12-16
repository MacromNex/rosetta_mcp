# Rosetta MCP Examples

This directory contains example scripts demonstrating various Rosetta capabilities through PyRosetta.

## Available Use Cases

### 1. Protein Structure Refinement (`use_case_1_protein_refinement.py`)
- **Category**: Relax, Structure quality analysis
- **Description**: High-resolution structure refinement using small backbone perturbations, sidechain packing, and minimization
- **Usage**: `python examples/use_case_1_protein_refinement.py --input examples/data/test_in.pdb --demo`
- **Demo command**: `python examples/use_case_1_protein_refinement.py --demo`

### 2. Protein-Protein Docking (`use_case_2_protein_docking.py`)
- **Category**: Protein-protein docking
- **Description**: Rigid-body docking protocol with low-resolution sampling and high-resolution refinement
- **Usage**: `python examples/use_case_2_protein_docking.py --input examples/data/complex.pdb --chain_break 100`
- **Demo command**: `python examples/use_case_2_protein_docking.py --demo`

### 3. Loop Modeling (`use_case_3_loop_modeling.py`)
- **Category**: Loop modeling
- **Description**: Models protein loops using CCD loop closure and fragment insertion
- **Usage**: `python examples/use_case_3_loop_modeling.py --input examples/data/test_in.pdb --loop_start 50 --loop_end 55`
- **Demo command**: `python examples/use_case_3_loop_modeling.py --demo`

### 4. Ligand Docking (`use_case_4_ligand_docking.py`)
- **Category**: Ligand docking
- **Description**: High-resolution protein-ligand docking with rigid-body sampling and interface optimization
- **Usage**: `python examples/use_case_4_ligand_docking.py --input examples/data/protein_ligand.pdb --ligand_chain B`
- **Demo command**: `python examples/use_case_4_ligand_docking.py --demo`

### 5. ddG Calculations (`use_case_5_ddg_calculations.py`)
- **Category**: ddG calculations, Structure quality analysis
- **Description**: Alanine scanning mutagenesis to calculate binding energy changes
- **Usage**: `python examples/use_case_5_ddg_calculations.py --input examples/data/complex.pdb --chain1 A --chain2 B`
- **Demo command**: `python examples/use_case_5_ddg_calculations.py --demo`

## Demo Data

### PDB Files
- `test_in.pdb` - General test structure for refinement and loop modeling
- `design_in.pdb` - Test structure for protein design applications

### Parameter Files
- `ZN1.params` - Zinc parameter file for ligand applications

## Running Examples

### Prerequisites
1. Activate the conda environment:
   ```bash
   mamba activate ./env  # or: conda activate ./env
   ```

2. Ensure PyRosetta is installed:
   ```bash
   # PyRosetta requires special licensing - install according to RosettaCommons instructions
   # Alternatively, these scripts will gracefully handle missing PyRosetta
   ```

### Quick Start with Demo Data
Each script has a `--demo` flag that creates minimal test data and runs with reduced parameters:

```bash
# Run all demos
python examples/use_case_1_protein_refinement.py --demo
python examples/use_case_2_protein_docking.py --demo
python examples/use_case_3_loop_modeling.py --demo
python examples/use_case_4_ligand_docking.py --demo
python examples/use_case_5_ddg_calculations.py --demo
```

### Using Your Own Data
```bash
# Protein refinement
python examples/use_case_1_protein_refinement.py --input your_protein.pdb --trajectories 10

# Protein docking (two-chain complex)
python examples/use_case_2_protein_docking.py --input your_complex.pdb --chain_break 150

# Loop modeling
python examples/use_case_3_loop_modeling.py --input your_protein.pdb --loop_start 45 --loop_end 55

# Ligand docking
python examples/use_case_4_ligand_docking.py --input your_complex.pdb --ligand_chain L --protein_chain P

# ddG calculations
python examples/use_case_5_ddg_calculations.py --input your_complex.pdb --chain1 A --chain2 B
```

## Output Files

Each script generates output files with descriptive prefixes:
- **Refinement**: `refined_*.pdb` - Refined structures
- **Docking**: `docked_*.pdb` - Docked conformations
- **Loop modeling**: `loop_model_*.pdb` - Loop models
- **Ligand docking**: `ligand_dock_*.pdb` - Ligand poses
- **ddG calculations**: `ddg_results.txt` - Mutation energy changes

## Common Options

All scripts support these common flags:
- `--input, -i`: Input PDB file
- `--output, -o`: Output prefix
- `--trajectories, -n`: Number of sampling trajectories
- `--demo`: Run with demo data and reduced parameters
- `--help`: Show full usage information

## Troubleshooting

1. **PyRosetta not found**: Install PyRosetta from RosettaCommons or run in demo mode to test script structure
2. **Memory issues**: Reduce `--trajectories` parameter
3. **Long runtime**: Use `--demo` for quick testing

## File Structure

```
examples/
├── README.md                          # This file
├── use_case_1_protein_refinement.py   # Structure refinement
├── use_case_2_protein_docking.py      # Protein-protein docking
├── use_case_3_loop_modeling.py        # Loop modeling
├── use_case_4_ligand_docking.py       # Ligand docking
├── use_case_5_ddg_calculations.py     # ddG calculations
├── data/                              # Demo data
│   ├── test_in.pdb                    # Test protein structure
│   ├── design_in.pdb                  # Design test structure
│   └── ZN1.params                     # Zinc parameters
└── models/                            # Pre-trained models (if any)
```