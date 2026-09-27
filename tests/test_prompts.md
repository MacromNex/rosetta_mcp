# Comprehensive Test Prompts for Rosetta MCP Server

## Tool Discovery Tests

### Prompt 1: List All Tools
"What MCP tools are available from the rosetta server? Give me a brief description of each."

### Prompt 2: Tool Details
"Explain how to use the refine_protein_structure tool, including all parameters."

### Prompt 3: Job Management Tools
"List all job management tools and explain how they work together for async operations."

## Sync Tool Tests

### Prompt 4: Basic Structure Refinement
"Use refine_protein_structure with input_file='examples/data/test_input.pdb' and default parameters."

### Prompt 5: Custom Refinement Parameters
"Run refine_protein_structure on examples/data/test_input.pdb with trajectories=3, cycles=50, and temperature=1.5"

### Prompt 6: DDG Calculation
"Calculate ΔΔG for mutation A10G using examples/data/test_input.pdb"

### Prompt 7: Structure Validation
"Validate the PDB structure examples/data/test_input.pdb"

### Prompt 8: List Example Structures
"List all available example structures for testing."

### Prompt 9: Error Handling - Invalid File
"Try running refine_protein_structure with a non-existent file '/fake/path.pdb'"

### Prompt 10: Error Handling - Invalid Parameters
"Run calculate_ddg with mutations='INVALID' on examples/data/test_input.pdb"

## Submit API Tests

### Prompt 11: Submit Protein Docking
"Submit a protein-protein docking job for examples/data/test_docking.pdb"

### Prompt 12: Submit Loop Modeling
"Submit loop modeling for examples/data/test_loop.pdb with loop_start=5 and loop_end=15"

### Prompt 13: Submit Ligand Docking
"Submit protein-ligand docking for examples/data/test_ligand.pdb"

### Prompt 14: Submit Large Refinement
"Submit a large-scale refinement with 20 trajectories for examples/data/test_input.pdb"

### Prompt 15: Check Job Status
"What's the status of job {job_id}?"

### Prompt 16: Get Job Results
"Show me the results of job {job_id}"

### Prompt 17: View Job Logs
"Show the last 30 lines of logs for job {job_id}"

### Prompt 18: List All Jobs
"List all jobs with status 'completed'"

### Prompt 19: Cancel Job
"Cancel the running job {job_id}"

## Batch Processing Tests

### Prompt 20: Batch Submit
"Process multiple files in batch: examples/data/test_input.pdb, examples/data/test_complex.pdb"

### Prompt 21: Batch Results
"Get all results from batch job {batch_job_id}"

## End-to-End Scenarios

### Prompt 22: Full Analysis Workflow
"Analyze the protein in examples/data/test_input.pdb:
1. First validate its structure
2. Then refine it with 3 trajectories
3. Finally calculate ΔΔG for mutation A10G"

### Prompt 23: Conditional Processing
"If the structure in examples/data/test_input.pdb has more than 100 residues,
submit it for large refinement with 10 trajectories. Otherwise, refine it directly with default parameters."

### Prompt 24: Error Recovery
"Submit a loop modeling job. If it fails, show me the error log
and suggest what might be wrong."

### Prompt 25: Complex Analysis
"Compare two refinement approaches for examples/data/test_input.pdb:
1. Quick refinement with 3 trajectories and 50 cycles
2. Submit a large refinement with 20 trajectories and 500 cycles
Show the differences in results when both complete."

### Prompt 26: Multi-Step Pipeline
"Create a complete analysis pipeline:
1. Validate examples/data/test_complex.pdb
2. If valid, submit protein docking
3. While docking runs, calculate ΔΔG for mutation A15G
4. When docking completes, show both results"

### Prompt 27: Resource Management
"List all currently running jobs, cancel any that have been running for more than 10 minutes,
and show me the logs of the most recent completed job."

### Prompt 28: Data Comparison
"For each PDB file in examples/data/, validate the structure and create a summary table
showing file name, atom count, chain count, and validation status."

## Advanced Test Scenarios

### Prompt 29: Stress Test
"Submit 5 different jobs in parallel:
- Protein docking for test_docking.pdb
- Loop modeling for test_loop.pdb (loop 5-15)
- Ligand docking for test_ligand.pdb
- Large refinement for test_input.pdb
- Batch refinement for test_complex.pdb
Monitor all jobs until completion."

### Prompt 30: Error Handling Comprehensive
"Test error handling by:
1. Trying to get status of non-existent job 'fake123'
2. Attempting to cancel a completed job
3. Getting results of a pending job
4. Viewing logs of invalid job ID
Show how the system handles each error case."