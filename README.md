# Rosetta MCP

**Protein modeling and design using the Rosetta suite via Docker**

An MCP (Model Context Protocol) server for molecular modeling with 7 core tools:
- Protein structure refinement using Relax protocol
- Mutation stability analysis (ΔΔG calculations)
- Submit protein-protein docking jobs with async tracking
- Submit loop modeling jobs
- Submit protein-ligand docking jobs
- Monitor and retrieve job results
- List available example structures

## Quick Start with Docker

### Approach 1: Pull Pre-built Image from GitHub

The fastest way to get started. A pre-built Docker image is automatically published to GitHub Container Registry on every release.

```bash
# Pull the latest image
docker pull ghcr.io/macromnex/rosetta_mcp:latest

# Register with Claude Code (runs as current user to avoid permission issues)
claude mcp add rosetta -- docker run -i --rm --user `id -u`:`id -g` -v `pwd`:`pwd` ghcr.io/macromnex/rosetta_mcp:latest
```

**Note:** Run from your project directory. `` `pwd` `` expands to the current working directory.

**Requirements:**
- Docker
- Claude Code installed

That's it! The Rosetta MCP server is now available in Claude Code.

---

### Approach 2: Build Docker Image Locally

Build the image yourself and install it into Claude Code. Useful for customization or offline environments.

```bash
# Clone the repository
git clone https://github.com/MacromNex/rosetta_mcp.git
cd rosetta_mcp

# Build the Docker image
docker build -t rosetta_mcp:latest .

# Register with Claude Code (runs as current user to avoid permission issues)
claude mcp add rosetta -- docker run -i --rm --user `id -u`:`id -g` -v `pwd`:`pwd` rosetta_mcp:latest
```

**Note:** Run from your project directory. `` `pwd` `` expands to the current working directory.

**Requirements:**
- Docker
- Claude Code installed
- Git (to clone the repository)

**About the Docker Flags:**
- `-i` — Interactive mode for Claude Code
- `--rm` — Automatically remove container after exit
- `` --user `id -u`:`id -g` `` — Runs the container as your current user, so output files are owned by you (not root)
- `-v` — Mounts your project directory so the container can access your data

---

## Verify Installation

After adding the MCP server, you can verify it's working:

```bash
# List registered MCP servers
claude mcp list

# You should see 'rosetta' in the output
```

In Claude Code, you can now use all 7 Rosetta tools:
- `refine_protein_structure`
- `calculate_ddg`
- `submit_protein_docking`
- `submit_loop_modeling`
- `submit_ligand_docking`
- `get_job_status`
- `get_job_result`

---

## Next Steps

- **Detailed documentation**: See [detail.md](detail.md) for comprehensive guides on:
  - Available MCP tools and parameters
  - Local Python environment setup (alternative to Docker)
  - Example workflows and use cases
  - PyRosetta license information
  - Troubleshooting

---

## Usage Examples

Once registered, you can use the Rosetta tools directly in Claude Code. Here are some common workflows:

### Example 1: Protein Structure Refinement

```
I have a protein structure at /path/to/protein.pdb. Can you use refine_protein_structure to refine it with 5 trajectories and 100 cycles, saving results to /path/to/results/?
```

### Example 2: Mutation Stability Analysis

```
I want to analyze how mutations A10G, A10L, and A10F affect the stability of /path/to/protein.pdb. Can you use calculate_ddg with 5 trajectories and report which mutations are stabilizing vs destabilizing?
```

### Example 3: Protein-Protein Docking

```
I have a protein complex at /path/to/complex.pdb with chain break at residue 150. Can you submit a docking job using submit_protein_docking with 20 trajectories, save results to /path/to/docking/, and monitor progress until completion?
```

---

## Troubleshooting

**Docker not found?**
```bash
docker --version  # Install Docker if missing
```

**Claude Code not found?**
```bash
# Install Claude Code
npm install -g @anthropic-ai/claude-code
```

**PyRosetta not available?**
- All tools work in demo mode without PyRosetta
- For full functionality, a PyRosetta academic or commercial license is required
- See [detail.md](detail.md) for PyRosetta installation instructions

---

## License

Rosetta License — Based on [Rosetta Commons](https://www.rosettacommons.org/) software suite. Academic use is free with registration.
