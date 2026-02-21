FROM condaforge/miniforge3:latest

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Create conda environment matching quick_setup.sh (Python 3.10)
RUN mamba create -p /app/env python=3.10 -y && mamba clean --all -y

# Install Python dependencies
COPY requirements.txt .
RUN /app/env/bin/pip install --no-cache-dir -r requirements.txt
# Force-reinstall fastmcp to match quick_setup.sh behavior
RUN /app/env/bin/pip install --force-reinstall --no-cache-dir fastmcp

# -------------------------------------------------------------------
# Optional: Install PyRosetta (requires license from RosettaCommons)
#
# Build with PyRosetta (academic):
#   docker build --build-arg PYROSETTA_CHANNEL=rosettacommons -t rosetta-mcp .
#
# Build with PyRosetta (commercial):
#   docker build --build-arg PYROSETTA_CHANNEL=rosettacommons/label/commercial -t rosetta-mcp .
# -------------------------------------------------------------------
ARG PYROSETTA_CHANNEL=""
RUN if [ -n "$PYROSETTA_CHANNEL" ]; then \
        mamba install -p /app/env -c "$PYROSETTA_CHANNEL" pyrosetta -y && \
        mamba clean --all -y && \
        echo "PyRosetta installed successfully"; \
    else \
        echo "Skipping PyRosetta (no PYROSETTA_CHANNEL specified)"; \
    fi

# Pre-initialize PyRosetta database if installed (saves startup time)
RUN if /app/env/bin/python -c "import pyrosetta" 2>/dev/null; then \
        /app/env/bin/python -c "\
import pyrosetta; \
pyrosetta.init(extra_options='-mute all'); \
print('PyRosetta initialized and database cached')" ; \
    else \
        echo "PyRosetta not installed, skipping initialization"; \
    fi

# Copy application source
COPY src/ ./src/
RUN chmod -R a+r /app/src/
COPY scripts/ ./scripts/
RUN chmod -R a+r /app/scripts/
COPY configs/ ./configs/
RUN chmod -R a+r /app/configs/
COPY examples/ ./examples/
RUN chmod -R a+r /app/examples/

# Create required runtime directories
RUN mkdir -p jobs tmp/inputs tmp/outputs

# Verify core imports work
RUN /app/env/bin/python -c "import fastmcp; import loguru; print('Core packages OK')"

ENV PYTHONPATH=/app

CMD ["/app/env/bin/python", "src/server.py"]
