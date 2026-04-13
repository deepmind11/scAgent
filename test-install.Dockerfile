FROM python:3.12-slim

# Install Node.js 22 + build tools (g++ needed for annoy/scrublet)
RUN apt-get update && apt-get install -y curl git g++ && \
    curl -fsSL https://deb.nodesource.com/setup_22.x | bash - && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/*

# Verify prerequisites
RUN python3 --version && node --version && npm --version

# Clone from GitHub (as a real user would)
RUN git clone https://github.com/deepmind11/scAgent.git /app
WORKDIR /app

# Install exactly as README says
RUN python -m venv .venv && \
    . .venv/bin/activate && \
    pip install -e .

# Pre-install pi-coding-agent so first launch is fast
RUN npm install -g @mariozechner/pi-coding-agent

CMD ["/bin/bash", "-c", "source /app/.venv/bin/activate && exec bash"]
