FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libcairo2 \
    libpango-1.0-0 \
    libpangocairo-1.0-0 \
    curl \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install optional chemistry tools (non-fatal)
RUN apt-get update && apt-get install -y --no-install-recommends \
    openbabel libopenbabel-dev \
    2>/dev/null; rm -rf /var/lib/apt/lists/* || true

# Try to install autodock-vina (non-fatal if not in apt)
# Note: fpocket is not in Debian repos; pipeline uses P2Rank or geometric fallback
RUN apt-get update && apt-get install -y --no-install-recommends \
    autodock-vina \
    2>/dev/null; rm -rf /var/lib/apt/lists/* || true

# V3: Install smina for rapid scoring (Vinardo scoring function)
RUN wget -q https://sourceforge.net/projects/smina/files/smina.static/download -O /usr/local/bin/smina \
    && chmod +x /usr/local/bin/smina \
    || echo "WARNING: smina not available, using mock fallback"

# V5bis: Install Java JRE for P2Rank
RUN apt-get update && apt-get install -y --no-install-recommends \
    default-jre-headless \
    2>/dev/null; rm -rf /var/lib/apt/lists/* || true

# V5bis: Install P2Rank 2.4.2 (ML-based pocket detection)
RUN mkdir -p /opt && \
    wget -q https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz -O /tmp/p2rank.tar.gz && \
    tar -xzf /tmp/p2rank.tar.gz -C /opt/ && \
    rm /tmp/p2rank.tar.gz && \
    chmod +x /opt/p2rank_2.4.2/prank \
    || echo "WARNING: P2Rank not available, using mock fallback"

# V5bis: Install GNINA binary (CNN-scored docking)
RUN wget -q https://github.com/gnina/gnina/releases/download/v1.1/gnina -O /usr/local/bin/gnina \
    && chmod +x /usr/local/bin/gnina \
    || echo "WARNING: GNINA not available, using Vina/mock fallback"

WORKDIR /app

# Install core Python dependencies
COPY requirements.txt .

# Pre-install CPU-only PyTorch BEFORE other deps to avoid pulling 16GB CUDA version
# (admet-ai depends on torch, but this project runs without GPU)
RUN pip install --no-cache-dir torch torchvision --index-url https://download.pytorch.org/whl/cpu

RUN pip install --no-cache-dir -r requirements.txt

# Install optional packages (may fail on some platforms)
RUN pip install --no-cache-dir rdkit-pypi==2023.9.3 || \
    pip install --no-cache-dir rdkit || \
    echo "WARNING: RDKit not available, using fallbacks"

RUN pip install --no-cache-dir meeko || \
    echo "WARNING: Meeko not available, using fallbacks"

RUN pip install --no-cache-dir cairosvg==2.7.1 || \
    echo "WARNING: CairoSVG not available"

# V2: Install optional AI/chemistry packages (non-fatal)
RUN pip install --no-cache-dir admet-ai 2>/dev/null || \
    echo "WARNING: admet-ai not available, using mock fallback"

RUN pip install --no-cache-dir aizynthfinder 2>/dev/null || \
    echo "WARNING: aizynthfinder not available, using mock fallback"

# V5bis: Install ProLIF + MDAnalysis for interaction analysis (non-fatal)
RUN pip install --no-cache-dir prolif MDAnalysis 2>/dev/null || \
    echo "WARNING: ProLIF/MDAnalysis not available, using mock fallback"

# Copy application code
COPY . .

# Create data directories
RUN mkdir -p /data/jobs /data/zinc /data/aizynthfinder /data/chembl /data/off_targets

EXPOSE 8000

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "8000"]
