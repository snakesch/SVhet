# syntax=docker/dockerfile:1.4
FROM python:3.12.6-slim-bookworm

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /opt

# Install system dependencies in one layer
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y \
    build-essential \
    wget \
    curl \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    autoconf \
    automake \
    bzip2 \
    ca-certificates \
    tabix

# Install htslib, bcftools, samtools, and bedtools in parallel using multi-stage
FROM python:3.12.6-slim-bookworm AS builder
WORKDIR /build

RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y \
    build-essential wget zlib1g-dev libbz2-dev liblzma-dev \
    libncurses5-dev libcurl4-openssl-dev libssl-dev autoconf automake

# Build htslib
RUN --mount=type=cache,target=/build/cache \
    wget -q https://github.com/samtools/htslib/releases/download/1.22/htslib-1.22.tar.bz2 && \
    tar -xjf htslib-1.22.tar.bz2 && \
    cd htslib-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf htslib-1.22*

# Build bcftools
RUN --mount=type=cache,target=/build/cache \
    wget -q https://github.com/samtools/bcftools/releases/download/1.22/bcftools-1.22.tar.bz2 && \
    tar -xjf bcftools-1.22.tar.bz2 && \
    cd bcftools-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf bcftools-1.22*

# Build samtools
RUN --mount=type=cache,target=/build/cache \
    wget -q https://github.com/samtools/samtools/releases/download/1.22/samtools-1.22.tar.bz2 && \
    tar -xjf samtools-1.22.tar.bz2 && \
    cd samtools-1.22 && \
    ./configure --prefix=/usr/local && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf samtools-1.22*

# Build bedtools
RUN --mount=type=cache,target=/build/cache \
    wget -q https://github.com/arq5x/bedtools2/releases/download/v2.31.1/bedtools-2.31.1.tar.gz && \
    tar -xzf bedtools-2.31.1.tar.gz && \
    cd bedtools2 && \
    make -j$(nproc) && \
    make install && \
    cd .. && rm -rf bedtools-2.31.1* bedtools2

# Final stage
FROM python:3.12.6-slim-bookworm
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Hong_Kong
WORKDIR /opt

# Install only runtime dependencies
RUN --mount=type=cache,target=/var/cache/apt,sharing=locked \
    --mount=type=cache,target=/var/lib/apt,sharing=locked \
    apt-get update && apt-get install -y \
    libcurl4 \
    libssl3 \
    libncurses6 \
    zlib1g \
    libbz2-1.0 \
    liblzma5 \
    git \
    tabix \
    && rm -rf /var/lib/apt/lists/*

# Copy compiled binaries from builder
COPY --from=builder /usr/local/bin/* /usr/local/bin/
COPY --from=builder /usr/local/lib/* /usr/local/lib/
COPY --from=builder /usr/local/include/* /usr/local/include/

# Update library cache
RUN ldconfig

# Install Python packages with cache
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --no-cache-dir \
    numpy==1.26.4 \
    pysam==0.22.1

# Clone SVhet
RUN git clone https://github.com/snakesch/SVhet.git && \
    cd SVhet && \
    chmod +x svhet.py

ENV PATH="/opt/SVhet:${PATH}"

# Create test scripts (same as before)
RUN echo '#!/bin/bash\n\
set -e\n\
echo "SVhet Installation Test"\n\
echo "========================"\n\
python --version\n\
python -c "import numpy; print(f\"Numpy: {numpy.__version__}\")" \n\
python -c "import pysam; print(f\"Pysam: {pysam.__version__}\")" \n\
bcftools --version | head -n1\n\
bedtools --version\n\
samtools --version | head -n1\n\
python /opt/SVhet/svhet.py --help > /dev/null\n\
echo "✓ All tests passed!"\n\
' > /opt/test_installation.sh && chmod +x /opt/test_installation.sh

WORKDIR /opt/SVhet
RUN /opt/test_installation.sh

CMD ["/bin/bash"]