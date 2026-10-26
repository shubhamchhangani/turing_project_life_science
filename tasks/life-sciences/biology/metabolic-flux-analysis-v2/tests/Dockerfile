FROM python:3.11.9-slim-bookworm

RUN apt-get update && \
    apt-get install -y --no-install-recommends bash && \
    rm -rf /var/lib/apt/lists/*

# Pre-install test dependencies (no runtime installs).
# pytest and pytest-json-ctrf are pinned benchmark-wide to these exact versions.
RUN pip install --no-cache-dir pytest==8.4.1 pytest-json-ctrf==0.3.5 numpy==1.26.4 scipy==1.13.1 numba==0.60.0

# Pre-create parent directories for declared artifacts and verifier logs
RUN mkdir -p /root/results /logs/verifier

# Bake the test scripts and expected values into the verifier image
COPY . /tests/

ENV NUMBA_CACHE_DIR=/tmp/numba-cache OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
