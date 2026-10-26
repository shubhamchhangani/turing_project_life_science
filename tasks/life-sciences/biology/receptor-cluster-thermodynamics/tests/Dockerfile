FROM python:3.11.9-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends bash gfortran gcc libblas-dev liblapack-dev pkg-config && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==1.26.4 scipy==1.13.1 pfapack==1.0.2 pytest==8.4.1 pytest-json-ctrf==0.3.5
RUN mkdir -p /root/results /logs/verifier
COPY . /tests/
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
