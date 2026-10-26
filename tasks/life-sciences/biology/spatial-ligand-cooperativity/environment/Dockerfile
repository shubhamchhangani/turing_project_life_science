FROM python:3.11.9-slim-bookworm
RUN apt-get update && apt-get install -y --no-install-recommends bash gfortran gcc libblas-dev liblapack-dev pkg-config && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir numpy==1.26.4 scipy==1.13.1 pfapack==1.0.2
WORKDIR /root
COPY data /root/data
RUN mkdir -p /root/results
ENV OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1
