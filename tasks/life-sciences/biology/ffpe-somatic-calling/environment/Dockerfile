FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends bash ca-certificates python3 python3-pip samtools bcftools tabix && \
    rm -rf /var/lib/apt/lists/*

# Pinned scientific Python stack available to the agent
RUN pip install --no-cache-dir --break-system-packages \
    numpy==1.26.4 \
    scipy==1.13.1 \
    pandas==2.2.2 \
    pysam==0.22.1

WORKDIR /root

# Agent-visible inputs only: reference genome and the two CRAM files with indexes.
# No truth tables, expected values or solution code are copied into this image.
COPY data /root/data
RUN mkdir -p /root/results
