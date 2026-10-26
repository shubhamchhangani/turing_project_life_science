#!/bin/bash
set -euo pipefail
cd /root
python3 /solution/mfa_fit.py /root/data /root/results/fluxes.csv 24 7
