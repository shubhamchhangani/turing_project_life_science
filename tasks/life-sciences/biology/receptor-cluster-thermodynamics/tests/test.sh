#!/bin/bash
set -euo pipefail
mkdir -p /logs/verifier
if timeout --kill-after=5s 840s pytest --ctrf /logs/verifier/ctrf.json /tests/test_outputs.py -rA; then
    echo 1 > /logs/verifier/reward.txt
else
    echo 0 > /logs/verifier/reward.txt
fi
