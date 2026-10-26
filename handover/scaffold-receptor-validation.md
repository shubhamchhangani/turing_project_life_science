# Scaffold receptor candidate — rejected after Claude solve

Claude Opus 5.5 (low effort) earned reward 1 with all four verifier tests passing.
Agent elapsed time: 869.35558 seconds; 14 completed commands; no infrastructure
error or quota interruption. GPT was not launched and no ZIP was created.

The trajectory contains an exact planar dimer/Kasteleyn construction, recursive
Pfaffian elimination, a Walsh transform and reweighting of all nonlocal couplings
and fields. The solver recovered from a memory failure by chunking the final
sum. It also tested its implementation against exhaustive enumeration on a
small lattice. The solver reports approximately 25 seconds for the final numeric
run; 869.35558 seconds is the independently logged full agent time.

This is a substantive solve. Enlarging this same family is not evidence of a
harder reasoning benchmark. Preserve this package unchanged as a rejected trial.

Validation completed before the model trial:
- Reference: 4703.688881181 seconds in one CPU, 2 GiB, offline.
- Separate source implementation: 3411.5270232500043 seconds.
- Maximum probability disagreement: 1.2896628209801975e-13.
- Log-partition disagreement: 4.547473508864641e-13.
- Log-clamp disagreement: 1.4210854715202004e-13.
- Twelve exhaustive small-model comparisons and 1788 additional moment checks.
- Three valid and 280 invalid-output regression cases.
- Four separate offline verifier tests and static package checks passed.

Package fingerprint: 808e620476918bc4f8265255b881f43ea1341031772d6dbb62b6eecf724eb1f1
Trial: jobs/receptor-torus-frontier/opus-5-5-low-scaffold-receptor/scaffold-receptor-ensemble__NgbJBw8/result.json

No result here guarantees Turing portal acceptance or failure by other models.
