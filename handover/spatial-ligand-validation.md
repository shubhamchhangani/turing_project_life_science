# Spatial ligand cooperativity: rejected after a successful model solution

Candidate directory:
`tasks/life-sciences/biology/spatial-ligand-cooperativity/`.
This original synthetic biophysics task recovers 646 receptor interactions
and twelve local ligand biases from equilibrium observations. It is not a
claim to have discovered a previously unsolved biological problem.

## Completed validation

- Schema 1.4, allowed root contents, no canary/hidden/cache files, shell
  errexit, artifact paths and documented tests: 34/34 static checks passed.
- Every one of the 658 parameters and every supplied moment is checked.
- Exhaustive small-system validation includes all twelve fields; maximum
  discrepancy was 1.32e-13. Full-size independent agreement was 2.85e-13.
- Input-only packaged reference, offline with one CPU and 2 GiB:
  745.39 seconds, 56 function evaluations.
- Maximum reference recovery errors: 9.192e-7 for couplings and 7.96e-8
  for biases, versus the stated 5e-5 tolerance.
- Separate offline verifier: all four tests passed, reward 1, 16.31 seconds.
- Four valid output variants accepted; 34 invalid variants rejected,
  including whitespace numerics, missing/duplicate keys, the twelfth bias,
  nonregular files and a collective fit error within parameter tolerances.
- The automatic release conditions passed eleven targeted checks. Quota
  interruptions, incomplete runs, errors, model mismatches, changed task
  snapshots, machine-sleep interruptions and solved candidates cannot
  produce a release ZIP.

## Completed difficulty evaluation

Claude Opus 5.5, low effort, solved this candidate: reward 1 and all four
verifier tests passed. Agent execution took 933.908309 seconds. The trial
completed on 24 October 2026 at 04:03 UTC. No exception or quota interruption
was reported. The queued GPT trial was not launched; the worker exited after
the successful Claude result. No ZIP was created.

The trajectory shows a genuine input-only solution: a four-sector Pfaffian
calculation, an expansion over the 2048 even subsets of the twelve biased
sites, low-rank inverse updates, and bounded L-BFGS moment matching. Claude
checked its evaluator against exhaustive enumeration on small lattices.
Its submitted moments have maximum errors of approximately 4.8e-8 (edges)
and 2.3e-8 (sites). Increasing the count of local biases did not establish
sufficient difficulty for this model.

Evidence:
`jobs/receptor-torus-frontier/opus-5-5-low-spatial-ligand/spatial-ligand-cooperativity__pEUBHq4/`.
The frozen package manifest remains unchanged in the author records. This
candidate must not be released as an all-models-fail benchmark.

## Earlier rejected candidates

The three-bias inverse task was solved by GPT-6-sol: reward 1, four tests
passed, 1358.61 seconds of agent execution. Its Opus trial was cancelled
afterward and does not count as a model failure. The earlier forward
receptor and metabolic-flux candidates were also solved and remain rejected.

## Evidence

`tools/receptor-torus/draft/twelve-fields/` contains the generation provenance,
independent physics and enumeration checks, artifact regressions, offline
validation, verifier logs, frozen package manifest and release-condition
checks. These author records are outside the allowed task root and will not
be included in the submission ZIP.
