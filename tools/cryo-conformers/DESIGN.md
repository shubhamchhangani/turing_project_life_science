# Unknown-orientation structural reconstruction screening

Status: exploratory, not a task package; no ZIP and no frontier-model trials.

The previously tested polyploid-phase-posterior task was rejected for release:
Claude Opus solved it with reward 1 in 84.37186 seconds. This direction is a
separate synthetic inverse-problem experiment, not a claim to solve an open
biological research problem.

## Completed screen

`screen.py` constructs 33-cubed densities from 18 Gaussian components and
generates 128 clean random projections. The estimator gets an ArrayImageSource
containing pixels only, without simulation metadata or true orientations.
Ground-truth rotations are accessed afterwards solely to evaluate the estimate.
The seed is fixed. Dependencies are isolated in tools/.venv-cryo (ASPIRE 0.14.3).

On the local host:

| Density states | Orientation-estimation seconds | Registered rotation matrix MSE |
| --- | ---: | ---: |
| 1 | 1.18599 | 0.0000156388 |
| 2 | 0.85242 | 0.0941240 |

The second density displaces seven components by (3, -2, 1) voxels.
Both cases use the unmodified single-state CLSyncVoting baseline with 180
angular samples. This measures orientation recovery only, not successful
volume reconstruction or full task completion. Host timings do not establish
one-CPU container resource compliance. The two-state error does not establish
failure of an appropriately designed heterogeneous reconstruction method.

## Decision and outstanding evidence

The clean single-state orientation problem has a strong off-the-shelf shortcut;
do not release it as a demonstrated hard benchmark. The heterogeneous case has
not yet earned a difficulty claim. Before constructing a package it needs:

1. A complete input-only reference reconstruction with no true labels or poses.
2. Identifiability checks and acceptance invariant to equivalent global frames,
   reflections and state permutations, with full structural coverage.
3. Independent generation/provenance and calibrated scientific tolerances.
4. Resource-limited reference validation, adversarial verifier checks, and the
   previously required schema, shell, decimal and direct-root package audits.
5. Fresh authorized, isolated roster trials; any successful solve disqualifies
   the candidate from this user's requested release gate.

No all-model-failure guarantee can be inferred from these screens.

## Reproduction

From the workspace root, set ASPIREDIR to a writable scratch configuration
directory and MPLCONFIGDIR to a writable font-cache directory, then run:

```
tools/.venv-cryo/bin/python tools/cryo-conformers/screen.py
tools/.venv-cryo/bin/python tools/cryo-conformers/screen.py --states 2
```

Results and full logs are in draft/screen-1.* and draft/screen-2.*.

Method documentation: https://computationalcryoem.github.io/ASPIRE-Python/aspire.abinitio.html
