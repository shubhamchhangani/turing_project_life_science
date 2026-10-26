# Exact phase posterior: rejected after Claude solution

Candidate: `tasks/life-sciences/biology/polyploid-phase-posterior/`.
No ZIP has been created. Current state is written to
`handover/phase-posterior-progress.json`.

## Rejections confirmed this session

- Spatial ligand cooperativity: Claude Opus 5.5 low solved it, reward 1,
  four tests passed, 933.91 seconds of agent execution.
- Polyploid haplotype assembly: GPT-6-sol max solved it, reward 1,
  three tests passed, 236.89 seconds; correction cost 19878 <=19885.
  Its unfinished Claude trial was cancelled and is not a model failure.

Neither solved candidate was packaged.

## New scientific task

Compute normalized Bayesian evidence and both read-error posterior
quantities for all 114 reads of a synthetic hexaploid chromosome. The task
averages over phase and read-origin uncertainty; it does not request a
single best-fit assembly. The finite prior and every observation likelihood
are explicitly defined. This is an original benchmark, not a claim to have
discovered an unsolved biological problem.

## Completed validation

- 34/34 package checks: schema 1.4, correct namespace and allowed root
  structure, no canary markers or hidden/cache files, shell errexit,
  documented tests and explicit execution/output contracts.
- Six independent exhaustive labelled-state checks, each enumerating
  160000 phases: posterior agreement within 1.2e-12.
- Burnside's lemma independently confirms all 1822528 seven-column orbits;
  their multiplicities sum to 1280000000 labelled histories.
- Full-size scientific expectations were computed by a separate
  forward-only modified-partition evaluator, not by copying the reference
  solution's posterior outputs. Source and reproduction utility are included
  under `tests/provenance/` for review.
- All 228 posterior quantities agree with the packaged forward/backward
  reference within 8.95e-13; log evidence agrees within 5.69e-13.
- Input-only reference ran offline with one CPU and 2 GiB. The inference
  engine took 115.60 seconds, excluding compilation. Temporary checkpoints
  stay below the 10 GiB storage allowance.
- Separate offline verifier: 3/3 tests, reward 1, 0.01 seconds.
- Six valid representations accepted; 57 invalid cases rejected. Coverage
  includes both statistics for the final read, missing prior normalization,
  independent-read approximations, numeric whitespace/Unicode/nonfinite
  values, bounds, IDs, CSV shape, size, symlinks, FIFOs and directories.
- Twelve conditional release checks passed.

Frozen fingerprint:
`5253e4c18a71fb19664c5a5e12cfbd00bdf00ed59379d7d0a9d2c94a333bef4f`.
Author evidence: `tools/phase-posterior/draft/`.

## Authorized evaluation workflow

Claude Opus 5.5 low launched on 24 October 2026 at 15:10 IST, using the
existing Claude Pro login. GPT primary subscription usage was 89%, so the
GPT-6-sol max npx trial is queued no earlier than 19:22:35 IST that day,
after its reported reset plus one minute. Capacity is checked again before
launch. Both trials retain the full 18000-second task budget. There is no
paid API-key fallback.

Worker: `tools/phase-posterior/run_trials.py`, PID recorded in the progress
JSON. Log: `tools/phase-posterior/draft/trial-worker.log`.
Jobs are under `jobs/receptor-torus-frontier/`, named
`opus-5-5-low-phase-posterior` and `gpt6-sol-max-phase-posterior-npx`.

One successful model solution rejects the candidate. Quota, setup, runtime
and detected machine-sleep interruptions are not qualifying failures. Only
two normal completed failures on the unchanged, validated package permit
a ZIP in Downloads. Its files will be directly at the archive root; file
membership and bytes are audited against the frozen manifest. Existing
ZIPs are not overwritten.

The computer must remain running for active local trials to progress. The
workflow does not guarantee failure in a later Turing-platform run, and no
all-models-fail claim has been made.

## Completed result

Claude Opus 5.5 low solved this candidate: reward 1, all three verifier
tests passed, 84.37186 seconds of agent execution, three completed commands,
no exception and no quota interruption. The worker rejected the candidate
on 24 October 2026 at 09:47:16 UTC (15:17:16 IST).

The trajectory contains an input-only C++ implementation of orbit reduction
and forward/backward marginalization. It correctly preserves labelled-state
multiplicities and computes every required posterior. The queued GPT trial
was cancelled before launch. No ZIP was created and no trial remains queued
for this candidate. Earlier launch/queue details above are historical.
