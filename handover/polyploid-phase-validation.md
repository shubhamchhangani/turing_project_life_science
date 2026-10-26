# Polyploid haplotype assembly: rejected after GPT solution

Candidate: `tasks/life-sciences/biology/polyploid-haplotype-assembly/`.
No ZIP exists. This is an original synthetic computational-genomics task,
not an unsolved biological discovery or a claim that all models will fail.

The former spatial-ligand candidate is rejected: Claude Opus 5.5 low solved
it with reward 1, all four tests passing in 933.91 seconds of agent execution.
The new task uses discrete chromosome phasing rather than receptor fitting.

## Completed checks

- 32/32 package checks: schema 1.4, allowed root entries, no canary markers,
  no hidden/cache files, documented tests, explicit execution contract.
- Independent target: original synthetic chromosomes have weighted MEC
  19885. The task asks for any dosage-correct assembly at or below that cost.
  No unique hidden chromosome sequence or unproved global optimum is required.
- Input-only reference: cost 19873, 175.52 seconds, offline Docker with one
  CPU and 2 GiB memory; no original chromosomes or answer file supplied.
- Separate immutable-data offline verifier: 3/3 tests, reward 1, 0.10 seconds.
- Twelve exhaustive small-instance scoring/bookkeeping checks.
- Five valid artifact representations accepted and 26 invalid cases rejected.
  These cover whitespace, Unicode, missing/duplicate IDs, all variant dosages,
  poor phasing with correct dosages, objective equality at the boundary,
  symlinks, FIFOs, directories, size limits and malformed CSV.
- Twelve checks of the conditional release decision.

Evidence: `tools/polyploid-phase/draft/`. Package fingerprint:
`a7ce73b19acb76b8310eb092787f62b54479113651cf95588904f1c823397a37`.

## Frontier evaluation

Fresh isolated Opus 5.5 low and GPT-6-sol max trials launched on
24 October 2026 at approximately 14:39 IST. Existing subscription
credentials are used; no API-key billing fallback is enabled. Both have
the full 18000-second task budget and the same frozen package.

The first launcher attempt encountered a broken local Claude executable
symlink before any model run. Existing login was confirmed with the pinned
npx CLI and the launcher was repaired. That setup error is not a model
failure. Its records are retained separately.

Current state: `handover/polyploid-phase-progress.json`.
Worker: `tools/polyploid-phase/run_trials.py`.
Worker log: `tools/polyploid-phase/draft/trial-worker.log`.
Jobs: `jobs/receptor-torus-frontier/opus-5-5-low-polyploid-phase/` and
`jobs/receptor-torus-frontier/gpt6-sol-max-polyploid-phase-npx/`.

One solved result rejects the candidate and cancels unneeded work. Quota
interruptions, missing models, changed packages, abnormal execution, and
machine-sleep interruptions do not qualify for release. Only two normal
completed failures on this same validated package permit a direct-root ZIP
in Downloads. The archive's paths and bytes are checked against the frozen
source. Local trials cannot guarantee acceptance or failure in a later
Turing-platform run. The computer must remain running for trials to progress.

## Completed result

GPT-6-sol max solved the candidate: reward 1, all three tests passed,
236.891953 seconds of agent execution, 12 completed commands, no exception
or quota interruption. Submitted weighted MEC was 19878 <=19885.
The worker recorded rejection at 2026-10-24 09:17:11 UTC and stopped the
unneeded Claude run. The cancelled run is not a model failure. No ZIP
was created. The earlier running-state description is retained as history.
