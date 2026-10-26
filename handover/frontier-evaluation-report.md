# Local model trials: both candidates solved; no new ZIP

The requested fresh npx trials used GPT-6-sol at max reasoning effort through
the existing ChatGPT Codex login. Each ran in an isolated Harbor Docker
environment containing only the corresponding instructions and input data.
Neither trial received the reference solution, verifier, injected truth, or
previous trial conversation. Each had the declared 18000-second allowance.

| Candidate | Agent execution | Verifier result | Release decision |
|---|---:|---|---|
| Metabolic flux analysis v2: transient pulse-chase, one population | 11m 18s | 16/16 tests, reward 1 | Do not submit |
| New population deconvolution task: four tracer programs, two populations, 44 unknowns | 10m 14s | 19/19 tests, reward 1 | Do not submit |

The v2 agent reported weighted SSR 1961.0084 against limit 2110. The new
two-population agent reported 3896.9618 against limit 4090. These are real
completed solution attempts, not API quota failures. The earlier API-key run
that returned `ApiUsageLimitError` remains classified as inconclusive.

The user explicitly requested no ZIP for a solved candidate. No new ZIP was
created during these trials or for the new task. Existing older candidate ZIPs
were preserved, but should not be submitted as evidence of model difficulty.

## What the experiments established

The earlier portal errors were problems with packaging and verification.
Those can be tested directly. The user-provided later screenshot reports a
different failure: the target models solved the task. Improving package quality
does not itself make a task hard enough for that requirement.

The new candidate added a real scientific inference layer: population-specific
fluxes and pools inferred from molecular mixtures. It was scientifically
checked and has a genuine input-only reference solver. However, this GPT trial
implemented the coupled simulator and successfully optimized it. Adding this
particular complexity did not produce the required model failure.

Neither task qualifies as an all-models-fail submission. No claim is made that
an untested open research question would have a reliable correct-answer
verifier, or that any untested task is guaranteed to defeat frontier models.
The requested difficulty goal remains unmet.

## Validation retained for reuse

- `tools/preflight_task.py`: root layout, schema 1.4 fields, task namespace,
  metadata, artifacts, resource/network contracts, canary and hidden/cache
  files, Python syntax/test docstrings, and runner errexit/shell syntax.
  It can also inspect an existing ZIP for direct-root layout and source-byte
  equality. It never creates a ZIP. New candidate: 40/40 checks passed.
- New task artifact regressions: four valid cases accepted, 28 invalid cases
  rejected, including PGI exchange inflation, compensating population flux
  errors, whitespace/Unicode numeric fields, duplicate rows, FIFO files, and
  size limits.
- New packaged reference: two independent starts, one CPU, 2 GiB, no network;
  373.69 seconds, best SSR 3896.9940, all 19 verifier tests passed.
- Four additional noise realizations passed the target tolerances. These were
  local sensitivity fits initialized from an independently fitted estimate.
- Independently implemented pooled-ion predictions agreed within 5.11e-8.
- `release_decision.py` distinguishes solved tasks, substantive failures, and
  infrastructure/no-attempt errors. A successful trial blocks packaging.

These checks prevent the known tested regressions; they cannot promise that a
future platform review will never identify another issue.

## Evidence locations

- V2: `jobs/mfa-v2-frontier/gpt6-sol-max-v2-npx-chatgpt/evaluation-status.json`.
- New task: `jobs/mfa-mixture-frontier/gpt6-sol-max-mixture-npx-chatgpt/evaluation-status.json`.
- New source: `tasks/life-sciences/biology/metabolic-population-deconvolution/`.
- Validation and release decision: `tools/metabolic-mixture/draft/validation-status.json`
  and `tools/metabolic-mixture/draft/release-decision.json`.

Raw trial logs, output artifacts, and rewards are preserved. README status
updates after the trials only document the results; instructions, inputs,
solution code, and verifier rules were not changed during either trial.
