# Inverse receptor cooperativity: validation record

Candidate: `tasks/life-sciences/biology/inverse-receptor-cooperativity/`.
This original synthetic life-sciences benchmark estimates 646 interactions
and three local ligand biases from equilibrium observations. It is an
identifiable finite-system inverse problem. It is not a claim that an open
biological problem has been resolved.

## Completed checks

- Schema 1.4, exact allowed root contents, no hidden files or canary markers:
  33/33 static checks passed.
- Both physics implementations independently match exhaustive enumeration on
  nine small-system cases each, including inhibitory interactions and signed
  fields. Errors are below 2e-14.
- Exact analytic Fisher information matches exhaustive covariance to 2.45e-15.
- Full-size independently computed moments agree to 1.06e-13, and log Z to
  2.85e-13.
- Input-only reference fits every parameter within 1e-9 of its independently
  drawn generating value. It runs offline in 115.62 seconds on one CPU and
  2 GiB; no injected parameters or verifier code are accessible to it.
- Separate offline verifier: 4/4 tests pass, reward 1, 24.73 seconds.
- Artifact regressions: four valid variants accepted, 33 invalid variants
  rejected. Includes numeric whitespace, malformed coordinates, duplicates,
  missing data, nonfinite numbers, symlinks and nonblocking FIFO handling.
- A simultaneous 4e-5 coupling perturbation remains inside each parameter's
  tolerance but is rejected for failing the complete equilibrium-moment fit.
  Thus every parameter and the collective scientific fit are covered.

## Difficulty evaluation

A fresh isolated GPT-6-sol max-effort run via `npx @openai/codex@0.164.0`
first ran under
`jobs/receptor-torus-frontier/gpt6-sol-max-inverse-receptor-npx/`.
It sees only instructions and public inputs, with the full 18000-second
allowance. The current subscription window had 21% capacity remaining at
launch. That run was interrupted by the usage limit and is inconclusive. A fresh
retry is running under `gpt6-sol-max-inverse-receptor-npx-retry1`, launched
after the quota reset with 97% primary-window capacity available.

Claude Pro login is verified using macOS Keychain. An isolated Opus 5.5
low-effort trial is running under `opus-5-5-low-inverse-receptor`. The initial
launch stopped before creating a trial because its access token needed
refreshing; the current token was subsequently verified as valid. No
additional API-billed trial is authorized or launched.
The fresh GPT retry completed normally and solved the task: reward 1, all
four tests passed, 1358.61 seconds of agent execution. The Opus run was then
cancelled to conserve subscription quota; its cancellation does not count
as a model failure. This candidate is rejected for the requested difficulty
condition. No ZIP has been created.

The preceding forward receptor candidate was solved by GPT-6-sol in 475.75
seconds of agent execution. It is rejected for the requested difficulty
condition, like the earlier metabolic-flux candidates.

## Evidence

Author reports are in `tools/receptor-torus/draft/field-inverse/`:
`generation.json`, `newton-report.json`, `independent-validation.json`,
`regressions.json`, `package-before-trial.json`, `offline-verifier/` and
`validation-status.json`. Author-only reports are excluded from task roots
and any eventual directly rooted ZIP.

Local package checks cannot guarantee Turing portal acceptance. Even valid
failures from both local target models would not guarantee identical results
under another harness or on a later run.
