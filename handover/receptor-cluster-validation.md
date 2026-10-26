# New candidate: cooperative receptor thermodynamics

The candidate is in
`tasks/life-sciences/biology/receptor-cluster-thermodynamics/`.
It changes the scientific problem from flux fitting to exact, finite-system
statistical mechanics of a cooperative membrane-receptor lattice.

It is an original synthetic benchmark instance of an established idealized
biophysical model. It is not a claim to have discovered an unsolved biological
law. No all-models-fail claim has been established yet, and no ZIP has been
created for this candidate.

## What is validated

- 323 receptors, 646 heterogeneous interactions, periodic boundaries in both
  directions, four inhibitory interactions, six designated probe receptors.
- Full 64-state joint distribution and absolute natural log partition function.
- Genuine input-only reference; no expected probabilities in the agent image.
- Independent verifier uses a different graph ordering, a different Pfaffian
  reduction, different disorder paths, and a different expansion for log Z.
- Direct enumeration on twelve small-system/coupling cases checks both
  implementations. Maximum error: 1.60e-14.
- Full-size maximum difference: 4.97e-13 in joint probabilities and 5.69e-14
  in log Z. The required absolute accuracy is 1e-7.
- Packaged reference ran offline with one CPU and 2 GiB. Its outputs passed all
  six tests in a separate offline verifier in 215.44 seconds, reward 1.
- 32/32 static package checks passed.
- Four valid artifact variants were accepted; 29 invalid variants were rejected.
  A perturbation that preserves all correlations of orders below six is rejected,
  so partial-summary agreement cannot substitute for the full requested joint
  distribution.

## Difficulty trials

The first fresh GPT-6-sol max run stopped at a Codex usage limit after inspecting
the data and beginning an implementation. Its reward 0 is **inconclusive** and
does not count as model failure.

After checking that the usage window had reset and capacity was available, a
fresh run was started under
`jobs/receptor-torus-frontier/gpt6-sol-max-receptor-torus-npx-retry1/`.
It receives exactly the frozen instructions and input files, the full
18000-second allowance, and no reference code or previous-trial hints.
The retry completed normally: reward 1, all six verifier tests passed. Agent
execution took 475.75 seconds. GPT-6-sol recovered the full joint distribution
and absolute log partition function. This candidate fails the requested
difficulty condition and must not be packaged for submission.

The user has selected local Claude CLI login for the required Opus 5.5
low-effort trial. Login is not yet available. The prepared runner uses a Claude
subscription and has no API-key billing fallback. Credentials are not printed
or included in the task package.

## Evidence

- `tools/receptor-torus/draft/physics-validation.json`
- `tools/receptor-torus/draft/regressions.json`
- `tools/receptor-torus/draft/package-before-trial.json`
- `tools/receptor-torus/draft/offline-verifier/ctrf.json`
- `tools/receptor-torus/draft/validation-status.json`
- `tools/receptor-torus/draft/release-decision.json`

The prior successful MFA trials remain rejected candidates. Package quality
checks and a working expert reference do not establish frontier-model failure.
Even completed local trials would not guarantee the result of a different
Turing agent harness or a future run.
