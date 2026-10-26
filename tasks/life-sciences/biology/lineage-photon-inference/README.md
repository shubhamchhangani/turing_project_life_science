# Integrated fluorescence lineage inference

This is an original finite, synthetic stochastic gene-expression problem.
Capacities, bursts, feedback, observation timing and partition rules are explicit
model assumptions; no claim is made that the instance predicts a natural cell.

## Independent expectation source

`tests/provenance/generate.py` records the fixed model parameters.
`tests/provenance/independent.py` independently assembles a row generator and
uses positive uniformization recurrences for joint photon and transcription
counts. Conditional division messages enforce complementary inheritance.
`tests/data/expectations.json` comes from this computation, not an oracle run.
The genuine solution separately assembles a column generator and uses SciPy
sparse augmented matrix exponential actions. Both implement the same declared
stochastic model, with different numerical algorithms and separate assemblies.

Eight small instances agree between methods. Analytic limits with
state-independent photon emission and disabled transcription reproduce Poisson
observation masses, zero reaction counts and the unconditional mother prior.
These are cross-algorithm and analytic checks, not exhaustive path enumeration.

## Verification and tolerances

All mother states, all joint daughter RNA pairs, all event-count probabilities
including the tail, and the evidence are checked. Absolute tolerances are 2e-8
for every probability and 2e-7 for log evidence, with normalization within 1e-6.
They comfortably exceed independent numerical disagreement and decimal rounding
error. They are numerical accuracy tolerances, not biological noise estimates.
Contract validation precedes every accuracy check. The separate offline verifier
opens bounded regular output files without following symlinks and never executes
submitted code. Reference and provenance files are not copied to the agent image.

The reference does not fit hidden planted parameters: it computes the full
likelihood and posterior for the public finite model.

## Background

Finite chemical master equations: Munsky and Khammash (2006),
https://cnls.lanl.gov/~munsky/JChemPhys_124_044104.pdf
Sparse matrix exponential action implementation:
https://github.com/scipy/scipy/blob/main/scipy/sparse/linalg/_expm_multiply.py
These describe numerical background, not this task's answers.

## Completed validation

Reference computation completed in 0.7564 seconds under one CPU and 2 GiB RAM
with no network; all four separate offline verifier tests passed. Independent
full-instance maximum probability difference was 5.56e-16 (including overflow);
log-evidence difference was 1.78e-15. Three valid-output regression cases passed
and 419 invalid-output cases were rejected, including every numerical coordinate,
incorrect independent daughter marginals and removal of the event-count tail.
