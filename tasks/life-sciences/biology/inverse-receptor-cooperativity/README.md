# Inverse receptor cooperativity

This original synthetic computational-biophysics instance asks for 646
heterogeneous pair interactions and three local ligand biases in a periodic
17 by 19 membrane patch. It concerns an idealized two-state receptor model,
not experimental data or a claim to have solved an open biological question.

## Scientific basis and identifiability

Cooperative receptor lattices model collective membrane sensitivity; see
[Duke and Bray (1999)](https://pubmed.ncbi.nlm.nih.gov/10468569/).
Here the supplied moments are expectations under the explicitly stated finite
Gibbs distribution. Positive and negative interactions and nonzero local
biases are included. Periodicity and all signs are defined in the instructions.

The sufficient statistics are the 646 distinct edge spin products and three
distinct single-site spins. They are distinct nonconstant Walsh characters on
the full configuration space. No nonzero linear combination is constant.
Every configuration has positive probability at finite parameter values;
therefore their covariance matrix, the Hessian of log Z, is positive definite.
The parameter-to-mean map is injective. This gives an identifiable inverse
problem rather than requiring an arbitrary representative of a degenerate fit.

## Original generation and independent targets

The author generator is `tools/receptor-torus/field_inverse_prototype.py` in
the development workspace, seed 20261028. For each edge it draws a coupling
from a spatial sinusoid plus bounded variation and replaces a deterministic
sparse subset with inhibitory couplings. Couplings are rounded to nine decimal
places before generating observations. Three prescribed local biases are
used. Expectations are rounded to twelve decimal places. There is no sampling
noise. The private `tests/data/injected_parameters.json` contains the original
parameter draw, not fitted reference output.

The agent image contains only `model.json`, `edge_moments.csv` and
`site_moments.csv`. Neither the original parameters nor reference algorithms
are copied into that image. The verifier has immutable copies of its own data.

## Genuine reference and independent physics

The reference reads public inputs only. It expands the three local-field
factors around a zero-field toroidal partition function. Four even subsets
are needed. Signed Pfaffians with disorder lines give the terms; analytic
derivatives yield expectations and the Fisher information matrix. A damped
Newton method solves the strictly convex moment-matching problem.

The loop formulation follows
[Cimasoni (2010)](https://arxiv.org/abs/1004.3158) and
[Chelkak, Cimasoni and Kassel (2015)](https://arxiv.org/abs/1507.08242).
The verifier implements the decorated graph with independent vertex indexing,
a separate edge map, different disorder paths, and Householder Pfaffians.
It recomputes all 649 expectations from the submitted parameters. It does not
execute solver code or read reference artifacts.

Both physics implementations were checked separately against enumeration on
3 by 3, 3 by 4 and 4 by 4 tori in three interaction regimes, including
inhibitory couplings and signed fields. Maximum discrepancy was below 2e-14.
The analytic Fisher matrix agrees with exhaustive covariance to 2.45e-15.
At full size, independent moments agree to 1.06e-13 and log Z to 2.85e-13.
The input-only Newton fit recovered all parameters to below 1e-9 in its
author validation run, with maximum moment residual below 3.64e-10.

## Coverage and tolerance rationale

Every parameter and every supplied observation is checked. The 5e-5 parameter
tolerance is over 50,000 times the measured numerical recovery error. The
2e-6 moment tolerance allows substantial margin over twelve-decimal input
rounding and independent implementation error while preventing a collection
of individually small parameter errors from producing a poor equilibrium
fit. The smallest full-size Fisher eigenvalue at the recovered solution is
about 0.07568, consistent with the stated numerical identifiability.

The output contract is checked before physics: bounded regular files,
literal headers, complete unique keys, strict ASCII decimal fields, finite
values and explicit bounds. Every test has a docstring. The verifier is
separate, offline and bounded by a runner timeout. Accuracy is evaluated only
from submitted outputs, with no restrictions on the solver's method.

## Release status

Schema 1.4 and the agreed task-root structure are used. Behavioral regressions,
offline execution evidence and fresh-model trial results are recorded outside
the task package. Difficulty remains an empirical question. A candidate solved
by a target model must not be released as meeting the requested all-models-fail
condition. No such result is asserted here.

A fresh GPT-6-sol max-effort trial solved this candidate in 1358.61
seconds of agent execution: reward 1, all four verifier tests passed. It must
not be packaged as satisfying the requested model-failure condition. The
parallel Opus trial was cancelled after this result, not counted as a failure.
