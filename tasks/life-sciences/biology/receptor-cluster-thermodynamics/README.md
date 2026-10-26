# Cooperative receptor thermodynamics — candidate under validation

This original synthetic instance concerns the equilibrium joint activity of
six receptors embedded in a heterogeneous 17 by 19 periodic membrane patch.
There are 323 binary receptors and 646 undirected interactions. Both attractive
and inhibitory couplings are present. The task requests an absolute partition
function and all 64 joint probe probabilities, including information beyond
pairwise correlations.

This is a computational biophysics benchmark, not a claim to have discovered
an unsolved biological question or to predict a particular real membrane.
Exact algorithms exist for this zero-field model. Whether the instance is hard
for a particular agent must be established by completed trials.

## Scientific basis and original data

Cooperative two-state receptor lattices are established idealized models of
membrane sensitivity: [Duke and Bray, *Heightened sensitivity of a lattice of
membrane receptors* (1999)](https://pubmed.ncbi.nlm.nih.gov/10468569/).
The present heterogeneous coupling instance and probe placement are original;
they are not experimental measurements or a published benchmark answer.

The author generator is `tools/receptor-torus/build_candidate.py` in the
development workspace. It uses NumPy's default generator with seed 20261026,
a spatial sinusoid plus bounded coupling variation, and a deterministic sparse
inhibitory-edge rule. It writes the decimal model parameters directly. No
expected probability, fitted parameter, oracle result, or solver implementation
is included in the agent image.

## Genuine reference calculation

`solution/solve.py` reads only the two public input files. It constructs four
signed loop polynomials on the dual torus and obtains even spin moments through
disorder-line insertions. The 64 probabilities follow by a finite Walsh
transform. Signed, pivoted skew elimination avoids silently choosing the
positive square root of each determinant. The reference computes 32 even
moments; global inversion gives all odd moments exactly zero.

The toroidal loop/Pfaffian formulation is based on
[Cimasoni, *A generalized Kac–Ward formula* (2010)](https://arxiv.org/abs/1004.3158)
and [Chelkak, Cimasoni and Kassel, *Revisiting the combinatorics of the 2D Ising
model* (2015)](https://arxiv.org/abs/1507.08242).

## Independent verification and thresholds

The verifier does not load reference output or stored expected probabilities.
It recomputes all probabilities from its own immutable copy of the model.
Its decorated graph uses vertex ordering, a separately implemented edge map,
and PFAPACK Householder reduction. Disorder paths have different pairings and
routings from the reference. The absolute partition function is independently
obtained from the primal high-temperature expansion, while the reference uses
the dual low-temperature expansion.

Both implementations were separately checked against direct enumeration of
all receptor assignments on 3 by 3, 3 by 4, and 4 by 4 lattices, in four coupling
regimes including negative interactions. The largest observed discrepancy in
these checks was below 2e-14. Full-size cross-checks, offline execution, artifact
regressions, and model-trial evidence are recorded outside the task package.

The 1e-7 absolute accuracy requirement concerns deterministic thermodynamics
of exactly specified decimal inputs, not recovery of noisy experimental truth.
It leaves a substantial margin over observed numerical errors and decimal
rounding. Probabilities, normalization, every marginal, every complementary
pattern, and the full joint distribution are checked. Merely matching
single-site or pairwise summaries is insufficient.

Contract checks run before physics. Input artifacts are opened nonblockingly,
checked for regular-file type and bounded size, and parsed with an ASCII numeric
grammar before float conversion. The separate verifier installs nothing and
uses no network at runtime. All six tests are documented.

## Packaging and release status

Schema 1.4, `turing/` namespace, strict root layout, and a test runner beginning
with `set -euo pipefail` are used. A release ZIP, if produced, must expose the
task files and directories directly at its root. Author reports and model logs
belong outside the task package.

A fresh GPT-6-sol max-effort trial solved this task: reward 1 and all six
verifier tests passed. Agent execution took 475.75 seconds. This candidate
does not meet the requested difficulty condition and must not be packaged
as an all-models-fail submission. The first trial was interrupted by a usage
limit and is inconclusive; it is not counted as a model failure.
