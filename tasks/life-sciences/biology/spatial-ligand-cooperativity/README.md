# Spatial ligand cooperativity

This original synthetic biophysics benchmark reconstructs 646 heterogeneous
pair interactions and twelve spatially separated local biases from equilibrium
observations on a periodic 17 by 19 receptor lattice. The finite-system model
and the complete output contract are stated in `instruction.md`.

## Scientific basis and identifiability

Cooperative two-state receptor lattices are idealized models of membrane
sensitivity; see [Duke and Bray (1999)](https://pubmed.ncbi.nlm.nih.gov/10468569/).
Positive and negative local biases respectively favour active and inactive
states. The synthetic observations are not measurements from a real membrane,
and this task is not a claim to have discovered an unsolved biological law.

The 658 sufficient statistics are distinct nonconstant Walsh characters:
646 edge spin products and twelve single-site spins. No nonzero linear
combination is constant on the full configuration space. Every configuration
has positive probability for finite parameters. Consequently their covariance
matrix, the Hessian of log Z, is positive definite. The log-partition gradient
is injective, giving a unique parameter vector for exact moments.

## Original data generation

The development generator is `tools/receptor-torus/benchmark_sparse_fields.py`,
using NumPy's default generator with seed 20261030. It combines a spatial
sinusoid and bounded coupling variation with a deterministic sparse rule for
inhibitory edges. Twelve signed local biases are prescribed. Couplings are
rounded to nine decimal places before expectation calculation; observations
are rounded to twelve decimal places, with no added sampling noise.

The private injected parameters are the original draw, not fitted reference
outputs. The agent image contains only the three public input files. The
reference and verifier implementations and original parameters are absent
from that image. The verifier uses its own immutable copies of all data.

## Genuine reference and independent verification

`solution/solve.py` reads only public data and minimizes the strictly convex
moment-matching objective. The partition function expands into 2048 even
field subsets and four toroidal sectors. Signed Pfaffians and low-rank matrix
updates permit reuse across subsets. Analytic derivatives provide all 658
model expectations for optimization. No generating parameter is used to
initialize or constrain the fit beyond the published bounds.

The toroidal loop formulation follows
[Cimasoni (2010)](https://arxiv.org/abs/1004.3158) and
[Chelkak, Cimasoni and Kassel (2015)](https://arxiv.org/abs/1507.08242).
The independent verifier uses a separate vertex-ordered decorated graph,
different edge indexing and disorder paths, Householder Pfaffians, and linear
solves for the updates. It recomputes every supplied expectation from the
submitted parameter vector. It neither executes the submitted solver nor
compares its implementation with the reference.

Both implementations were checked against exhaustive enumeration on small
tori with signed couplings and three, five and eight fields. Errors were
below 7e-14. Full-size calculations with twelve fields agreed to 2.85e-13.
Three additional twelve-field cases on a 4 by 4 torus were exhaustively
enumerated in weak, stronger and mixed-sign coupling regimes. Maximum error
was 1.32e-13. The packaged input-only reference ran offline on one CPU with
2 GiB in 745.39 seconds. Maximum recovery errors were 9.20e-7 for couplings
and 7.97e-8 for biases. Its outputs passed all four tests in a separate
offline verifier in 16.31 seconds. Four valid artifact variants were accepted
and 34 invalid variants rejected, including an error in the twelfth bias and
a collective fit error despite individually tolerated parameter changes.

## Coverage, thresholds and execution

Every coupling and every bias must be within 5e-5 of its original generating
value. In addition, all 658 predicted equilibrium observations must be within
2e-6 of the supplied values. Both thresholds are stated in the instruction.
The moment tolerance allows substantial margin over twelve-decimal rounding
and independent numerical error. Checking the collective fit prevents many
small, individually tolerated parameter errors from producing an inaccurate
equilibrium distribution. The parameter threshold is over 54 times the largest
measured reference recovery error. It concerns a deterministic, identifiable
synthetic problem, not uncertainty in experimental measurements.

Bounded regular-file, literal-header, unique-key, ASCII-number, finiteness and
bound checks precede physics. All parameters and all observations are covered.
Every test is documented. The verifier runs separately without networking or
runtime installation, with a fixed timeout and binary result.

Schema 1.4 and the agreed root layout are used. Author reports stay outside
the task root. Any eventual ZIP must expose the allowed files and directories
directly, without an enclosing directory. This candidate was solved by Claude Opus 5.5 (low effort): reward 1,
all four tests passed, and 933.91 seconds of agent execution. It is rejected
for the requested difficulty criterion and must not be packaged for submission.
The GPT trial was cancelled before launch. See the author handover report
for the completed trial evidence.
