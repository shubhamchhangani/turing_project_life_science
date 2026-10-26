# Frustrated receptor ensemble

Original synthetic finite model of competing conformational interactions in a
receptor sheet. Positive and negative pair couplings represent preferences for
equal and opposite conformations. This is a coarse-grained statistical-biophysics
benchmark, not a prediction about any named natural receptor or an unsolved
biological discovery. The finite model, all signs, clamps and normalization
conventions are explicitly stated in instruction.md.

## Data and independently generated expectations

`tests/provenance/generate.py` defines the original seeded graph and signed
couplings. It embeds no solutions. The expectation source is
`tests/provenance/independent_source.py`: it constructs a complex skew matrix
from directed edge geometry and evaluates signed moment ratios using sequential
rank-two inverse updates. It uses a different root and path ordering from the
reference, long-double update arithmetic and explicit character sums. It does
not import the reference or read reference output.

The input-only reference in solution/ constructs a skew Kac-Ward matrix, solves
for a restricted inverse and evaluates small Pfaffian updates by pivoted skew
elimination. A Walsh transform recovers the joint distribution. Both methods
retain correlation signs; taking the positive square root of a determinant is
insufficient with frustrated interactions. The two implementations share the
same exact graph identity and are not independent physical theories.

Mathematical background: Aizenman and Warzel, *Kac-Ward formula and its extension
to order-disorder correlators through a graph zeta function*,
https://arxiv.org/abs/1709.06052 ; Wimmer, *Efficient numerical computation of the
Pfaffian for dense and banded skew-symmetric matrices*,
https://arxiv.org/abs/1102.3440 . Neither supplies this instance or its answers.

For an even set of spins, a path set with precisely these odd-degree vertices
converts the spin moment into a modified even-subgraph polynomial. Reciprocal
edge-weight substitutions become skew rank-two updates. A Pfaffian ratio gives
the signed polynomial ratio; the accompanying edge-product factor gives the
moment. Odd moments vanish by global spin inversion. Character inversion and
conditioning then produce every requested probability.

## Verification and numerical evidence

`tests/provenance/check_small.py` independently enumerates all 512 spin states
on each of twelve signed 3-by-3 grids and compares both implementations, including
both logarithms and every conditional probability. All cases pass at 1e-10.
The full-instance cross-implementation difference is below 2.6e-16 for every
probability and 5.7e-13 for either logarithm. Expectations in tests/data are
produced by the independent source, never copied from the reference output.
The source can be rerun by creating a draft directory beside it, copying
`tests/data/model.json` there and running independent_source.py.

All 256 probabilities and both logarithms are checked. Probability tolerance
5e-8 and logarithm tolerance 2e-6 allow substantially more floating-point and
serialization error than observed. Normalization tolerance is 1e-6. Correct
single-site marginals alone do not satisfy the complete joint-law contract.
Adversarial checks accepted three valid format variants and rejected 280
incorrect outputs, including individual probability perturbations preserving
normalization and a product distribution with correct individual marginals.

The verifier runs offline, opens only bounded regular files and never executes
agent code. Its strict ASCII decimal grammar rejects whitespace, Unicode digits,
underscores and nonfinite values. Contract validation precedes scientific checks.
Private expected values and provenance are not copied into the agent image.

## Release status

The reference completed offline with one CPU and 2 GiB RAM and passed all four
tests in the separate offline verifier. Candidate awaiting ordered frontier trials. Claude must fail a valid attempt before GPT is launched. A successful
solve, quota interruption or infrastructure failure cannot justify ZIP release.
Local trial outcomes do not guarantee future portal acceptance.
