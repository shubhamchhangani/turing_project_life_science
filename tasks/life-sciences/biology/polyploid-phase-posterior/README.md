# Exact polyploid phase uncertainty

Original synthetic computational-genomics benchmark. Its finite Bayesian
model averages over chromosome phases and read origins. It does not ask for
an optimal phase assembly, nor claim to solve an open biological problem.
The prescribed prior is uniform over dosage-compatible labelled allele
columns, independent across variants. Read origins are independent and
uniform over six copies; sequencing errors are symmetric substitutions.

Scientific motivation: phase ambiguity in polyploid genomes requires
propagating uncertainty in read assignment. See Hosseini et al.,
*pHapCompass: Probabilistic Assembly and Uncertainty Quantification of
Polyploid Haplotype Phase*, https://arxiv.org/abs/2512.04393.
The particular finite prior, dataset and requested observables here are
original; this benchmark does not claim to reproduce that paper's method.

## Original data

Author generator: `tools/phase-posterior/generate.py`, seed 20261102,
64 variants and read span 8. Each of the 64 genotype columns independently
chooses three alternate copies from six uniformly. Two independently
sampled reads end at each index 7--63. Each read chooses one homolog;
its per-call error probabilities are sampled from 0.025, 0.08 and 0.19,
and errors are then sampled from those probabilities. Only the observed
calls, exact dosages and model parameters are public. No generating phase
or read-origin labels are supplied to the reference.

## Independent scientific expectations

The private expectation fixture is computed by the author program
`tests/provenance/partition_source.cpp`, directly from the public
probability model. That program separately evaluates modified partition
functions. Replacing a read's emission by its zero-mismatch contribution
and dividing by the ordinary partition function gives its error-free
probability. Replacing it by the mismatch-count-weighted emission gives
its expected number of errors. These identities follow from linearity in
each read's factor. No fitted or submitted solution is used as ground truth. Its independently
written forward-only implementation and rebuild utility are included under
`tests/provenance/` for review; they are not run by the verifier.

The packaged reference instead uses forward and backward messages and
local posterior expectations. It does not read the private fixture or
import the author expectation program. Both calculations sum complete
finite-state probability distributions, including the normalized priors.

Homolog relabellings can be grouped into equivalence classes, provided
their multiplicities are retained. The author independently checked the
state count using Burnside's lemma and the normalized prior's stationarity.
For seven remembered columns, 1,280,000,000 labelled histories correspond
to 1,822,528 equivalence classes. Exhaustive enumeration of all 160,000
labelled phases in six small instances agrees with the reference posterior
calculation to within 1.2e-12 and with log evidence to within 1.8e-15.

## Verification

The separate offline verifier checks the literal output contract before
all scientific values. It compares the normalized log evidence and both
posterior quantities for every one of the 114 reads to the independently
computed fixture. It never executes agent code or imports solution code.
It uses immutable verifier data and fixed tolerances. Decimal syntax is
checked before numeric conversion; whitespace, Unicode numerals, infinities
and NaNs are rejected. Files are opened nonblocking with no symlink
following, checked to be regular, and bounded to 1 MiB each.

The absolute tolerances are 2e-6 for log evidence and 5e-7 for each posterior
quantity. They accommodate decimal rounding and finite-precision numerical
accumulation, not Monte Carlo sampling error. A nine-decimal output has
rounding error at most 5e-10. The full-size independent partition-ratio calculation differs from the
packaged posterior reference by less than 9e-13 for every quantity. The
reference inference ran offline with one CPU and 2 GiB in 115.60 seconds
(excluding compilation). The separate offline verifier passed all three
tests in 0.01 seconds. Six valid artifact representations and 57 invalid
cases behaved as specified, including each statistic on the final read.

## Status

Draft. No model-failure claim or submission ZIP exists. The reference,
independent expectations and artifact regressions are validated. Fresh
frontier trials remain required before considering this candidate for release.
