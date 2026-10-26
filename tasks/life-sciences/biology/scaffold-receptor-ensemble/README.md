# Scaffold-mediated receptor ensemble

Original synthetic finite statistical-biophysics model of a 1024-subunit
two-conformation receptor sheet, with competing local couplings, 36 effective
long-range scaffold links and single-site biases at 26 sites. It is not a
molecular prediction for a named natural receptor or a claim to solve an open
biological research problem. All interactions and units are stated explicitly.
The full graph is nonplanar; the output is the complete conditional law of
eight reporter conformations after six distinct receptor states are clamped.

## Independently generated expectations

`tests/provenance/generate.py` defines the original seeded instance. Expected
values are computed directly from the model by `independent_source.py` and
`source_kernel.py`, without importing the reference solver or reading its
outputs. They construct their own skew matrix, use a different path root and
ordering, sequential rank-two inverse updates and chunked explicit Boltzmann
reweighting. The reference instead evaluates pivoted Pfaffian updates and uses
compiled streaming aggregation. Both implement the same exact mathematical
identity; they are not independent physical theories. Expectations are never
copied from the reference run.

The local signed planar subgraph supplies a joint marginal on the 26 sites
touched by the additional factors, clamps or reporters. Its even spin moments
follow from path-modified high-temperature graph polynomials; skew-matrix
Pfaffian ratios retain the signs of these polynomials. Odd moments vanish in
this zero-field subproblem. Character inversion yields the subproblem marginal.
Explicit scaffold and field Boltzmann factors then reweight all configurations;
normalizing and conditioning give the full, nonplanar model's requested law.
The final model has fields and does not have global spin-inversion symmetry.

Method background: Aizenman and Warzel, *Kac-Ward formula and its extension to
order-disorder correlators through a graph zeta function*,
https://arxiv.org/abs/1709.06052 ; Wimmer, *Efficient numerical computation of the
Pfaffian for dense and banded skew-symmetric matrices*,
https://arxiv.org/abs/1102.3440 . Neither paper supplies this instance or answers.

## Scientific and verifier checks

`check_small.py` sums every one of the 512 spin states on twelve small graphs
with signed local couplings, additional nonlocal interactions and fields. Both
independent implementations agree with these direct enumerations at 1e-10 for
both logarithms and every conditional probability. Full-instance comparison
and resource-limited reference results are recorded below before task freezing.

The verifier checks every one of the 256 probabilities and both logarithms,
plus normalization. Tolerances are 5e-8 absolute for each probability, 2e-6 for
each logarithm and 1e-6 for normalization. These are substantially above the
measured numerical differences, including accumulation over the larger marginal.
Three valid formatting variants are accepted; 280 malformed or scientifically
wrong submissions are rejected by author regression checks, including changes
to each joint probability that preserve normalization and an independent
product law having correct single-site marginals.

The verifier is offline and never executes agent code. It opens only bounded
regular files with symlink and nonblocking checks. Strict ASCII decimal syntax
rejects surrounding whitespace, Unicode digits, underscores and nonfinite
numbers. Contract validation precedes all scientific checks. Private expected
values and provenance are not included in the agent image.

To reproduce provenance, create a draft directory beside independent_source.py,
copy tests/data/model.json there, and run it with NumPy, SciPy and Numba installed.
No model-specific answer files are needed.

## Release policy

Claude Opus 5.5 is tried first. GPT is launched only after a normal, substantive
Claude failure; ZIP release requires normal failures of both models on identical
frozen task bytes. Quota, authentication and infrastructure errors are
inconclusive. No local result guarantees future portal acceptance.

## Completed full-instance validation

Reference: 4703.69 seconds in a one-CPU, 2 GiB offline container; all four separate-verifier tests passed.

Maximum probability difference between implementations: 1.29e-13. Log-partition difference: 4.55e-13; log-clamp difference: 1.42e-13.
