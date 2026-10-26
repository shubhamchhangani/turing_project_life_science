# Integrated photon observations across a cell division

Status: rejected. Claude Opus 5.5 low-effort solved the validated candidate in
49.361426 seconds (reward 1, two commands). GPT and ZIP creation were skipped.

The receptor family was solved repeatedly. This design instead asks for exact
Bayesian smoothing in a finite stochastic gene-expression model, with photon
counts integrated over time and anticorrelated molecular inheritance at division.
The aim is to test construction and composition of path likelihoods, not merely
increase a lattice size. Difficulty remains an empirical question.

A cell has a binary promoter and finite RNA/protein capacities. Birth, death,
burst transcription, feedback switching and photon emission are continuous-time
events. Rates and capacity boundary rules are public and explicit. These are
finite coarse-grained kinetics, not an approximation to an unspecified infinite
chemical network. A photon emission does not change the molecular state.

One mother divides into two daughters at a known time. Both inherit the mother's
promoter state. Each existing RNA and protein independently goes to daughter A
with stated probabilities; daughter B receives the complement. The two daughters
then evolve independently conditional on the partition. Observations are exact
photon counts in consecutive, finite exposures, including different stimulation
conditions. The output will include the normalized posterior at division, the
joint daughter RNA law at birth, and the posterior distribution of the total
number of transcription events across the entire observed lineage, with an
explicit overflow bin. All observations condition every reported posterior.

Author validation plan:
1. Sparse augmented-generator exponential reference.
2. Separately assembled uniformization recurrence, without copying reference
   output; dense small-state exponential checks and probability identities.
3. Full-instance agreement, resource and tolerance calibration.
4. Complete numerical coverage and malformed-output regressions, including
   incorrect independent daughter partitioning and endpoint photon shortcuts.
5. Frozen schema-1.4 package; Claude first, GPT only after a valid Claude failure;
   direct-root ZIP only after both fail and validation passes.

Sources for the underlying numerical framework (not task answers):
- Munsky & Khammash (2006), finite-state chemical master equation:
  https://cnls.lanl.gov/~munsky/JChemPhys_124_044104.pdf
- SciPy's sparse matrix-exponential action implementation:
  https://github.com/scipy/scipy/blob/main/scipy/sparse/linalg/_expm_multiply.py

No claim is made that this is a newly discovered unsolved biological problem.
