# Lineage photon inference — rejected after Claude solve

Candidate: tasks/life-sciences/biology/lineage-photon-inference.
A different scientific task family replaces the repeatedly solved receptor
lattice. It asks for exact smoothing from integrated photon counts over a
mother and two daughter cells, including complementary molecular inheritance
and the posterior number of transcription events across the full lineage.

Completed validation:
- Sparse augmented matrix exponential reference and independently assembled
  positive uniformization source agree within 5.56e-16 for all probabilities
  and 1.78e-15 for log evidence.
- Eight small cross-algorithm cases plus analytic state-independent-photon and
  zero-transcription limits pass. These are not exhaustive path enumeration.
- One CPU / 2 GiB offline reference: 0.7564 seconds.
- Four separate offline verifier tests pass.
- Three valid formats accepted; 419 malformed or scientifically incorrect
  outputs rejected, including mutations to every output coordinate.
- 36/36 static package checks pass. Portal acceptance is not guaranteed.

Claude Opus 5.5 low-effort earned reward 1 in 49.361426 seconds, with two
completed commands and no quota or infrastructure error. All four tests passed.
The solver built an augmented sparse generator, integrated photon/count
likelihoods, retained an absorbing overflow bin, and contracted complementary
daughter inheritance exactly. This is a genuine solve. GPT was not launched.
Strict worker state:
handover/lineage-photon-progress.json. A genuine Claude failure is required
before GPT-6-sol max can start. Both must fail on the same frozen package,
without quota/setup errors, before a direct-root ZIP may be released.

This candidate fails the requested difficulty criterion. No ZIP was created.
The frozen package is retained solely as evidence of a rejected candidate.
Frozen fingerprint: 50d5815bdb2285689ebdd02944e126f50b625d80a74b314cf8b7cd695db31e96
