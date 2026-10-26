# Exact phase uncertainty, after two rejected candidates

The spatial-ligand task was solved by Claude Opus 5.5 low. A new weighted-MEC
polyploid assembly task was then solved by GPT-6-sol max in 236.89 seconds,
reward 1, cost 19878 versus a public budget of 19885. Neither was packaged.

The current task asks for normalized Bayesian evidence and per-read
posterior uncertainty, not another feasible assembly or larger MEC instance.
It sums the stated finite model over all phases and read origins. Basic
independent-read estimates differ from full posterior probabilities by as
much as 0.6385 and from expected error counts by 1.7041 on this dataset.
A MAP phase also does not determine the requested marginal probabilities.

The source paper motivating uncertainty-aware polyploid phasing is
https://arxiv.org/abs/2512.04393 (Hosseini et al.). The model and data are
original; no claim is made that the task reproduces that paper's method.

Exact reference algorithm: canonical sorted homolog allele-history rows,
with uniform new dosage-compatible columns and correct orbit multiplicity.
The seven-column history has 1,822,528 orbits. Its count is independently
verified by Burnside's lemma, and orbit masses sum to 20^7 labelled histories.
Forward/backward messages yield local read-error posterior quantities.
Disk checkpoints prevent all forward messages from occupying RAM at once.

Independent expected-value source: a separately written forward-only
partition evaluator. For each read it replaces the emission by either the
zero-error contribution or the mismatch-weighted contribution, then computes
the ratio of that partition function to the ordinary evidence. This source
has no reference posterior or backward messages. Sharing the mathematically
defined transition catalog does not supply it with any expected answer.

Six exhaustive labelled-state checks (160000 phases each) agree within
1.2e-12 for posterior quantities. The full reference takes 11.39 seconds on
the host; isolated resource reproduction and full independent expectations
must finish before model trials. No frontier difficulty is established.
