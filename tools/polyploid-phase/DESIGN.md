# Polyploid phasing candidate

This is an original, finite synthetic computational-genomics instance, not
an unsolved biological discovery. The former receptor candidates were
solved using exact Pfaffian methods and convex fitting. The new candidate
requires discrete chromosome phasing under allele-dosage constraints.

Scientific basis: weighted minimum error correction for polyploid assembly.
A read originates from one homolog and its mismatch cost is evaluated
against a whole homolog. For a call with weight w, the synthetic symmetric
substitution probability is 1/(1+exp(w)). These weights are log odds;
minimizing their sum corresponds to joint optimization of read assignments
and haplotypes up to constants, not marginalizing the read origins.

References:
- Berger et al., HapTree (2014): https://doi.org/10.1371/journal.pcbi.1003502
- Schrinner et al., Haplotype threading (2020): https://doi.org/10.1186/s13059-020-02158-1

The task explicitly asks for a feasible assembly with weighted MEC <=19885,
not the true biological sequence or a certified global optimum. The cutoff
comes from scoring the originally simulated chromosomes with an independent
direct calculation. It does not come from the reference search output.
The original chromosome witness is author-only. The verifier needs no
reference answer: it recomputes all dosages and the full integer objective.
Global label permutations and alternative valid assemblies are accepted.

Screening (public-input solver, no truth initialization):
- 8 homologs, 227 retained variants, 1800 reads: solved in seconds; rejected.
- 12 homologs, 360 requested variants, 1800 reads: solved in 11.46 seconds.
- 16 homologs, 480 variants, 2400 reads: four simple restarts missed the
  19885 target, best 21542. This does not establish model difficulty.
- An iterated population search on the same public inputs reached 19873
  after 28 restart/perturbation iterations, 177.85 seconds on the host.
  Packaged offline reproduction is required before any frontier trial.

The input-only reference uses genotype-preserving allele and suffix swaps,
read reassignment, and repeated segment perturbations. Its fixed search
seeds are independent of the data-generation seed; it contains no original
chromosomes or precomputed answer.

No model trial has yet established difficulty. Do not ZIP if either required
model solves it; quota interruptions and harness failures are inconclusive.
