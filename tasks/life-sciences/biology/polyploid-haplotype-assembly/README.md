# Polyploid haplotype assembly

Original synthetic computational-genomics task: reconstruct sixteen phased
homolog sequences from overlapping noisy reads and exact allele dosages.
The contract asks for a feasible weighted minimum-error-correction assembly,
not recovery of a uniquely identifiable hidden genome or proof of optimality.

## Scientific basis

Polyploid chromosome phasing with ambiguous read origins and identical
shared tracts is a studied assembly problem. This task uses weighted MEC:
each read is assigned to one entire homolog, with minimum weighted mismatch
cost. For a call of integer weight w, the symmetric substitution-error
probability is 1/(1+exp(w)), so the mismatch weight is the log likelihood
ratio. This is joint read-assignment optimization, not a mixture likelihood.

- Berger et al., HapTree (2014), https://doi.org/10.1371/journal.pcbi.1003502
- Schrinner et al., Haplotype threading (2020), https://doi.org/10.1186/s13059-020-02158-1

## Data generation and independent target

Author generator: `tools/polyploid-phase/prototype.py`, function `generate`,
seed 20261031, ploidy 16, 480 variants and 2400 reads. Independent random
ancestral allele sequences supply mosaic 30-variant tracts to each homolog.
Monomorphic columns are removed (none for this instance). Reads span
25--65 ordered variants, omit calls independently with probability 0.18,
and carry integer weights sampled from 1, 2, 3 with probabilities
0.25, 0.50, 0.25. Calls are flipped with the stated error probability.

The original simulated chromosome phase has weighted correction cost 19885,
computed directly from all calls. That independent feasibility witness
sets the public correction budget. The threshold is an integer, so there
is no floating-point tolerance. The source phase is not a required answer:
any assembly satisfying every dosage and the public budget is valid.
Original chromosomes and generation records are retained outside this task.

The input-only reference reads only the public files. It uses dosage-
preserving local moves, exact read reassignment and iterated segment
perturbations. It neither reads the original phase nor embeds a solution.

## Verification and isolation

The separate offline verifier uses its own immutable copies of all inputs.
It never executes agent code and does not import reference search routines.
It first validates the full literal output contract, then every dosage, then
the complete weighted correction cost of every read. Cost computation uses
exact Python integers. The read's best homolog is selected after summing
all its calls, not independently for each call.

Output opening is nonblocking, disallows symlinks, checks the file type and
limits reading to 1 MiB. Raw sequence fields permit only the specified
number of ASCII binary characters; whitespace and Unicode digits fail.
The runner enables errexit and has a 100-second verifier invocation limit.

Author checks include exhaustive small-instance scoring comparisons, valid
homolog-label permutations and CSV representations, malformed artifact
rejection, and dosage-correct assemblies with poor read fit. The original
witness exercises acceptance exactly on the specified objective boundary.

## Candidate status

The packaged input-only reference ran offline with one CPU and 2 GiB in
175.52 seconds, producing cost 19873. The separate offline verifier passed
all three tests in 0.10 seconds. Twelve exhaustive small scoring checks,
five valid artifact variants and 26 invalid artifact cases passed.

Draft: not approved for submission. Fresh isolated frontier-model trials
are required. A model solving this
candidate disqualifies it for the requested difficulty criterion. No claim
is made that this is an unsolved biological problem or that every model
will fail it. No submission ZIP has been created.

## Completed model result

Rejected for the requested difficulty criterion: GPT-6-sol max solved this
candidate, reward 1, all three tests passed in 236.89 seconds of agent
execution. Its correction cost was 19878. Claude's unfinished trial was
cancelled after GPT succeeded and does not count as a failure. No ZIP was
created. The pre-trial package manifest remains in the author evidence.
