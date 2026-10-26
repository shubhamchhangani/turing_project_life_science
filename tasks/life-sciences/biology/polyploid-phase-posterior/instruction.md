# Quantify chromosome-phase uncertainty

A synthetic hexaploid individual has six homologous chromosome copies and
64 ordered biallelic variants. Every variant has alternate dosage three:
three homologs carry allele 1 and three carry allele 0. No chromosome phase
is known. Short overlapping DNA reads provide uncertain information about
phase and about which calls are sequencing errors. Compute the Bayesian
quantities below under the fully specified finite model.

## Probability model

Let $ H_{h,j} $ denote the binary allele on homolog $ h $ at variant $ j $.
At each variant, the six-component allele column is uniform over the
20 binary vectors with exactly three ones. Columns are independent in the
prior. Homologs are labelled 0 through 5 for defining this distribution;
these labels are unobserved. There is no additional linkage prior.

Each read independently chooses an origin $ Z_r $ uniformly from the six
homologs. All eight calls of that read refer to the same origin. Origins
are not observed, and different reads may originate from the same homolog.
For observed call $ a_{r,j} $ with supplied error probability $ e_{r,j} $,

$$ P(a_{r,j}\mid H,Z_r)=\begin{cases}1-e_{r,j},&a_{r,j}=H_{Z_r,j},\\e_{r,j},&a_{r,j}\ne H_{Z_r,j}.\end{cases} $$

Calls are conditionally independent given the complete phase and all read
origins. Errors are symmetric binary substitutions. There are no indels,
unknown alignments, reverse complements, extra chromosome classes, or
parameters to fit. The total probability averages over every allowed phase
and all read origins, with the normalized priors stated above.

## Input

`/root/data/model.json` contains:

- `ploidy`: 6; `variants`: 64; `dosages`: one value, 3, per variant;
- `span`: 8 consecutive variants per read;
- `reads`: 114 records with unique `id` fields `r0001` through `r0114`.
  Each record has `end`, its final zero-based variant index, and `calls`.
  A call has `variant`, `allele` (0 or 1), and `error_probability`.

Variant indices run from 0 through 63. Exactly two reads end at each index
7 through 63. Within each read, all variants from `end-7` through `end` are
observed once. The supplied decimal error probabilities are exact model
parameters. Reads have already been aligned. All supplied calls, including
calls in overlapping reads, are observations to condition on.

## Required quantities

Let $ D $ denote the complete set of observed calls and let $ G $ denote
the supplied dosages. Compute:

1. **Log evidence:** the natural logarithm $ \log P(D\mid G) $, including
   the normalized phase prior and the factor $ 1/6 $ for each read origin.
2. For every read, **error-free probability**
   $ P(N_r=0\mid D,G) $, where
   $ N_r=\sum_{j\text{ called in }r}\mathbf{1}[a_{r,j}\ne H_{Z_r,j}] $.
3. For every read, **expected error count** $ E[N_r\mid D,G] $.

An error is a mismatch to the true allele on the read's latent origin, not
to a consensus sequence. Both posterior quantities condition on the entire
dataset. Homolog-label permutations represent distinct labelled states in
the stated prior; any symmetry reduction must preserve their probability
mass. No particular computational method or explicit phase assembly is
required.

## Outputs

Write `/root/results/evidence.csv` with the literal header `metric,value`
and exactly one data row, whose metric is `log_evidence`.

Write `/root/results/read_posteriors.csv` with the literal header:

```text
read_id,p_error_free,expected_errors
```

Include each of the 114 read IDs exactly once, in any row order. No other
columns, rows, or blank rows are allowed. All numeric fields must be finite
ASCII decimal numbers, optionally with a sign, decimal point and/or an
`e` or `E` exponent. Spaces, thousands separators, Unicode digits, `NaN`
and infinity are not allowed. Probabilities must lie in [0,1] and expected
counts in [0,8]. Use UTF-8 CSV with LF or CRLF endings; standard CSV quoting
is accepted. Both files must be regular files, not symbolic links, each no
larger than 1 MiB.

Absolute errors must be at most $ 2\times10^{-6} $ for log evidence and
$ 5\times10^{-7} $ for every posterior probability and expected count.
Writing at least nine decimal places avoids significant rounding loss.

The environment has one CPU, 2 GiB memory and 10 GiB storage. Python,
NumPy, SciPy, Numba and a C++ compiler are available.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
