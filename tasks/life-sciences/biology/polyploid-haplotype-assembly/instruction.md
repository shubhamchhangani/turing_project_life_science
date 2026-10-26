# Assemble a polyploid chromosome phase

A synthetic polyploid individual has sixteen homologous copies of one
chromosome. The supplied reads cover 480 ordered biallelic variants. Several
homologs share long sequence tracts, and read calls contain substitution
errors. Assemble the homologous allele sequences using the known alternate
allele dosage at every variant.

Use the weighted minimum-error-correction criterion defined below. This is
a constrained assembly problem; the requested result is a feasible assembly
meeting the stated correction budget. Different assemblies can be equally
acceptable. A chromosome label has no biological meaning, and repeated
homolog sequences are allowed.

## Inputs

All inputs are in `/root/data/`:

- `model.json`: the ploidy, number of variants, number of reads, and maximum
  total weighted corrections, **19885**.
- `variants.csv`: columns `variant_id,alt_dosage`. File order is chromosome
  order, from `v0001` to `v0480`. Each dosage is the number of homologs whose
  allele is 1 at that variant; the other homologs carry allele 0.
- `read_calls.csv`: columns `read_id,variant_id,allele,weight`. Every row is
  one observed call. All calls with the same read ID originate from one
  homolog. A read does not switch homologs internally. Calls absent from
  this file are unobserved and contribute no cost. There are 2400 reads.

Alleles 0 and 1 denote the reference and alternate allele. Weights are
positive integers and apply to individual mismatching calls. They are
supplied values, not quantities to infer. The synthetic substitution-error
probability for a call of weight $ w $ is $ 1/(1+\exp(w)) $; consequently,
$ w $ is its log likelihood ratio in favour of a matching rather than a
mismatching call. Reads have already been aligned to the ordered variants;
there are no indels, reverse-complement operations, or unknown alignments.

## Assembly criterion

Let $ H_{h,j} $ be the binary allele on homolog $ h $ at variant $ j $.
Every column must have exactly its supplied alternate dosage:

$$ \sum_{h=1}^{16} H_{h,j}=d_j. $$

For read $ r $, let $ C_r $ be its observed calls. Its correction cost is
its weighted Hamming distance to the best matching *single* homolog:

$$ c_r(H)=\min_{1\leq h\leq16}\sum_{(j,a,w)\in C_r} w\,\mathbf{1}[H_{h,j}\ne a]. $$

The total is $ C(H)=\sum_r c_r(H) $. Supply an assembly satisfying every
dosage constraint and $ C(H)\leq19885 $. The minimum is taken separately for
each read, with ties allowed. Each observed call contributes once to its
read's distance. Do not replace this minimum with a mixture likelihood,
majority vote, or independent per-call homolog choice. No prior, smoothness
penalty, or additional biological constraints enter the criterion.

## Output

Write `/root/results/haplotypes.csv` as a UTF-8 CSV with the literal header:

```text
haplotype_id,alleles
```

Include exactly sixteen data rows, with each ID `H01` through `H16` appearing
once. The `alleles` field must be a string of exactly 480 ASCII characters,
each either `0` or `1`, in the order of `variants.csv`. No whitespace,
separators, missing symbols, extra columns, or blank rows are allowed.
Rows may be in any order. LF or CRLF line endings and standard CSV quoting
are accepted. The output must be a regular file, not a symbolic link, and
no larger than 1 MiB. No read-assignment file or claimed cost is required;
the assembly itself determines every read's cost.

You may use the installed Python, NumPy, SciPy and Numba, or another method.
The available environment has one CPU and 2 GiB memory.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
