A patient with pleomorphic mantle cell lymphoma underwent an excisional lymph-node biopsy, and a blood sample was drawn on the same day as the matched normal. Both were whole-genome sequenced. Determine the somatic single-nucleotide variants (SNVs) of the lymphoma and three sample-level quantities: tumour purity, tumour ploidy and the fraction of malignant cells in the normal sample.

## Specimens

- **Tumour** (`tumor`): formalin-fixed, paraffin-embedded (FFPE) lymph-node biopsy. The block was stored for 8 years before extraction. The DNA was not treated with uracil-DNA glycosylase, and the library has no unique molecular identifiers.
- **Normal** (`normal`): fresh-frozen peripheral-blood mononuclear cells. The absolute lymphocyte count at the draw was 21 × 10⁹/L. Flow cytometry of the same draw detected a monoclonal CD5⁺ cyclin D1⁺ B-cell population.
- Both libraries are Illumina paired-end 2 × 150 bp. Reads are aligned to the reference below and coordinate-sorted. Duplicates are **not** marked. Base qualities are binned.

## Input files (all in `/root/data/`)

- `reference.fa` and `reference.fa.fai`: the reference genome. It is synthetic and has four chromosomes, `chr1` (300,000 bp), `chr2` (250,000 bp), `chr3` (250,000 bp) and `chr4` (200,000 bp), 1,000,000 bp in total. The patient's germline genome is diploid on all four chromosomes.
- `tumor.cram` and `tumor.cram.crai`: tumour alignments (read group sample `tumor`).
- `normal.cram` and `normal.cram.crai`: normal alignments (read group sample `normal`).

Decode both CRAM files with `/root/data/reference.fa`.

## Biological model of the samples

The data were simulated from a physical model with known truth. These statements hold exactly:

- Every cell in either specimen is either a non-malignant diploid cell or a malignant cell of this lymphoma. Non-malignant cells carry the germline genome only, with no somatic mutations.
- All malignant cells share one clonal copy-number profile, with no subclonal copy-number changes. Copy number is constant within each segment and segment boundaries are not given. Non-malignant cells have no copy-number changes.
- Malignant cells may carry subclonal SNVs as well as clonal ones. The malignant cells in the blood have the same clonal composition as the malignant cells in the lymph node.
- The only genomic differences are germline SNVs, somatic SNVs and the tumour copy-number changes. There are no indels or other structural variants, and all reads are correctly placed.

## Output 1: `/root/results/somatic_snvs.tsv`

A UTF-8, tab-separated text file of at most 1 MiB. The first line is exactly `chrom	pos	ref	alt` (four tab-separated column names). Each further line is one somatic SNV, in any order, with these four fields:

- `chrom`: one of `chr1`, `chr2`, `chr3`, `chr4`.
- `pos`: the 1-based position on that chromosome, written as a plain decimal integer with no sign, leading zeros or whitespace.
- `ref`: the reference base at that position in `reference.fa`, as one uppercase letter from `A`, `C`, `G`, `T`.
- `alt`: the somatic allele, as one uppercase letter from `A`, `C`, `G`, `T` that differs from `ref`.

Report between 1 and 5,000 SNVs, with each (`chrom`, `pos`) at most once. A somatic SNV is a single-nucleotide change that is carried by malignant cells, clonally or subclonally, and is absent from the patient's germline. Do not report germline variants, DNA damage, or library or sequencing errors.

## Output 2: `/root/results/sample_estimates.json`

A UTF-8 JSON object of at most 1 MiB with exactly these three keys. Each value is a finite JSON number, not a string or boolean.

- `tumor_purity`: the fraction of all cells in the tumour specimen that are malignant, in (0, 1].
- `tumor_ploidy`: the mean total copy number of a malignant cell, averaged over all 1,000,000 reference positions with equal weight. It must be in (0, 10].
- `normal_tumor_fraction`: the fraction of all cells in the normal (blood) specimen that are malignant, in [0, 1).

## How the result is judged

The outcome is pass or fail. Every requirement below must hold, and the truth is the set of somatic SNVs and the sample parameters used to simulate the data.

1. **Contract.** Both files exist and follow the formats above exactly.
2. **SNV precision.** At least 90 % of the reported rows match a true somatic SNV. A match requires `chrom`, `pos`, `ref` and `alt` to all be equal.
3. **SNV recall.** At least 85 % of all true somatic SNVs, clonal and subclonal, are reported with matching `chrom`, `pos`, `ref` and `alt`.
4. **Purity.** `tumor_purity` is within 10 % of the true value (relative error |reported − true| / true ≤ 0.10).
5. **Ploidy.** `tumor_ploidy` is within 10 % of the true value (relative error ≤ 0.10).
6. **Normal contamination.** `normal_tumor_fraction` is within 15 % of the true value (relative error ≤ 0.15).

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
