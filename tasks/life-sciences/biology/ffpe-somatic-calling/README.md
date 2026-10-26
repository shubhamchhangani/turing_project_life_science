# turing/ffpe-somatic-calling

Somatic SNV calling plus purity, ploidy and tumour-in-normal estimation from a whole-genome-doubled FFPE lymphoma and a matched blood normal that contains circulating lymphoma cells. The data are synthetic: CRAM files are simulated read by read from an explicit physical model with known truth.

## Difficulty

The task copies a real problem in clinical cancer genomics. A lymph-node biopsy is fixed in formalin for years, and the only "normal" available is blood drawn while lymphoma cells circulate (leukaemic-phase mantle cell lymphoma). Analysts in molecular pathology and cancer-genome groups do this analysis, and every standard tool here fails silently on at least one of the following. All five problems must be solved at once, and an early mistake carries through to everything downstream.

1. **FFPE cytosine deamination.** About 20,000 C>T/G>A artefact sites sit at allele fractions of 1–15 %, some supported by 10 or more fragments. Strand-bias filters do not remove them, because each damaged molecule is read on both strands. They are confined to one read-pair orientation (F1R2 for C>T, F2R1 for G>A). True subclonal mutations are themselves C>T-rich (an SBS1-like CpG signature), so dropping low-fraction C>T calls loses recall.
2. **Tumour-in-normal contamination.** Twenty percent of the blood cells are malignant, so every somatic SNV also shows up in the normal. The usual "no alt reads in the normal" rule removes about 88 % of true SNVs. A germline-only normal filter instead lets in the next problem.
3. **Recurrent library and sequencing errors.** 260 sites carry the same error in both libraries at 3–16 %. Of these, 107 are per-read errors that differ between overlapping mates. The other 153 are library-preparation errors fixed into the molecule, so they agree between mates and across PCR duplicates. These sites look exactly like subclonal mutations that are "enriched in blood". Telling them apart needs the copy-number-aware normal/tumour allele-fraction relationship that a fixed tumour-in-normal fraction implies, plus the homopolymer context of the errors.
4. **Whole-genome-doubled tumour with allele-specific copy number.** The segments are 2+2, 2+1, 2+0 (LOH), 2+2 and 3+2, giving ploidy 3.71. Clonal mutations have multiplicity 1, 2 or 3. Diploid assumptions give the wrong purity (about 0.54 from "2 × VAF peak"), the wrong ploidy (2) and a wrong tumour-in-normal fraction. Germline SNPs in LOH and gain regions also leak through naive normal-VAF filters. A doubled solution (ploidy ≈ 7.4, purity ≈ 0.45) fits depth and BAF almost equally well, and only SNV multiplicities reject it.
5. **Unmarked PCR duplicates** (31 % of tumour reads) and overlapping mates (short FFPE fragments). Without collapsing them, single damaged molecules look like multi-read support.

The instruction gives only the clinical and laboratory facts an expert would read: FFPE block age, no UDG, no UMIs, lymphocytosis with a monoclonal CD5⁺ cyclin D1⁺ population in the blood. It never names these effects or how to handle them.

Expert time: about 8–12 hours for an experienced cancer-genomics analyst. The five sub-analyses each need a custom model, because no single off-the-shelf caller handles tumour-in-normal contamination, FFPE orientation bias and allele-specific copy number together.

## Reference solution

`solution/call_somatic.py` (pysam, numpy, scipy; about 3 minutes on 1 CPU and 2 GB) works only from the CRAM files and the reference.

- **Counting.** Duplicate families are collapsed on contig, ends and pair orientation. Overlapping mates are counted once. Allele counts are kept separately for F1R2 and F2R1, and mate-discordant overlap observations are recorded.
- **Copy number.** Tumour depth in 5 kb bins and germline-het BAF are segmented by optimal partitioning. Purity and ploidy come from an ASCAT-style grid fit with integer allele-specific states and a self-consistent ploidy constraint. Competing integer solutions are separated by how close the clonal-SNV multiplicities, VAF·(P·n + 2(1 − P))/P, sit to integers. The doubled solution loses decisively (score 19 vs 33).
- **FFPE damage.** A nonparametric (grid EM) prior on the per-site deamination rate is fitted on background C/G sites. For C>T and G>A, the clean orientation gives the main evidence and the damage-prone orientation is marginalised over that prior.
- **Classification.** A two-class mixture compares a somatic model with a library-independent recurrent-error model. In the somatic model, normal VAF = T·x/(T·n + 2(1 − T)), with x the mutant copies per malignant cell. In the recurrent-error model, normal VAF ≈ tumour VAF up to a lognormal factor. Both models are integrated over the true allele fraction. T is the maximiser of the mixture likelihood. Homopolymer context enters as a class-specific feature learned by EM. Mate-discordant sites are removed as per-read errors.

Result on the shipped data: 404 calls, precision 0.958, recall 0.968, purity 0.620, ploidy 3.710, normal tumour fraction 0.195.

## Verification

`tests/test_outputs.py` runs in a separate, offline container. Its seven tests, each with a docstring, read the two artifacts.

- **Contract checks first.** Each failure names the violated field:
  - The SNV file must be a regular UTF-8 file of at most 1 MiB, with the exact header and four fields per row.
  - `chrom` must be in a closed set, and `pos` must be a plain integer within the chromosome.
  - `ref` must equal the base in `reference.fa`, and `alt` must differ from it.
  - There must be no duplicate positions and 1–5,000 rows.
  - The JSON must have exactly three finite numeric keys in their stated ranges.
- **Accuracy checks.**
  - SNV precision ≥ 0.90 and recall ≥ 0.85, with exact chrom/pos/ref/alt matching against all 400 planted SNVs.
  - Purity and ploidy within 10 % relative error.
  - Normal tumour fraction within 15 % relative error.
  - The reward is 1 only if every test passes (`reward.txt`). Pytest-ctrf writes a per-test report.

**Independence of expected values.** The truth set and parameters are the values injected by the generator (truth by construction). None of them come from the reference solution.

**Tolerance calibration.** The table shows the reference solution on four independent round-2 simulations (seed 8 is the shipped data):

| seed | precision | recall | purity err | ploidy err | T err |
|---|---|---|---|---|---|
| 8 | 0.958 | 0.968 | 0.0 % | 0.0 % | −2.5 % |
| 61 | 0.979 | 0.955 | −1.6 % | 0.0 % | −3.7 % |
| 71 | 0.970 | 0.955 | 0.0 % | 0.0 % | −1.2 % |
| 81 | 0.955 | 0.960 | −1.6 % | 0.0 % | −1.2 % |

Thresholds leave at least 0.055 margin in precision, 0.10 in recall, 8 % in purity and 11 % in T. Recall 0.85 allows for subclone-B SNVs (cancer-cell fraction 0.22, VAF about 0.035–0.06) that no method can separate from noise. Each simplified pipeline below fails at least one threshold (seed 8):

| baseline | precision | recall |
|---|---|---|
| strict normal (≤ 1 alt read), no FFPE handling | 0.005 | 0.125 |
| orientation filter + strict normal | 0.923 | 0.120 |
| orientation + germline-only normal (VAF < 0.3) | 0.631 | 0.932 |
| T/N ratio filter, no FFPE handling | 0.028 | 0.785 |
| T/N ratio, drop low-VAF C>T/G>A | 0.900 | 0.542 |
| T/N ratio + orientation (diploid thresholds) | 0.917 | 0.772 |
| orientation + mate-discordance + pooled support (no ratio model) | 0.696 | 0.983 |

Diploid purity/ploidy estimates (ploidy 2, T ≈ 0.25) fail the ploidy and T tests.

**Other checks.**
- **Harbor:** oracle reward 1.0, nop reward 0.0.
- **Offline oracle run:** I ran the oracle offline (`docker run --network none --cpus 1 --memory 2g` on the agent image). It finished in 174 s, wrote both artifacts and passed all seven tests.
- **Adversarial regression suite (35 cases):** all 4 valid variants (CRLF, no final newline, integer-valued numbers, the oracle) pass. All 31 invalid ones fail, including truth plus recurrent artefacts, truth plus germline SNPs, clonal-only calls, wrong alt alleles, a diploid ploidy, the doubled solution, NaN, strings or booleans, unicode digits, leading zeros, duplicates, FIFO files and files over 1 MiB.

**Agent trials.** All trials used Claude Opus 5.5 in isolated containers (1 CPU, 2 GB, internet, 18,000 s budget), with solution, tests and truth unreachable during the runs.

- **Round 1 (earlier data, all recurrent errors per-read).** Trial 1 passed. Trial 2 failed with precision 0.57 (212 recurrent errors and 63 germline SNPs in copy-number-altered regions reported).
- **Hardening.** Library-preparation (molecule-level) recurrent errors were added, giving the shipped data.
- **Round 2 (shipped data).** Both trials failed, each on precision only:
  - Trial 3: precision 0.871, recall 0.958, all three parameters correct. It reported 53 recurrent errors.
  - Trial 4: precision 0.609, recall 0.960, all three parameters correct. It reported 247 recurrent errors as "blood-enriched subclonal mutations".
- **Transcript reading.** Both agents understood the task, used every input and solved the copy-number, FFPE and contamination parts. They failed on the scientific judgement of what a mutation that is "too frequent in the blood" means under a fixed tumour-in-normal fraction.
- **Instruction change after round 2.** After round 2 one sentence was added ("Non-malignant cells carry the germline genome only, with no somatic mutations"). It closes the explanation that trial 3 offered (mutations in non-malignant blood cells), which was already excluded by the requirement that the reported SNVs be carried by malignant cells.

## Data generation

The generator is outside the package. It uses one cell (and one haplotype copy) per DNA fragment, with fragments drawn in proportion to DNA content c·CN + 2(1 − c).

| Component | Settings |
|---|---|
| Genome | 4 synthetic chromosomes, 1 Mb in total, with CpG depletion |
| Germline SNVs | heterozygous 1/1100 bp, homozygous 1/3500 bp |
| Sample fractions | purity P = 0.62; blood malignant fraction T = 0.20 with the same clonal composition |
| Somatic SNVs | 230 clonal (multiplicity 1 or n_hap); 110 in subclone A (CCF 0.35); 60 in subclone B (CCF 0.22, nested in A) |
| FFPE damage | tumour only; per-site rates ~ Gamma(0.45, mean 0.01) on single original strands |
| PCR duplicates | geometric family sizes (tumour mean 1.45, normal 1.11) |
| Fragment lengths | lognormal; tumour median 210 bp, normal median 320 bp |
| Base qualities | binned Illumina qualities (2/12/23/37) with per-cycle decay and matching error rates |
| Unique depth | tumour mean 170×, normal mean 90× |
| Recurrent errors | at bases following homopolymer runs ≥ 5 |

The generator places all reads correctly and writes them as CRAM 3.0. The mutation density (400/Mb) is deliberately higher than in real mantle cell lymphoma, to give the precision and recall statistics power.

## References

- Costello M. et al. (2013) Discovery and characterization of artifactual mutations in deep coverage targeted capture sequencing data due to oxidative DNA damage during sample preparation. *Nucleic Acids Res* 41:e67.
- Do H., Dobrovic A. (2015) Sequence artifacts in DNA from formalin-fixed tissues: causes and strategies for minimization. *Clin Chem* 61:64–71.
- Taylor-Weiner A. et al. (2018) DeTiN: overcoming tumor-in-normal contamination. *Nat Methods* 15:531–534.
- Van Loo P. et al. (2010) Allele-specific copy number analysis of tumors (ASCAT). *PNAS* 107:16910–16915.
- Carter S.L. et al. (2012) Absolute quantification of somatic DNA alterations in human cancer (ABSOLUTE). *Nat Biotechnol* 30:413–421.
- Cibulskis K. et al. (2013) Sensitive detection of somatic point mutations in impure and heterogeneous cancer samples (MuTect). *Nat Biotechnol* 31:213–219.
- Nakamura K. et al. (2011) Sequence-specific error profile of Illumina sequencers. *Nucleic Acids Res* 39:e90.
- Ott G. et al. (1997) Blastoid variants of mantle cell lymphoma: frequent bcl-1 rearrangements at the major translocation cluster region and tetraploid chromosome clones. *Blood* 89:1421–1429.
- Alexandrov L.B. et al. (2020) The repertoire of mutational signatures in human cancer. *Nature* 578:94–101.
