
# turing/metabolic-flux-analysis

Steady-state ¹³C metabolic flux analysis (13C-MFA) of a cancer-cell-like central carbon network from two parallel tracer experiments, [1,2-¹³C₂]glucose and [U-¹³C₅]glutamine, measured as raw GC-MS mass isotopomer distributions of TBDMS-derivatized fragments.

## Difficulty

13C-MFA is the standard way to measure intracellular fluxes in metabolic engineering and cancer metabolism (for example reductive carboxylation, oxidative PPP use and anaplerosis). In practice it is done by specialist postdocs and PhD students with dedicated software (INCA, 13CFLUX2, OpenFLUX). The agent receives only the reaction list, the measured ions and the raw spectra. It must build the whole analysis itself, and a mistake at any step biases the fluxes:

1. **Atom transitions.** These are not provided. The agent must derive the carbon fate for 31 reactions, including the triose split, the transketolase and transaldolase carbon transfers, the decarboxylations (PDH, oxidative PPP, IDH, AKGDH, malic enzyme) and the stereospecific citrate synthase mapping.
2. **Symmetric intermediates.** Succinate and fumarate scramble label, while citrate does not.
3. **Natural-isotope correction for the derivative.** About 75 % of the carbon in each ion comes from TBDMS groups, plus Si, O, N and H isotopes. Raw spectra cannot be compared with backbone labeling directly. Unlabeled substrate carbon also carries natural ¹³C, which must not be "corrected away".
4. **Parallel-tracer fitting.** Neither tracer alone constrains every graded flux (see the negative controls below). The fit is a non-convex weighted least-squares problem with 18 parameters, including poorly identifiable exchange fluxes, so it needs multiple starting points.
5. **Reductive carboxylation.** Detecting it requires keeping IDH reversible and recognizing the M+5 citrate signature in the glutamine tracer.

**Data:** synthetic, generated from known fluxes (see *Data generation*), so the true answer is known exactly. The network, tracers, fragment ions and noise level mirror published parallel-labeling studies of cancer cells.

**expert_time_estimate_hours = 12.** An experienced 13C-MFA practitioner who writes the model from scratch, or encodes it in INCA, would need about a day and a half to build, check and fit it.

## Reference solution

`solution/mfa_fit.py` encodes the standard atom transitions for every reaction in `network.tsv` and simulates both experiments with a full-isotopomer steady-state model. Every metabolite carries a probability vector over its 2ⁿ carbon isotopomers, solved by fixed-point iteration with renormalization, and Suc/Fum are symmetrized. The simulated backbone MID of each fragment is convolved with the natural-isotope distribution of all non-backbone atoms in the ion formula, truncated to M0…M(n+2) and renormalized, and then compared with the raw data. Fitting raw data this way is equivalent to correcting it first.

Net fluxes are parameterized on the stoichiometric null space, with 10 free fluxes plus glucose = 100. The 8 reversible reactions have bounded exchange fluxes, and irreversibility is enforced by a penalty term. The weighted residuals (MIDs and extracellular rates) are minimized with `scipy.optimize.least_squares` from 24 seeded starting points. On the task data, all starting points converge to SSR = 133.1 for 156 measurements and 18 parameters, which lies inside the 95 % χ² range for 138 degrees of freedom.

## Verification

`tests/test_outputs.py` runs in a separate no-network container. It reads `fluxes.csv`, injected truth in `expected_fluxes.json`, and trusted copies of the network and measured data baked into the verifier image. `forward_model.py` comes from the independently cross-checked data generator, not the reference solution. It predicts measurements without fitting parameters or loading a reference flux map.

1. **Contract:** literal header; one row for each of the 31 reactions; finite ASCII decimal fields without whitespace or underscores. Signs, decimal points, and decimal exponents are allowed.
2. **Consistency:**
   - GLCup = 100.
   - forward − reverse = net within 0.05.
   - Irreversible reactions have reverse ≈ 0 and net ≥ 0.
   - Every intracellular metabolite balance closes within 0.1.

   These checks make sure the agent submits a real steady-state flux map. Both
   one-way fluxes must be nonnegative without a negative allowance. The 0.05
   directional allowance accommodates independent rounding of net, forward and
   reverse fluxes: reporting each to two decimal places gives at most 0.015
   error in their three-term identity. The 0.1 balance allowance accommodates
   accumulation of that rounding across reactions: the largest sum of absolute
   stoichiometric coefficients for an intracellular metabolite is 6 (T3P),
   so two-decimal rounding contributes at most 6 × 0.005 = 0.030.
   The 0.1 allowance leaves margin above this rounding bound. These
   absolute allowances concern reporting precision, not measurement noise or
   relaxation of the steady-state model. Glucose normalization separately uses
   an absolute tolerance of 0.0001.

   Output reads use a nonblocking descriptor and require a regular file of at
   most 1 MiB, preventing FIFOs and unbounded files from hanging the verifier.
   The test runner has a 100-second deadline and a 5-second kill grace period,
   within the verifier's 120-second budget.
3. **Accuracy:** six fluxes are compared with the **injected true values** (truth by construction, not oracle output): G6PDH, PDH, PC, ME, CS net, and IDH reverse (reductive carboxylation).

4. **Full measurement fit:** weighted SSR across all 148 raw MID fractions and
   8 positive-SD extracellular rates must be at most 180. Predictions use the
   reported net fluxes and exchange rates reconstructed as the smaller one-way
   rate, resolving independent rounding consistently with the instructions.
   The forward prediction runs in a separate Python process with a 60-second
   timeout; it executes only trusted verifier code and receives validated JSON.

**Measurement-fit threshold rationale.** With 156 independent Gaussian
observations and 18 fitted parameters, the nominal residual degrees of freedom
are 138. The 99th percentile of chi-square with 138 degrees of freedom is
179.5611, rounded up to 180. This is an explicit statistical goodness-of-fit
criterion, not a guarantee of a unique maximum-likelihood map. The existing
six injected-truth tolerances remain separate requirements. The threshold was
derived from the noise model rather than from reference-solution output.
Reversible exchange rates can be weakly identified: increasing PGI exchange by
1000 in the saved reference map changes SSR from 133.09 to 134.89, which remains
acceptable under this declared criterion. Larger changes in other exchanges
can make the same net map inconsistent with the observations. This distinction
is why the instructions allow statistically equivalent maps instead of
requiring exact equality to a single fitted map.

**Verifier regression checks (2026-10-23, version 1.1.0).** The rebuilt verifier
image with pinned NumPy 1.26.4 was run without network, with 1 CPU and 2048 MiB.
The saved genuine reference output passed all 14 tests with reward 1; the shell
runner completed in 0.19 seconds. Its forward-model SSR was 133.0892. The
injected true map also passed. Author-side regressions are retained outside the
task in `tools/metabolic-flux-analysis/check_verifier_regressions.py`.
Checks cover literal LF/CRLF headers, quoted and whitespace-altered headers,
negative one-way fluxes, missing or duplicate rows, nonfinite and overflowing
numbers, numeric spaces and underscores, Unicode digits, valid signed decimals
and exponents, FIFOs, oversized files, and absent or empty output.
Adding 1000 to both one-way fluxes of TK1, TK2 or GLUDH preserves every net flux
and all six accuracy targets, but the new measurement-fit requirement rejects
these maps. The scientific optimizer was not rerun; its saved genuine output
was revalidated against the complete revised verifier.

**Independent validation of the data generator.**
- The isotopomer simulator (`tools/metabolic-flux-analysis/mfa_model.py`, kept outside the task) was compared at the true fluxes with **FreeFlux 0.3.8** (Wu et al. 2023), an independent EMU-based package. All 22 fragment backbone MIDs agree to within 4 × 10⁻⁶.
- The derivatization natural-abundance model was compared with **IsoCor 2.2.2** (Millard et al. 2019). Correcting the simulated raw MIDs recovers the backbone MIDs to within 4 × 10⁻⁸.

**Tolerance calibration.** 48 independent noise realizations were generated and fitted with the reference method. Each tolerance is at least 4 SD of the correct-method scatter and at least 1.3 × the worst error observed:

| Flux | True | SD over 48 fits | Worst error | Tolerance |
|---|---|---|---|---|
| G6PDH net | 15 | 1.74 | 4.6 | 45 % |
| PDH net | 42 | 1.12 | 2.9 | 12 % |
| PC net | 11 | 0.26 | 0.87 | 12 % |
| ME net | 40 | 0.62 | 1.75 | 8 % |
| CS net | 30 | 0.57 | 1.34 | 10 % |
| IDH reverse | 8 | 0.32 | 0.86 | 20 % |

**Negative controls.** Previously saved outputs from fits with one realistic mistake were reevaluated using the independent forward model. The current measurement-fit requirement rejects all eight mistaken fits:

| Analysis | Measurement SSR | Required maximum |
|---|---:|---:|
| Correct | 133.09 | 180 |
| No derivatization correction | 23824.81 | 180 |
| Unlabeled carbon treated as pure ¹²C | 527.91 | 180 |
| No succinate/fumarate symmetry | 256.00 | 180 |
| [1,2-¹³C]glucose data only | 257.91 | 180 |
| [U-¹³C]glutamine data only | 3007.45 | 180 |
| Naive citrate synthase atom map | 2785.07 | 180 |
| IDH irreversible | 1042.20 | 180 |
| Tracer purity 99 % taken as 100 % | 260.70 | 180 |

The tracer-purity simplification passed the earlier six-target-only verifier,
but fails the new full-data fit requirement.

**Historical checks before version 1.1.0 (Harbor 0.24.0, 2026-10-22):** oracle reward 1.0 and nop reward 0.0. **Offline oracle run:** `solve.sh` ran in the built agent image with `docker run --network none --cpus 1`. It finished without error and its `fluxes.csv` passed all 13 verifier tests (reward 1). This check passed.

**Frontier-agent run:** gpt-5.6-sol scored 0/2 for scientific reasons; see *Agent trials* below.

## Data generation

`tools/metabolic-flux-analysis/make_dataset.py`, seed 20261022:
- True net fluxes are in `true_fluxes.py`, with glucose = 100. Key values: oxPPP 15, PDH 42, PC 11, ME 40, CS 30, glutamine uptake 40, lactate 167.
- Exchange fluxes: PGI 60, TK1 15, TA 25, TK2 10, IDH 8, FH 80, MDH 40, GLUDH 50.
- MIDs are simulated exactly, then Gaussian noise (SD 0.003 per fraction) is added. Extracellular rates get 5–10 % relative noise.
- Fragment ions are the standard TBDMS [M−57]⁺, [M−85]⁺, [M−159]⁺ and [M−189]⁺ fragments (Antoniewicz et al. 2007).
- Isotope abundances are IUPAC representative values.

## Agent trials

Codex CLI with **gpt-5.6-sol**, 2 trials (Harbor 0.24.0, 2026-10-22). **Mean reward 0.0.** Both trials passed every contract and consistency check and failed on flux accuracy.

| Trial | Reward | Graded fluxes failed | Agent CS / PDH | SSR of the agent's map (exact model) |
|---|---|---|---|---|
| 1 | 0 | PDH, CS | 37.3 / 48.0 | 15,895 |
| 2 | 0 | CS | 34.6 / 45.8 | 2,561 |
| Reference | 1 | none | 29.4 / 40.3 | 133 |

**Transcript analysis.**
- **Trial 1.** The agent read all inputs and correctly removed the derivative's natural isotopes. It then skipped a statistical fit: it estimated oxPPP, malic enzyme and reductive IDH by inspection and closed the network by fixing every extracellular rate at its measured value. Measured lactate secretion (160.3 ± 8.4) is lower than the true 167, so treating it as exact pushed the missing carbon into PDH and the TCA cycle.
- **Trial 2.** The agent built a full joint isotopomer model with correct derivatization handling and correct atom transitions, apart from a C1/C2 swap in PGI that is harmless for this tracer. Two problems remained:
  - Its isotopomer solver used capped Jacobi iteration without renormalizing, so it simulated labeling inaccurately.
  - It accepted a fit whose own SSR was 331 with a large gradient (optimality ≈ 10⁵), above the χ² acceptance bound (≈ 171 for 138 degrees of freedom).

  Evaluated exactly, its flux map is rejected by the data (SSR 2,561 compared with 133).

**Task-side causes ruled out.** Neither agent reported missing information or misread the contract. The container and data worked in both trials. The failed quantity (CS) is tightly determined by the data: any map near the agent's value fits 19× worse than the optimum, so the tolerance did not reject a defensible answer. The failures come from 13C-MFA methodology: weighting measurements by their errors, converging the isotopomer balances, and goodness-of-fit testing.
## References

- Antoniewicz, Kelleher & Stephanopoulos (2007) Elementary metabolite units (EMU). *Metab. Eng.* 9:68–86.
- Antoniewicz, Kelleher & Stephanopoulos (2007) Accurate assessment of amino acid mass isotopomer distributions for metabolic flux analysis. *Anal. Chem.* 79:7554–7559.
- Crown & Antoniewicz (2013) Parallel labeling experiments and metabolic flux analysis. *Metab. Eng.* 16:21–32.
- Metallo et al. (2012) Reductive glutamine metabolism by IDH1 mediates lipogenesis under hypoxia. *Nature* 481:380–384.
- Young (2014) INCA: a computational platform for isotopically non-stationary metabolic flux analysis. *Bioinformatics* 30:1333–1335.
- Wu et al. (2023) FreeFlux: a Python package for time-efficient isotopically nonstationary metabolic flux analysis. *ACS Synth. Biol.* 12:2707–2714.
- Millard et al. (2019) IsoCor: isotope correction for mass spectrometry labeling experiments. *Bioinformatics* 35:4484–4487.
