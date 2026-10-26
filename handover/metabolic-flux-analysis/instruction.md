A proliferating cancer cell line was grown in two parallel cultures that are identical except for the isotopic tracer, and harvested at metabolic and isotopic steady state. Intracellular metabolites were measured by GC-MS as TBDMS derivatives. Determine the intracellular flux map of central carbon metabolism by ¹³C metabolic flux analysis using both experiments together. The fluxes are the same in both cultures.

## Tracer experiments

- `glc12`: medium glucose is 100 % [1,2-¹³C₂]glucose. C1 and C2 are each independently 99.0 atom % ¹³C; C3–C6 have natural ¹³C abundance. Glutamine is unlabeled.
- `gln5`: medium glutamine is 100 % [U-¹³C₅]glutamine, with each carbon independently 99.0 atom % ¹³C. Glucose is unlabeled.

"Unlabeled" carbon always carries natural-abundance ¹³C at the fraction given for `C` in `/root/data/isotope_abundances.csv`. This includes all carbon in unlabeled glucose and glutamine, and the CO₂ fixed by carboxylation reactions.

## Input files (all in `/root/data/`)

- `network.tsv` (tab-separated; columns `reaction_id`, `equation`, `reversible`): the complete metabolic network of 31 reactions. Use exactly these reactions, metabolites and reversibilities.
  - Species ending in `.ext` are extracellular sources or sinks. `CO2` is an external pool and is not balanced.
  - Carboxylating reactions (`PC`, and `IDH` in its reverse direction) fix CO₂ with natural-abundance labeling. CO₂ released inside the cell does not re-enter.
  - `T3P` is the lumped, isotopically equilibrated triose-phosphate pool (DHAP + GAP). `P5P` is the lumped pentose-5-phosphate pool (Ru5P, R5P, X5P). `PG3` is 3-phosphoglycerate.
  - The model has a single compartment and no kinetic isotope effects. Succinate (`Suc`) and fumarate (`Fum`) are rotationally symmetric molecules. No other metabolite is treated as symmetric.
  - Carbon atom transitions are not given; they follow the standard biochemistry of each reaction.
- `fragments.csv` (columns `fragment`, `metabolite`, `backbone_carbons`, `ion_formula`): one row per measured GC-MS fragment ion.
  - `backbone_carbons` lists, separated by `;`, which carbons of the metabolite are contained in the ion. Carbons are numbered by IUPAC convention, so C1 is the carboxyl carbon of the amino and carboxylic acids.
  - `ion_formula` is the full elemental formula of the measured ion. It includes the metabolite backbone and all derivatization atoms.
- `isotope_abundances.csv` (columns `element`, `mass_shift`, `abundance`): natural isotope abundances to use for C, H, N, O and Si. These values are authoritative.
- `mid_measurements.csv` (columns `experiment`, `fragment`, `isotopologue`, `fraction`, `sd`): raw mass isotopomer distributions (MIDs), **not corrected for natural isotopes**.
  - For a fragment with $ n $ backbone carbons, isotopologues $ M_0 $ through $ M_{n+2} $ are reported as fractions of the sum over those $ n+3 $ peaks.
  - Each fraction carries additive Gaussian measurement noise with the given standard deviation. Values were not renormalized after the noise was added.
- `extracellular_rates.csv` (columns `reaction_id`, `rate`, `sd`): measured net rates of uptake and secretion reactions, with standard deviations.
  - All fluxes are normalized to glucose uptake, `GLCup` = 100 exactly.

## Output — `/root/results/fluxes.csv`

A regular UTF-8 comma-separated file no larger than 1 MiB, with header `reaction_id,net_flux,forward_flux,reverse_flux` and exactly 31 data rows, one for every reaction in `network.tsv`, in any order.

- All flux values are finite ASCII decimal numbers in the same units (glucose uptake = 100), with no whitespace inside a field. An optional sign, decimal point, and base-10 exponent are allowed (for example `100`, `-0.05`, `.5`, or `1e2`); underscores and other numeric formats are not allowed.
- `net_flux` is the net steady-state flux in the direction the equation is written.
- `forward_flux` and `reverse_flux` are the one-way fluxes, both nonnegative, with $ forward\_flux - reverse\_flux = net\_flux $.
- For irreversible reactions, `reverse_flux` is zero within the 0.05 consistency tolerance below.
- Report a flux map satisfying the consistency, accuracy, and measurement-fit requirements below. Uncertainty estimates are not required. Flux maps that are indistinguishable within the stated tolerances are acceptable; a unique global optimum is not required.

## Output requirements

The output is valid only when all of the following requirements hold.

1. **Contract.** The file exists with the exact header, exactly one row for each of the 31 reactions, and finite numbers.
2. **Consistency.**
   - `GLCup` `net_flux` is 100, within an absolute tolerance of 0.0001.
   - Both one-way fluxes must be nonnegative. For every reaction, $ |forward - reverse - net| \leq 0.05 $.
   - Irreversible reactions have `reverse_flux` within 0.05 of 0 and $ net\_flux \geq -0.05 $.
   - The reported net fluxes close the steady-state balance of every intracellular metabolite (all species except `.ext` species and `CO2`) to within 0.1.
3. **Accuracy.** For each of the following six quantities, $ |reported - true| \leq tolerance \times |true| $, where $ true $ is the flux value used to produce the data and $ tolerance $ is the stated relative tolerance:

| Quantity | Column | Relative tolerance |
|---|---|---|
| `G6PDH` (oxidative pentose phosphate pathway) | `net_flux` | 45 % |
| `PDH` | `net_flux` | 12 % |
| `PC` (pyruvate carboxylase) | `net_flux` | 12 % |
| `ME` (malic enzyme) | `net_flux` | 8 % |
| `CS` (citrate synthase) | `net_flux` | 10 % |
| `IDH` reductive carboxylation (AKG + CO₂ → Cit) | `reverse_flux` | 20 % |

4. **Measurement fit.** The full map must reproduce both tracer experiments and the measured extracellular rates with a total weighted squared residual of at most 180:

   $ S = \sum_i ((\hat{y}_i-y_i)/\sigma_i)^2 \leq 180 $.

   Here $ y_i $ is an observed raw MID fraction or extracellular rate, $ \sigma_i $ is its supplied `sd`, and $ \hat{y}_i $ is the corresponding prediction from the reported flux map. Include every MID fraction from both experiments and every rate with positive `sd`; glucose uptake is handled by the normalization requirement. Predict raw MID fractions using the reaction atom transitions and the stated tracer purities, natural isotope abundances, symmetry, derivatization, peak truncation, and normalization. For isotope predictions, reconstruct the two one-way rates from the reported net flux and the exchange flux $ x = \min(forward, reverse) $ as $ v_f = x + \max(net,0) $ and $ v_r = x + \max(-net,0) $ for reversible reactions; use $ v_f = \max(net,0) $ and $ v_r = 0 $ for irreversible reactions. This convention resolves independent rounding of the three reported flux columns.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
