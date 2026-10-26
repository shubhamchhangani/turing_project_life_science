Infer the central-carbon flux map and four intracellular metabolite pool sizes from two parallel ¹³C pulse-chase experiments in a proliferating cell culture. Metabolism and metabolite amounts remain at steady state, but isotope labeling changes throughout each experiment. Both cultures share the same constant fluxes and pool sizes. Intracellular metabolites were measured by GC-MS as TBDMS derivatives.

## Pulse-chase experiments

Before time zero, every intracellular pool and extracellular substrate has natural isotope abundance. At time zero, `glc12` receives [1,2-¹³C₂]glucose and natural-abundance glutamine; `gln5` receives [U-¹³C₅]glutamine and natural-abundance glucose. The labeled substrate is subsequently replaced with natural-abundance substrate. Each experiment starts independently from the same natural-abundance initial state.

`/root/data/labeling_schedule.csv` is authoritative for every switch time and per-carbon ¹³C fraction. Each row applies from its `start_min` until the next row for that experiment, or indefinitely for its last row. Carbon labeling is independent across positions within each substrate molecule. Substrate concentrations, uptake rates, and metabolic fluxes do not change at a switch. Delivery and extracellular mixing are instantaneous; there is no transport delay or instrument response lag. CO₂ is a continuously clamped external pool with the listed natural ¹³C fraction.

Intracellular pools are homogeneous, perfectly mixed, constant in size, and shared between the experiments. Their isotopomer fractions evolve according to carbon-atom balances with pool turnover; they are not instantaneously equilibrated to the new tracer. Label is lost through all consumption reactions at the labeling distribution of the consumed pool. There is no additional unlabeled dilution, inactive pool, or compartment. Succinate and fumarate are rotationally symmetric; no other metabolite is symmetrized.

## Input files (all in `/root/data/`)

- `network.tsv` (tab-separated; columns `reaction_id`, `equation`, `reversible`): the complete metabolic network of 31 reactions. Use exactly these reactions, metabolites and reversibilities.
  - Species ending in `.ext` are extracellular sources or sinks. `CO2` is an external pool and is not balanced.
  - Carboxylating reactions (`PC`, and `IDH` in its reverse direction) fix CO₂ with natural-abundance labeling. CO₂ released inside the cell does not re-enter.
  - `T3P` is the lumped triose-phosphate pool (DHAP + GAP); its components equilibrate isotopically with each other while the labeling of the combined pool changes over time. `P5P` is the lumped pentose-5-phosphate pool (Ru5P, R5P, X5P). `PG3` is 3-phosphoglycerate.
  - The model has a single compartment and no kinetic isotope effects. Succinate (`Suc`) and fumarate (`Fum`) are rotationally symmetric molecules. No other metabolite is treated as symmetric.
  - Carbon atom transitions are not given; they follow the standard biochemistry of each reaction.
- `fragments.csv` (columns `fragment`, `metabolite`, `backbone_carbons`, `ion_formula`): one row per measured GC-MS fragment ion.
  - `backbone_carbons` lists, separated by `;`, which carbons of the metabolite are contained in the ion. Carbons are numbered by IUPAC convention, so C1 is the carboxyl carbon of the amino and carboxylic acids.
  - `ion_formula` is the full elemental formula of the measured ion. It includes the metabolite backbone and all derivatization atoms.
- `isotope_abundances.csv` (columns `element`, `mass_shift`, `abundance`): natural isotope abundances to use for C, H, N, O and Si. These values are authoritative.
- `mid_measurements.csv` (columns `experiment`, `time_min`, `fragment`, `isotopologue`, `fraction`, `sd`): raw time-course mass isotopomer distributions (MIDs), **not corrected for natural isotopes**. Times are minutes after the initial tracer switch.
  - For a fragment with $ n $ backbone carbons, isotopologues $ M_0 $ through $ M_{n+2} $ are reported as fractions of the sum over those $ n+3 $ peaks.
  - Each fraction carries additive Gaussian measurement noise with the given standard deviation. Values were not clipped or renormalized after noise was added; small negative measured fractions are possible.
- `extracellular_rates.csv` (columns `reaction_id`, `rate`, `sd`): measured net rates of uptake and secretion reactions, with standard deviations.
  - Flux units are nmol per minute per mg protein; `GLCup` is 100 exactly.
- `labeling_schedule.csv`: `experiment`, `start_min`, `co2_13c`, `glucose_C1` through `glucose_C6`, and `glutamine_C1` through `glutamine_C5`. Fractions are absolute ¹³C probabilities, including natural abundance.
- `pool_sizes.csv` (columns `metabolite`, `pool_size`): fixed intracellular pool sizes in nmol per mg protein, except `OAA`, `Cit`, `AKG`, and `Mal`, whose values are marked `unknown` and must be inferred. Known entries are exact. External substrate and CO₂ pools are clamped sources and have no intracellular pool-size parameter.

## Output — `/root/results/fluxes.csv`

A regular UTF-8 comma-separated file no larger than 1 MiB, with header `reaction_id,net_flux,forward_flux,reverse_flux` and exactly 31 data rows, one for every reaction in `network.tsv`, in any order.

- All flux values are finite ASCII decimal numbers in the same flux units (glucose uptake = 100), with no whitespace inside a field. An optional sign, decimal point, and base-10 exponent are allowed (for example `100`, `-0.05`, `.5`, or `1e2`); underscores and other numeric formats are not allowed.
- `net_flux` is the net steady-state flux in the direction the equation is written.
- `forward_flux` and `reverse_flux` are the one-way fluxes, both nonnegative, with $ forward\_flux - reverse\_flux = net\_flux $.
- For irreversible reactions, `reverse_flux` is zero within the 0.05 consistency tolerance below.
- Report a flux map satisfying the consistency, accuracy, and measurement-fit requirements below. Uncertainty estimates are not required. Flux maps that are indistinguishable within the stated tolerances are acceptable; a unique global optimum is not required.

## Additional output — `/root/results/pools.csv`

A regular UTF-8 CSV file no larger than 1 MiB with literal header `metabolite,pool_size` and exactly four data rows: `OAA`, `Cit`, `AKG`, and `Mal`, each once, in any order. Pool sizes are in nmol per mg protein and must be positive finite numbers using the same decimal syntax as fluxes. Report only these four inferred pools; use the fixed amounts from the input for every other intracellular metabolite.

## Output requirements

The output is valid only when all of the following requirements hold.

1. **Contract.** Both output files have their specified literal headers, row sets, decimal syntax, and finite values. Fluxes include exactly one row for each of the 31 reactions; pools include exactly one row for each of the four unknown metabolites.
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

4. **Pool accuracy.** Each of the four inferred pool sizes must be within 25 % of its value used to produce the data: $ |reported - true| \leq 0.25 \times true $.

5. **Measurement fit.** The full flux map and pool sizes must reproduce all sampled times in both pulse-chase experiments and the measured extracellular rates with a total weighted squared residual of at most 2110:

   $ S = \sum_i ((\hat{y}_i-y_i)/\sigma_i)^2 \leq 2110 $.

   Here $ y_i $ is an observed raw MID fraction or extracellular rate, $ \sigma_i $ is its supplied `sd`, and $ \hat{y}_i $ is the corresponding time-dependent prediction from the reported flux map and pool sizes. Include every MID fraction from both experiments and every rate with positive `sd`; glucose uptake is handled by the normalization requirement. Predict raw MID fractions using the reaction atom transitions and the stated tracer purities, natural isotope abundances, symmetry, derivatization, peak truncation, and normalization. For isotope predictions, reconstruct the two one-way rates from the reported net flux and the exchange flux $ x = \min(forward, reverse) $ as $ v_f = x + \max(net,0) $ and $ v_r = x + \max(-net,0) $ for reversible reactions; use $ v_f = \max(net,0) $ and $ v_r = 0 $ for irreversible reactions. For isotopomer loss in each intracellular pool, use its total incoming one-way material flux, which equals its total consumption flux for an exactly balanced map. These conventions resolve independent rounding of the reported flux columns while preserving normalized isotopomer fractions.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
