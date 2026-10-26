Recover the central-carbon metabolic flux map and four intracellular metabolite pool sizes for each of two cell populations, `A` and `B`, from four parallel ¹³C pulse-chase experiments. GC-MS measures the pooled cell extract. The two populations have different metabolic phenotypes; every experiment contains the same populations in the same proportions.

## Biological and measurement model

Each population is homogeneous and perfectly mixed. Its fluxes and metabolite amounts are constant, but its isotope fractions change throughout an experiment. Fluxes and pool sizes may differ between populations and are shared across the four experiments. Each experiment starts independently with every intracellular pool at the supplied natural isotope abundance. Population protein fractions remain constant. There is no metabolite exchange between populations, unlabeled dilution, inactive pool, compartmentation within a population, or kinetic isotope effect.

The external medium is continuously clamped to the supplied isotope fractions. Delivery and mixing are instantaneous. Substrate concentrations and uptake rates do not change at tracer switches. Species ending in `.ext` and `CO2` are external sources or sinks. Released CO₂ does not re-enter metabolism; carboxylation uses the continuously clamped external CO₂ pool. Intracellular labeling evolves by carbon-atom balances and pool turnover. Consumption removes isotopomers in proportion to their fractions in the consumed pool.

`T3P` lumps DHAP and GAP, whose isotope distributions equilibrate within each population. `P5P` lumps Ru5P, R5P, and X5P. `PG3` is 3-phosphoglycerate. Succinate and fumarate are rotationally symmetric; no other metabolite is symmetrized. Carbon-atom transfers follow the standard biochemistry of the supplied reactions.

Fluxes and pool sizes are expressed per mg protein **of the corresponding population**. Each population has glucose uptake exactly 100 nmol per minute per mg of its own protein. The protein fractions in `populations.csv` are exact. Bulk extracellular rates are normalized per mg of total culture protein:

$ v_{bulk,r} = w_A v_{A,r} + w_B v_{B,r} $.

Extraction recovers every metabolite equally from both populations, and the same metabolite has the same instrument response in both. For metabolite $ m $, let $ p_{c,m}(t) $ be the full carbon-isotopomer probability vector in population $ c $, and $ P_{c,m} $ its pool size. The extracted probability vector is

$ p_{bulk,m}(t) = (w_A P_{A,m}p_{A,m}(t) + w_B P_{B,m}p_{B,m}(t))/(w_A P_{A,m}+w_B P_{B,m}) $.

Use this molecular mixture to calculate each fragment's backbone mass distribution. Convolve it with natural isotopes of the remaining ion atoms, retain the specified peaks, and then normalize their sum. Pool amounts therefore affect both labeling turnover and each population's contribution to the measured ion signal. The supplied population labels are anchored by their fixed pool sizes.

## Inputs in `/root/data/`

- `network.tsv`: tab-separated columns `reaction_id`, `equation`, `reversible`. Both populations use this same complete 31-reaction network. Balance every intracellular metabolite separately in each population; external species and CO₂ are not balanced.
- `populations.csv`: columns `population`, `protein_fraction`, giving the fixed fractions of total culture protein contributed by `A` and `B`.
- `pool_sizes.csv`: columns `population`, `metabolite`, `pool_size`. Known values are exact nmol per mg protein of that population. `OAA`, `Cit`, `AKG`, and `Mal` are marked `unknown` in each population, giving eight unknown pool sizes.
- `labeling_schedule.csv`: columns `experiment`, `start_min`, `co2_13c`, `glucose_C1` through `glucose_C6`, and `glutamine_C1` through `glutamine_C5`. Each row applies from its start time until the next row for that experiment, or indefinitely for its last row. The four experiments are `glc12`, `gln5`, `glcU`, and `dual`; the table defines every switch and purity. Fractions are absolute ¹³C probabilities, including natural abundance. Carbon labeling is independent across positions within each substrate molecule. The same external labeling applies to both populations.
- `fragments.csv`: columns `fragment`, `metabolite`, `backbone_carbons`, `ion_formula`. Backbone carbon positions are separated by semicolons and use IUPAC numbering; C1 is the carboxyl carbon of amino and carboxylic acids. The formula includes the entire measured TBDMS ion, including backbone and derivatization atoms.
- `isotope_abundances.csv`: authoritative natural isotope abundances, with columns `element`, `mass_shift`, `abundance`, for C, H, N, O, and Si.
- `mid_measurements.csv`: columns `experiment`, `time_min`, `fragment`, `isotopologue`, `fraction`, `sd`. These are raw pooled mass-isotopomer fractions, including natural isotopes. Times are minutes after the initial switch. For $ n $ backbone carbons, peaks $ M_0 $ through $ M_{n+2} $ are normalized by their sum before independent additive Gaussian noise is applied. No clipping or renormalization follows noise addition, so small negative measured fractions are possible.
- `extracellular_rates.csv`: columns `reaction_id`, `rate`, `sd`. These are bulk net rates in nmol per minute per mg total culture protein. Positive-SD rates have independent additive Gaussian noise. Glucose normalization is exact.

## Outputs

Create two regular UTF-8 comma-separated files, each no larger than 1 MiB.

`/root/results/fluxes.csv` must have the literal header `population,reaction_id,net_flux,forward_flux,reverse_flux` and exactly 62 data rows: every combination of `A` or `B` with the 31 reaction identifiers, once each, in any order. Net flux is in the written reaction direction. Forward and reverse fluxes are nonnegative one-way rates. All rates use the population-specific protein normalization above.

`/root/results/pools.csv` must have the literal header `population,metabolite,pool_size` and exactly eight data rows: every combination of `A` or `B` with `OAA`, `Cit`, `AKG`, and `Mal`, once each, in any order. Report only those unknown pools, in nmol per mg protein of the named population.

All numeric fields must be finite ASCII decimals without whitespace. Optional signs, decimal points, and base-10 exponents are allowed, including `100`, `-0.05`, `.5`, and `1e2`. Underscores and other numeric formats are not allowed. Pool sizes must be strictly positive.

## Required accuracy, consistency, and measurement agreement

1. Both files satisfy their exact headers, row sets, size limits, and decimal syntax.
2. In each population, `GLCup` net flux is 100 within absolute tolerance 0.0001. Both one-way rates are nonnegative and $ |forward-reverse-net| \leq 0.05 $. Irreversible reactions have reverse flux within 0.05 of zero and net flux at least −0.05. The net-flux balance residual of each intracellular metabolite is at most 0.1 in absolute value.
3. For **each population separately**, the following six quantities must satisfy $ |reported-true| \leq tolerance\times|true| $, where $ true $ is the value used to generate the measurements:

| Reaction quantity | Column | Relative tolerance |
|---|---|---|
| `G6PDH` oxidative PPP | `net_flux` | 45% |
| `PDH` | `net_flux` | 12% |
| `PC` | `net_flux` | 12% |
| `ME` | `net_flux` | 8% |
| `CS` | `net_flux` | 10% |
| `IDH` reductive carboxylation | `reverse_flux` | 20% |

4. Each of the eight inferred pools must be within 25% of its own population-specific value used to generate the measurements.
5. The full pair of flux maps and pool vectors must reproduce all MID fractions and all positive-SD bulk rates with weighted squared residual

   $ S=\sum_i((\hat y_i-y_i)/\sigma_i)^2\leq4090 $.

   Use every supplied observation, its SD, the tracer schedules, isotope abundances, symmetry, fragment chemistry, and molecular mixing rule. For predictions, reconstruct reversible one-way rates from the submitted net rate and $ x=\min(forward,reverse) $ as $ v_f=x+\max(net,0) $ and $ v_r=x+\max(-net,0) $. For irreversible reactions use $ v_f=\max(net,0) $ and $ v_r=0 $. Isotopomer loss from an intracellular pool uses its total incoming one-way material flux, equal to consumption for an exactly balanced map. These conventions resolve independent rounding while preserving normalized isotopomer fractions.

Any flux maps and pool vectors satisfying these conditions are acceptable. A unique global optimum or uncertainty estimates are not required.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
