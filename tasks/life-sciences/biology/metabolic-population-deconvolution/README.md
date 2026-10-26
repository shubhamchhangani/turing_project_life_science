# Population-resolved metabolic flux inference — solved, do not submit

This synthetic task asks for two population-specific flux maps and eight pool
sizes from pooled transient GC-MS measurements. Four tracer programs share
the same biological parameters. There are 44 unknowns, 3,848 noisy MID fractions,
and eight noisy extracellular rates. Exact glucose uptake fixes the scale.

This directory has not been zipped. A fresh GPT-6-sol max-effort npx trial
solved it and received reward 1, so it does not meet the user's required
difficulty. The previous single-population transient candidate was also solved
and held from further packaging. Both remain available for inspection only.

## Scientific model and provenance

Both populations use the specified central-carbon reaction network. They are
noninteracting under the stipulated clamped-medium conditions. Known protein
fractions and different fixed pool sizes anchor their identities. Bulk labeling
is weighted by each population's contribution of metabolite molecules; bulk
rates use protein-fraction weights. Derivative isotope convolution and peak
normalization follow molecular mixing. This distinction is explicit in the
instructions, together with all boundary, purity, symmetry, and noise rules.

Mixed-culture flux inference is supported by
[Gebreselassie and Antoniewicz (2015)](https://pubmed.ncbi.nlm.nih.gov/26219674/).
Pool-dependent deconvolution has also been used for compartmentalized metabolism
by [Stern et al. (2023)](https://www.nature.com/articles/s41467-023-42824-z).
The present dataset is an original synthetic experiment, not a reproduction of
either study or a claim to have solved an open scientific question.

The author-only generator `tools/metabolic-mixture/prototype.py generate`
uses seed 20261024 and injects fluxes independently of the reference fit.
Population A has G6PDH 15, PDH 42, PC 11, ME 40, CS 30, IDH reverse 8; population
B has 36, 30, 28, 65, 22, and 24 respectively. Unknown pools OAA/Cit/AKG/Mal are
8/45/22/30 in A and 20/20/42/14 in B. MID SD is 0.001; noisy bulk rates have
3% relative SD. Neither the generator nor the full injected vector is installed
in the agent image.

## Solution and verifier

The reference derives steady-state fluxes from a stoichiometric parameterization,
integrates the population-specific isotopomer balances, combines molecular ion
signals, and fits both populations jointly from independent starting points.
It reads only supplied inputs. An exact-key cache reuses an unchanged
population's dynamics during finite differences; it never rounds parameters.

The verifier reads only the two declared regular files with nonblocking opens
and a 1 MiB size bound. It checks both contracts before scientific computation:
literal headers, unique population/reaction or population/metabolite pairs,
ASCII decimal syntax, finite values, positive pools, and nonnegative one-way
rates. It then checks both populations' balances, all 12 injected flux targets,
all eight pools, and full-data measurement fit. It does not execute agent code.
Expected targets come from injected parameters, not stored reference output.

The 0.05 directional tolerance covers up to 0.015 error from independently
rounding three columns to two decimals. The 0.1 balance tolerance exceeds the
maximum 0.03 contribution from the largest stoichiometric row sum of six.
Glucose has its separate absolute tolerance of 0.0001. The convention for
rounding in isotope predictions is explicit in the instructions.

The fit limit 4090 is the upward-rounded 99.9th percentile, 4087.5345, of a
nominal chi-square distribution with 3856−44 = 3812 residual degrees of freedom.
It follows the Gaussian measurement model rather than a fitted reference score.
The flux and pool tolerances require separate noise-sensitivity validation;
they are not a guarantee of coverage for every possible noise realization.

The predictor has a 60-second subprocess timeout. The whole test suite has a
100-second timeout with five seconds of kill grace, within a 120-second verifier
allocation. The agent and verifier have pinned scientific dependencies; the
verifier runs without network access in a separate environment.

## Validation status

The first independently initialized author fit reached SSR 3897.2860 and met
all 12 flux and eight pool tolerances. Its Jacobian had 44 nonzero singular
values, with minimum 0.0989 and condition number 762.6. This supports local
identifiability at that solution, not a proof of global uniqueness.

An independently written pooled-ion-count calculation with RK45 integration
agreed with the production observation model to 5.11e-8 per raw MID fraction.
Both share the biochemical reaction specification. Pool normalization error
was below 8e-16. Omitting pool-amount weights at the injected parameters gives
MID-only SSR 90074.7; this is a forward negative control, not a refitted
misspecified model or evidence of frontier difficulty.

The complete packaged reference ran without network access on one CPU and
2 GiB RAM. Two independent starts completed in 373.69 seconds, reaching best
SSR 3896.9940. Its artifacts passed all 19 tests in the separate offline
verifier with reward 1; verification took 0.82 seconds.

Four additional noise realizations, seeds 1000–1003, met every target tolerance.
Their SSR values ranged from 3708.93 to 3853.85. These sensitivity fits started
at the independently fitted estimate from the original dataset, not at injected
parameters. They assess local noise sensitivity, not global optimization
robustness. The largest flux-target relative error was 11.70% (G6PDH, allowed
45%); the largest PDH error was 11.65% (allowed 12%). The largest pool error
was 16.94% (AKG, allowed 25%).

The behavioral regression suite accepted four valid artifact cases and rejected
28 invalid cases. These include the previous PGI exchange loophole, malformed
headers and numbers, missing or duplicated population pairs, nonregular and
oversized files, and opposing population PPP shifts that preserve all balances,
bulk rates, and target tolerances but violate full isotope-data agreement.
All 40 static package checks passed. These checks cover known failure modes;
they do not guarantee an external review outcome.

The fresh GPT-6-sol max trial used `npx @openai/codex@0.164.0` and the existing
ChatGPT login in an isolated Harbor container. It saw only the instructions
and input files. Agent execution took 613.69 seconds. Its 62 flux rows and
eight pools passed all 19 tests, receiving reward 1; the agent reported
weighted SSR 3896.9618, below 4090. Original records are retained under
`jobs/mfa-mixture-frontier/gpt6-sol-max-mixture-npx-chatgpt`.

No Opus trial was launched: the successful GPT trial already disproves the
all-models-fail condition for this candidate. No submission ZIP was created.
