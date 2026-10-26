# turing/metabolic-flux-analysis — version 2.0.0 candidate

Infer a common central-carbon flux map and four metabolite pool sizes from two
parallel isotope pulse-chase experiments. The 31-reaction network has 10 free
net fluxes, eight reversible exchanges, and four unknown pool sizes. Observations
include 1,924 raw GC-MS MID fractions at 13 sampling times per experiment, plus
eight noisy extracellular rates. Glucose uptake is fixed independently.

This version extends the steady-state task to isotopically nonstationary MFA.
The scientific objective requires reconstructing carbon transitions, dynamic
pool turnover, pulse and chase boundary conditions, rotational symmetry,
tracer purity, derivative isotope contributions, and a joint nonlinear fit.
The extra complexity comes from physical modeling rather than hidden formatting
rules or arbitrary accuracy tightening. Version 1.1.0 remains available
separately; the user's external evaluation showed both listed frontier models
solved that version. No frontier-model failure is claimed for version 2.0.0.

## Scientific basis

Metabolism is stationary but isotope fractions change. Intracellular pools are
well mixed and constant, sources follow the supplied piecewise constant
labeling schedule, and consumption removes isotopomers proportionally to their
pool fractions. Products follow carbon-atom transfers and molecule symmetry.
The two experiments reset independently to natural abundance. This is a
synthetic experimental design, not a claimed reproduction of a biological study.

Isotopically nonstationary MFA jointly fits fluxes and pool-size parameters from
labeling dynamics. The model class is described in [Young (2014), INCA](https://pmc.ncbi.nlm.nih.gov/articles/PMC3998137/)
and [INCA 2.0](https://pmc.ncbi.nlm.nih.gov/articles/PMC8789327/).
The inherited carbon-transition and derivatization models were previously
cross-checked with FreeFlux and IsoCor for the original task. New transient
checks are recorded separately below; those earlier checks are not represented
as a validation of new time-course data.

## Independent expected values and data generation

The author-only `tools/metabolic-flux-analysis-v2/make_dataset.py` generates data
with seed 20261023. It injects known fluxes, not fitted reference output:
G6PDH net 15, PDH net 42, PC net 11, ME net 40, CS net 30, and IDH reverse 8.
Unknown pool sizes are OAA 8, Cit 45, AKG 22, and Mal 30 nmol per mg protein.
The complete source flux vector remains outside the task package; only the
stated accuracy targets are placed in verifier-only expected data.

The glucose pulse ends at 4 minutes; the glutamine pulse ends at 6.5 minutes.
Sampling covers early turnover, later labeling, and washout. MID noise is
independent Gaussian with SD 0.002 per reported peak; extracellular rate noise
has 3% relative SD. Fractions are neither clipped nor renormalized after noise.
All input measurements are independent of the numerical reference fit.

## Reference solution

`solution/dynamic_fit.py` derives a stoichiometric null-space flux
parameterization, solves the transient isotopomer balances, and fits the
22 parameters jointly from four deterministic random starting points.
It reads only supplied measurements, known pools, and tracer schedules; it does
not read injected true values. `network_model.py` contains the biochemical atom
mapping and stoichiometric parameterization. `dynamic_model.py` evaluates the
ODE balances, with a compiled sparse product-formation kernel; `forward_model.py`
contains the carbon and GC-MS isotope chemistry.

Numerical normalization of substrate probability vectors suppresses round-off
drift in repeated bimolecular convolution. It leaves exact normalized
isotopomer balances unchanged. No numerical integration method or fitting
algorithm is required of submissions.

## Verification and tolerance rationale

The separate verifier reads only two submitted data artifacts, `fluxes.csv`
and `pools.csv`. It never executes submitted code. It first checks exact headers,
row sets, ASCII decimal syntax, finite values, positive pools, and nonnegative
one-way fluxes. Both artifact reads use nonblocking descriptors, a regular-file
check, and a 1 MiB limit, preventing FIFO hangs and unbounded reads.

The 0.05 directional and 0.1 balance allowances accommodate rounding. Reporting
three flux columns to two decimals introduces at most 0.015 directional error.
The largest absolute stoichiometric row sum is six (T3P), giving a maximum
balance-rounding contribution of 0.03. These allowances do not permit negative
one-way fluxes. Glucose normalization has an independent absolute tolerance
of 0.0001. The prediction convention for rounded fluxes is explicit in the
instructions.

The six flux accuracy tolerances are stated in the instructions and remain
conservative; each of the four unknown pools must be within 25% of its injected
value. The noise sensitivity checks below assess these tolerances on new
independent realizations; they do not claim uniform coverage over every possible
noise realization.

The full map must also achieve weighted measurement SSR at most 2110. There are
1,932 noisy observations and 22 fitted parameters, giving a nominal 1,910
residual degrees of freedom. The chi-square 99.9th percentile is 2106.7084,
rounded upward to 2110. This threshold follows the declared Gaussian noise
model, not a stored reference fit. A unique global optimum is not required;
statistically equivalent maps are allowed.

Full-data predictions use trusted model equations and input measurements, not a
fitted oracle vector. The model runs in a subprocess with a 60-second timeout;
the entire verifier has a 100-second deadline plus a 5-second kill grace, inside
its declared 120-second allocation. NumPy, SciPy, and Numba are preinstalled
at pinned versions. Reference execution and verifier checks use 1 CPU and 2 GiB.

## Validation

The complete reference solution was run offline in the built agent image,
limited to 1 CPU and 2 GiB RAM. Four starts completed in 120.5 seconds; the best
weighted SSR was 1961.013, below 2110. Its two artifacts passed all 16 tests
in the separate no-network verifier with reward 1; verification took 0.60 seconds. The solution inferred OAA 9.442, Cit
43.971, AKG 21.330, and Mal 28.137, all within the declared truth tolerances.

Six additional noise realizations (seeds 1000–1005) were generated and fitted
from two starts each. All met the target tolerances. The largest observed
relative errors, compared with the declared tolerances, were:

| Quantity | Largest error across six seeds | Tolerance |
|---|---:|---:|
| G6PDH net | 11.34% | 45% |
| PDH net | 2.40% | 12% |
| PC net | 2.05% | 12% |
| ME net | 1.10% | 8% |
| CS net | 1.17% | 10% |
| IDH reverse | 4.86% | 20% |
| OAA pool | 16.60% | 25% |
| Cit pool | 1.23% | 25% |
| AKG pool | 4.13% | 25% |
| Mal pool | 3.71% | 25% |

Their fitted SSR values ranged from 1823.08 to 1910.60. These checks support
solvability and tolerance selection; they are not measurements of frontier
model difficulty. The slightly broader OAA error in the submitted realization
(18.02%) also remains within the stated 25% tolerance.

An independently written vector-form RHS using RK45 agreed with the compiled
RHS using DOP853 to within 6.27e-10 in predicted raw MID fractions. This checks
the new integration and indexing implementation; both implementations share
the previously validated biochemical reaction specification.

At the injected true fluxes and pool sizes, these deliberately incorrect
forward predictions have MID-only SSR far above the full-data limit of 2110:

| Incorrect prediction | MID-only SSR |
|---|---:|
| Treat every sample as isotopic steady state during the pulse | 8,591,171 |
| Ignore the chase and leave tracer on | 1,260,658 |
| Set all four unknown pools to 1 | 378,748 |

These are forward-model negative controls, not claims that every optimized
misspecified model or every frontier agent must fail.

The author-only regression suite tests malformed headers, negative directional
fluxes, missing and duplicate rows, NaN, missing/empty/oversized files, FIFOs,
whitespace and underscore numeric fields, valid decimal exponents, and pool
header, row, numeric, positivity, accuracy, and FIFO defects. The injected true
map passes; every tested invalid artifact is rejected.

Authoring and diagnostic scripts stay outside the uploaded task package. The
agent image contains only the seven input files under `/root/data`; solution
code, expected values, and verifier code are isolated from it.

## Frontier evaluation status — solved, do not submit

A fresh GPT-6-sol max-effort trial using `npx @openai/codex@0.164.0` and the
existing ChatGPT login solved this candidate. The model received only the
instructions and agent-visible inputs in an isolated Docker environment.
Agent execution took 11 minutes 18 seconds. Its artifacts passed all 16 tests
and received reward 1, with a reported weighted SSR of 1961.0084. This candidate
does not meet the user's required difficulty; no new ZIP was created after
the result. The previous candidate ZIP is not recommended for submission.

Original trial logs and artifacts are retained under
`jobs/mfa-v2-frontier/gpt6-sol-max-v2-npx-chatgpt`. This is a local Codex-harness
result, not an assertion about an unrun Opus or Terminus trial.

The earlier API-key attempt under `jobs/mfa-v2-frontier/gpt6-sol-max-v2` ended
with `ApiUsageLimitError` before a substantive attempt. Its recorded reward 0
remains classified as an infrastructure error, not evidence of model failure.
