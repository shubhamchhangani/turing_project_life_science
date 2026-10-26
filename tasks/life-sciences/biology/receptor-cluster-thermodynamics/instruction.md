# Cooperative membrane-receptor thermodynamics

A membrane patch is represented by a finite lattice of interacting receptors. Each receptor has an inactive state $ s=-1 $ and an active state $ s=+1 $. The coupling map is heterogeneous: positive couplings favor agreement between neighboring receptors, while negative couplings favor opposite states. There is no single-site bias or applied ligand field. This is an idealized equilibrium model of receptor cooperativity, with dimensionless parameters already expressed in units of thermal energy.

Compute the partition function and the full joint activity distribution of the six designated probe receptors. The joint distribution is needed to quantify responses to combinations of locally imposed active or inactive states. Retain the finite lattice and both periodic boundaries in the calculation.

## Inputs and equilibrium model

The inputs are in `/root/data`:

- `model.json` gives `width`, `height`, and six probe records in the ordered list `probes`. Each record has an identifier and integer coordinates `x,y`.
- `couplings.csv` has columns `x,y,direction,K`. A row with direction `E` connects $ (x,y) $ to $ ((x+1)\bmod W,y) $; direction `N` connects $ (x,y) $ to $ (x,(y+1)\bmod H) $. Coordinates start at zero. Every undirected nearest-neighbor edge appears exactly once, including the wraparound edges. The listed decimal values of `K` define the model; there is no additional temperature factor.

For a complete receptor-state assignment $ s $, define

$$ B(s)=\sum_{e=\{u,v\}} K_e s_u s_v, \qquad Z=\sum_{s\in\{-1,+1\}^{WH}} \exp(B(s)). $$

Its equilibrium probability is $ \exp(B(s))/Z $. The six-character probe pattern uses the order of the `probes` list, with `0` meaning inactive and `1` meaning active. Its joint probability is the sum of equilibrium probabilities over all assignments with that probe pattern; all other receptors remain free to fluctuate. This is an unconditional joint distribution, not a product of individual activation probabilities. Use the natural logarithm for $ \log Z $.

## Deliverables

Write these two regular UTF-8 CSV files:

1. `/root/results/joint.csv`, with the literal header `pattern,probability` and exactly 64 data rows, one for every six-character binary pattern. Preserve leading zeros. Rows may be in any order.
2. `/root/results/thermodynamics.csv`, with the literal header `quantity,value` and exactly one data row whose quantity is `log_partition` and whose value is $ \log Z $.

Every numeric field must be a finite ASCII base-10 number, without surrounding whitespace; a decimal point and a base-10 exponent are allowed. Each file must be at most 1 MiB. Include no extra columns or blank rows. LF and CRLF line endings are accepted.

Each joint probability must be within absolute error $ 10^{-7} $ of the finite-model value, and `log_partition` within absolute error $ 10^{-7} $. Probabilities must lie in $ [0,1] $ and sum to one within $ 2\times10^{-7} $. Global state-inversion symmetry requires complementary probe patterns to have equal probabilities within $ 2\times10^{-7} $, and each probe's active-state marginal to equal $ 1/2 $ within $ 2\times10^{-7} $. Retain enough digits in the CSV files to meet these requirements.

Python 3, NumPy, and SciPy are available. You may use any computational method and install additional tools. No particular code structure or algorithm is required.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
