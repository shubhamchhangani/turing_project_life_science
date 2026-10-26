# Recover membrane-receptor cooperativity

A synthetic membrane patch contains two-state receptors on a rectangular lattice
with periodic boundaries in both directions. Receptor state $ s_{x,y} $ is +1
(active) or -1 (inactive). Twelve specified receptors experience unknown local
ligand biases; every other receptor has zero local bias. Recover every interaction
and the twelve biases from the supplied equilibrium observations.

The finite-system equilibrium distribution is

$$ P(s)=Z^{-1}\exp\left(\sum_{e=(u,v)} K_e s_u s_v + \sum_{a=1}^{12} h_a s_{r_a}\right). $$

All parameters are dimensionless: thermal energy has already been absorbed into
$ K $ and $ h $. Each undirected edge occurs once in the sum. An edge labelled
`E` at `(x,y)` joins `(x,y)` to `((x+1) mod width,y)`; an edge labelled `N`
joins `(x,y)` to `(x,(y+1) mod height)`. Positive interaction favours equal
states; negative interaction favours opposite states. Positive bias favours
the active state. A biased receptor remains a fluctuating variable, not a
clamped boundary. There are no additional interactions or adjustable parameters.

## Inputs

All inputs are in `/root/data/`:

- `model.json`: lattice dimensions, coordinates and IDs of the twelve biased
  receptors, and inclusive bounds for every coupling and bias.
- `edge_moments.csv`: one row per edge, with columns
  `x,y,direction,mean_product`. The observation is the raw equilibrium moment
  $ \langle s_u s_v\rangle $, not a connected covariance.
- `site_moments.csv`: columns `receptor_id,mean_state`. The observation is
  $ \langle s_{r_a}\rangle $, not the probability that the receptor is active.

The observations all describe the same finite-system equilibrium distribution.
They were rounded to twelve decimal places; there is no additional sampling
noise. The generating parameters lie strictly inside the stated bounds. Infer
the unique parameter vector consistent with these moments to the accuracy below.

## Required outputs

Write two UTF-8 regular files, each at most 1 MiB:

- `/root/results/couplings.csv`, with the literal header `x,y,direction,K`,
  followed by exactly one row for each of the 646 input edges.
- `/root/results/biases.csv`, with the literal header `receptor_id,h`, followed
  by exactly one row for each of the twelve input receptor IDs.

Row order is arbitrary. Use the input edge labels and receptor IDs exactly.
Coordinates must be nonnegative ASCII integers without leading zeroes or signs.
Parameter fields must be finite ASCII decimal numbers: an optional sign, a
decimal mantissa with at least one digit, and an optional `e` or `E` exponent
with an optional sign and at least one digit. Whitespace, underscores, NaN and
infinity are not allowed in numeric fields. CSV quoting is allowed in data rows.
Use LF or CRLF line endings; do not add columns, blank rows, or extra rows.
All parameters must respect their inclusive bounds in `model.json`.

Every inferred coupling and bias must be within $ 5\times10^{-5} $ of its
generating value. In addition, the finite-system equilibrium distribution
defined by your entire parameter vector must reproduce every supplied edge
and site moment within absolute error $ 2\times10^{-6} $. Report enough decimal
places to satisfy both requirements. Python, NumPy, SciPy and PFAPACK are installed;
you may implement or install other numerical tools.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
