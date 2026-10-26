A synthetic receptor sheet contains 1024 two-state subunits. Neighbor interactions
represent competing conformational preferences: positive couplings favor equal
spins and negative couplings favor opposite spins. Both signs are intentional. Six experimentally controlled
subunits are clamped, and eight other subunits are read out simultaneously.
Compute the complete joint conformational law of the eight reporters under this
intervention, retaining the correlations mediated by the rest of the sheet.

The input `/root/data/model.json` defines the model completely:

- `coordinates`: an array of two-dimensional locations; the zero-based array
  index is the site ID. Locations describe the sheet geometry.
- `edges`: records `[i, j, J]`, each specifying one undirected interaction.
  Each listed interaction occurs once. There are no other edges, periodic
  connections, single-site fields, or unlisted interactions.
- `ligand_sites` and `ligand_spins`: the six site IDs and their imposed spins.
- `reporter_sites`: the eight site IDs, in the required output order.

Every subunit has spin $ s_i \in \{-1,+1\} $. Couplings already include the
inverse temperature. For a complete spin assignment the dimensionless weight is
$ w(s)=\exp(\sum_{(i,j,J)\in E} J s_i s_j) $.
The unconstrained partition function is $ Z=\sum_s w(s) $, with the sum over
all assignments at all 1024 sites. This is a counting sum, not a uniform-prior
average. Let $ C $ denote the event that every ligand site has its specified
spin. Then $ P(C)=\sum_{s\in C} w(s)/Z $. The requested reporter law is
$ P(s_R=r\mid C) $. Clamping conditions on the stated spins; it does not remove
the sites or their incident interactions. This is a finite synthetic equilibrium
model, not a molecular prediction for a specific natural receptor.

Write these two UTF-8 CSV files, with LF or CRLF line endings:

1. `/root/results/thermodynamics.csv`, with literal header `metric,value` and
   exactly two rows, keyed by `log_partition` and `log_clamp_probability`.
   Use natural logarithms of $ Z $ and $ P(C) $, respectively.
2. `/root/results/reporter_distribution.csv`, with literal header
   `pattern,probability` and exactly 256 rows, one for each eight-character
   binary string. Character position $ j $ refers to `reporter_sites[j]`;
   `0` means spin +1 and `1` means spin -1. For example, `10000000` specifies
   only the first reporter in state -1. Each probability is conditional on $ C $.

Rows may appear in any order. Include each key exactly once and no extra
columns or rows. Numeric fields must be finite ASCII decimal numbers, optionally
in scientific notation, without whitespace, grouping separators, or underscores.
Probabilities must lie in [0,1]. Supply every probability with absolute error at
most $ 5\times10^{-8} $, their sum within $ 10^{-6} $ of 1, and each logarithm
with absolute error at most $ 2\times10^{-6} $. These accuracy requirements concern
the full finite model, including unobserved spins and all listed edges.

The working directory is `/root`. Python, NumPy, SciPy and a C/C++ compiler are
available. You may use any computational method and install additional tools.
Only the two output CSVs need to be retained in `/root/results`.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
