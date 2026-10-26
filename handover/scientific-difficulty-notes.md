# What the difficulty trials mean

A quality-review failure and a difficulty-gate failure are different findings.
The original review found incomplete verification and permissive numeric
parsing. Those defects allowed incorrect output to pass and have concrete
regression tests. The later screenshots instead show models producing
accepted solutions. Fixing the former does not establish the latter.

The receptor candidates use an established idealized biophysical model:
[Duke and Bray, 1999](https://pubmed.ncbi.nlm.nih.gov/10468569/).
The inputs and requested numerical instances are original synthetic data.
They should not be described as unsolved biological discoveries.

Exact inference methods exist for particular Ising graph classes; see
[Schraudolph and Kamenetsky, 2008](https://arxiv.org/abs/0810.4401).
The paper treats planar interactions and does not establish a difficulty
result for this task's toroidal model with three local fields. It does show
why exponentially many receptor configurations alone are not sufficient
evidence that a benchmark will defeat an expert model.

The current inverse task requires reconstructing the entire parameter vector
from collective observables, handling toroidal topology and nonzero local
fields. Its numerical solution and identifiability are independently checked.
Its difficulty remains subject to fresh, isolated trials with the actual
target models and the full declared time allowance.

An interrupted, unauthenticated or quota-limited run is inconclusive. A single
valid failure is evidence about that run, not a guarantee about every future
run or a different Turing harness. No benchmark in this workspace currently
has evidence establishing that all required models fail.
