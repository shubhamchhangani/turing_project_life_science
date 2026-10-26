A fluorescence experiment follows one mother cell and both of its daughters.
Infer the hidden molecular state at division and the transcription history from
photon counts integrated over finite exposure windows. Daughter observations
carry information about the mother, and the daughters share a conserved pool
of molecules at birth.

The file `/root/data/model.json` specifies a finite continuous-time Markov model.
Each cell has state $ (s,m,p) $, where the promoter state $ s $ is 0 or 1,
RNA count $ m $ is an integer from 0 through `rna_capacity`, and protein count
$ p $ is an integer from 0 through `protein_capacity`. These capacities are
part of the model. An event whose destination exceeds a capacity is disabled:
it has zero rate, does not occur and does not contribute to any event count.
The model is a synthetic, finite coarse-grained gene-expression system.

The `rates` object gives the following parameters. All times and rates use one
consistent unit, and the stimulus $ u $ is constant within each exposure:

| Event | Change | Rate |
| --- | --- | --- |
| Promoter activation | $ (0,m,p)\to(1,m,p) $ | $ \text{on}\,u\,[1+\text{feedback}\,p^2/(\text{half}^2+p^2)] $ |
| Promoter deactivation | $ (1,m,p)\to(0,m,p) $ | `off` |
| One-RNA transcription | $ (s,m,p)\to(s,m+1,p) $ | $ \text{transcription}\,g(s)(1-\text{burst_two}) $ |
| Two-RNA transcription | $ (s,m,p)\to(s,m+2,p) $ | $ \text{transcription}\,g(s)\,\text{burst_two} $ |
| RNA degradation | $ (s,m,p)\to(s,m-1,p) $ | $ \text{rna_decay}\,m $ |
| Translation | $ (s,m,p)\to(s,m,p+1) $ | $ \text{translation}\,m $ |
| Protein degradation | $ (s,m,p)\to(s,m,p-1) $ | $ \text{protein_decay}\,p $ |

Here $ g(0)=\text{leak} $ and $ g(1)=1 $. At zero count a degradation rate
is zero. Both transcription reactions count as **one transcription event**,
even though a two-RNA event produces two molecules. These are the only
molecular reactions. There are no events at exposure boundaries other than
the stated stimulus change.

Conditional on a molecular trajectory, photons form an inhomogeneous Poisson
process with instantaneous intensity
$ \text{photon_background}+\text{photon_brightness}\,p(t) $.
Photon emission leaves the molecular state unchanged. All emissions are
detected; there is no censoring, saturation, extra measurement noise, or dead
time. An exposure's `photons` value is its exact count over the entire interval,
not a measurement of its endpoint state. Conditional on the trajectories,
emissions in disjoint intervals and in different cells are independent.

The mother's state starts exactly at `initial = [s,m,p]` with probability one.
The ordered list `mother` specifies consecutive exposures starting at time zero.
Each record gives `duration`, `stimulus`, and observed `photons`. The mother
divides immediately after her last exposure and ceases to exist as a cell.

If the predivision mother has state $ (s,m,p) $, each of her $ m $ RNA molecules
independently goes to daughter A with probability `partition_rna`, and each of
her $ p $ protein molecules independently goes to A with probability
`partition_protein`. These assignments are mutually independent conditional
on the mother. Daughter B receives every molecule not assigned to A. Thus if
A receives $ (a,b) $, the birth states are $ (s,a,b) $ and
$ (s,m-a,p-b) $. Both daughters inherit the mother's promoter state unchanged;
neither promoter is resampled. Division is instantaneous and produces no
photons and no transcription events.

The ordered `daughter_a` and `daughter_b` exposure lists both start at division.
Each daughter subsequently follows the same reaction and photon rules with her
own stimulus schedule. The daughters evolve independently conditional on their
birth states. Their observation durations may differ. Do not extend a daughter
beyond the end of her own list when counting events. There are no further
divisions or observations.

Let $ D $ be all the observed photon counts in all three lists, and let $ N $
be the total number of transcription events in the mother and in both daughters
during their respective observed lifetimes. Write four UTF-8 CSV files under
`/root/results`, with LF or CRLF line endings:

1. `/root/results/evidence.csv`, literal header `metric,value`, with the single row key
   `log_evidence`: the natural logarithm of $ P(D) $ under the full model.
   Include the photon-count probability masses and all partition probabilities.
2. `/root/results/division_posterior.csv`, literal header `state,probability`, with every
   predivision mother state keyed as `s:m:p`, using canonical nonnegative
   integers without leading zeros. The value is $ P(s,m,p\mid D) $ immediately
   before division. Include all states within both capacities.
3. `/root/results/daughter_rna.csv`, literal header `counts,probability`, with every key `a:b`
   for RNA counts in A and B **immediately after division**, with each index
   from 0 through `rna_capacity`. The value is $ P(m_A=a,m_B=b\mid D) $.
   Include impossible pairs with probability zero. Marginalize all promoter
   and protein states and subsequent trajectories.
4. `/root/results/transcription_counts.csv`, literal header `count,probability`, with keys
   `0` through the input `count_limit`, plus the literal key `overflow`.
   Integer key $ k $ means $ P(N=k\mid D) $; `overflow` means
   $ P(N>\text{count_limit}\mid D) $. Do not condition on staying below this
   limit and do not discard the overflow probability.

Rows may appear in any order. Include each required key exactly once, with no
extra rows or columns. All numeric values must be finite ASCII decimal numbers,
optionally in scientific notation, without whitespace, grouping separators or
underscores. Probabilities must lie in [0,1]. The absolute error of every
probability must be at most $ 2\times10^{-8} $ and of `log_evidence` at most
$ 2\times10^{-7} $. Each of the three probability tables must sum to one within
$ 10^{-6} $. Supply probabilities for the full stated model, including hidden
trajectories and all observations.

The working directory is `/root`. Python, NumPy, SciPy and a C/C++ compiler are
available. You may choose any computational method and install additional tools.
Only the four output CSVs need to be retained in `/root/results`.

You have 18000 seconds to complete this task. Do not cheat by using online solutions or hints specific to this task.
