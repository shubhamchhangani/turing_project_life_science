# Experimental candidate: population-resolved metabolic flux inference

Status: author prototype under validation; no task ZIP, no frontier failure claim.

The previous single-population transient task was solved by GPT-6-sol (max)
in an isolated npx trial. This experiment adds a distinct scientific inference
problem: recover two population-specific metabolic maps from bulk measurements.

Two noninteracting populations share a clamped external medium. Their protein
fractions are known. Each has its own constant fluxes, reversible exchanges,
and intracellular metabolite amounts. Four independent pulse-chase tracer
programs probe the same populations. Measurements pool both populations.

For metabolite m, the isotopomer distribution in the extract is the sum of
each population's distribution weighted by its contribution of m molecules:
protein_fraction × metabolite_pool_per_population_protein. Normalize that sum
by the total amount of m. Apply derivatization isotope convolution and peak
truncation after mixing. Averaging the separately normalized raw MID vectors
does not implement this observation process. Extracellular rates instead use
protein-fraction weighting. All conventions will be explicit in instructions.

The draft has 44 unknown parameters, 62 reported reaction rows, eight inferred
pool sizes, and 3,848 MID measurements plus eight noisy extracellular rates.
Injected parameters generate measurements; the reference fitter reads only
input files and starts independently of the injected vector. Different known
pool sizes anchor population identities. Local Jacobian rank, multiple starts,
noise realizations, independent forward checks, and negative controls must be
reviewed before treating this as a scientifically valid task.

Population-specific flux inference from mixed cultures is an established
research problem: [Gebreselassie and Antoniewicz (2015)](https://pubmed.ncbi.nlm.nih.gov/26219674/).
Related compartment deconvolution combines pool-size-weighted labeling with
metabolic inference: [Stern et al. (2023)](https://www.nature.com/articles/s41467-023-42824-z).
This draft is an original synthetic extension; neither citation validates its
specific dataset, identifiability, thresholds, or difficulty. It is not claimed
to be a previously unsolved scientific discovery.

Release remains blocked until reference and adversarial verifier checks pass,
the documented contract is complete, and substantive model trials support the
required difficulty. Any successfully solved candidate must remain unpackaged.
