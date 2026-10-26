# Joint structure–kinetics inference across cryo-EM and photon-resolved smFRET

Status: prototype data, forward checks, and input-only inverse baselines now
implemented. There is no full joint inverse solution, validated task package,
frontier-model trial, or ZIP. No model-failure claim is justified.

## Scientific question

Recover a synthetic multidomain protein's conformational ensemble and kinetic
energy landscape jointly from heterogeneous cryo-EM projections acquired under
two cooling protocols and photon-resolved smFRET records at two temperatures.
The same conformations and transition parameters must explain both modalities.

This is a bounded synthetic inverse problem motivated by real structural
biology. It is not a claim to solve an open natural-protein research question.

## Why this family deserves a different screen

Previous receptor tasks reduced to a solvable planar core and exceptional-site
reweighting. The lineage task reduced to one augmented sparse generator and
message contraction. Both supplied the mechanism's numerical parameters.

Here unknown structures determine image projections AND fluorescence emissions.
Unknown kinetics determine photon histories AND cooling-distorted cryo-EM
mixture weights. Unknown image poses, class identities, conformational states
and dye states must be marginalized or optimized consistently. The central
bottleneck is recovering these unknowns from indirect measurements, rather than
evaluating a fully specified forward model once.

## Provisional finite specification

These are author design choices, not a finalized task contract:

- Four conformations of a coarse-grained multidomain protein.
- Public rigid-domain templates; unknown relative domain placements in each
  conformation. An asymmetric core fixes the coordinate gauge where possible.
- Approximately 512 noisy 2D projections, split between two public cooling
  protocols. Each image has unknown conformation, rotation and in-plane shift.
- Explicit, calibrated image formation: sums of projected Gaussian densities,
  known per-image transfer functions, and a specified noise model. Independent
  numerical integration must validate the analytic projection implementation.
- Multiple labeled FRET constructs observing the same conformational process,
  at two temperatures. Public labeling sites, Förster distances and observation
  calibration; raw donor/acceptor photon times rather than inferred states.
- A calibrated dye blinking process, so apparent fluorescence switching need
  not mean a protein transition. Include an acceptor-excitation stream only if
  the ablation checks establish that it resolves a relevant ambiguity.
- Conformational transition rates
  q(i->j,T) = nu * exp[(E_i-B_ij)*Tref/T], with B_ij=B_ji,
  fixed public nu, E_0=0, and public bounds B_ij>max(E_i,E_j).
  Thus pi_i(T) is proportional to exp(-E_i*Tref/T).
- A public cooling trajectory T(t) maps room-temperature populations to frozen
  populations through dp/dt = Q(T(t))p. No assumption that frozen populations
  equal room-temperature or final-temperature equilibrium is allowed.
- The number of states is stated. Do not add model selection merely to increase
  the number of requirements before identifiability is demonstrated.

The exact domain templates, counts, noise and parameter bounds remain undecided
until a genuine reconstruction and uncertainty analysis succeeds.

## Shared likelihood and dependence

Image likelihood sums over conformations with weights from the cooling ODE and
integrates/optimizes pose according to the explicitly chosen statistical contract.
For a fixed conformation, projections depend on that conformation's domain
placements. Those same placements determine FRET distances and photon rates.
Photon likelihoods marginalize the conformation-by-dye continuous-time chain,
including no-photon intervals, observations, and any censoring at trace end.

An image-only reconstruction cannot automatically identify the kinetic network.
Histogram-only fluorescence inference can confuse dye switching and molecular
switching. Separate modality fits require a consistent mapping of states and
populations; simply concatenating unrelated answers does not solve the joint
model. Whether the chosen synthetic instance actually requires joint fitting
must be established by ablation, not assumed.

## What exponential means here

For K states and N images, explicit class assignment has K^N combinations.
An explicitly enumerated hidden trajectory of length L has K^L state strings;
unknown rotations add a continuous search space. These are sizes of naive
searches, NOT a proof that every algorithm requires exponential time.
Hidden Markov summation can collapse the trajectory enumeration efficiently.
The research challenge is the joint nonconvex reconstruction, parameter
identifiability and competing solutions. Do not advertise an exponential lower
bound or model-breaking difficulty without evidence.

## Reference and validity gates before any model trial

1. Generate public observations and retain truth separately. Never use true
   poses, classes or parameter values to initialize the reference inverse solver.
2. Build a genuine input-only reconstruction using global/multistart geometry
   optimization, pose refinement, photon likelihood and constrained kinetic fit.
3. Validate separate forward implementations: analytic versus integrated
   projection, adaptive ODE versus time-ordered exponentials, small-state photon
   likelihood versus independently enumerated/uniformized calculations.
4. Check gauge freedoms, state permutations, reflections, unobserved directions
   and practical identifiability. Do not reject scientifically equivalent fits.
5. Demonstrate all required observables recover within tolerances calibrated by
   independent simulation/noise replicates. Do not demand exact planted
   parameters when finite data support multiple explanations.
6. Compare simple baselines: single-state reconstruction, image-only fitting,
   FRET histograms, independent modality fits, no-blinking model, and ignoring
   cooling. If a straightforward baseline succeeds, reject the candidate.
7. Verify the authentic reference under the task's actual CPU/RAM/time limits.
   Do not make infeasible computation the intended source of failure.
8. Freeze a complete schema-1.4 package, full-coverage outcome-only tests, safe
   regular-file reads, strict decimal grammar, independent expectation
   provenance, and adversarial wrong-solution tests.
9. Claude Opus 5.5 low first. GPT-6-sol max only after a normal substantive
   Claude failure. Direct-root Downloads ZIP only after both fail on unchanged
   bytes and every scientific/format/offline-verifier gate passes.

## Initial executable screen

`screen_cooling.py` checks a four-state detailed-balance model with nine energy
parameters (three independent well energies and six barriers). It compares
adaptive integration against successively refined time-ordered exponentials
and measures the deviation of frozen populations from equilibrium assumptions.
This tests only the forward coupling. It supplies no inverse solution, no
evidence about identifiability from noisy observations, and no model-failure
evidence. Results are recorded in `draft/cooling-screen.json`.

## Primary research grounding

- Levy et al., *Amortized Inference for Heterogeneous Reconstruction in Cryo-EM*:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10392957/
  Supports the joint unknown conformation/pose reconstruction setting.
- Harris et al. (2022), *Multi-parameter photon-by-photon hidden Markov modeling*:
  https://www.nature.com/articles/s41467-022-28632-x
  Demonstrates that dye photophysics and conformational transitions can be
  intertwined and that additional photon streams can distinguish them.
- *Effects of cryo-EM cooling on structural ensembles* (2022):
  https://www.nature.com/articles/s41467-022-29332-2
  Studies how cooling and kinetics alter the structural ensemble.

These papers motivate component physics. They do not establish that this newly
proposed combined benchmark is identifiable, feasible, or hard for the roster.

## Implementation progress

- `generate.py`: reproducible public observations and separately stored private
  truth. Initial pilot: 128 images, 24 traces, 31,891 photons. A longer-trace
  calibration pilot contains 317,957 photons with byte-identical image data.
- `model.py`: parameterized structure, projection, temperature-dependent kinetics,
  cooling and conformation-by-dye photon physics; no planted parameters.
- `check_forward.py`: analytic projection versus independent 3D quadrature;
  photon likelihood via reversible eigendecomposition, direct matrix
  exponentials and positive uniformization, including empty-trace survival.
  Differences: 6.94e-18 in projection values, 2.85e-14 in log likelihood.
- `fit_photons.py`: three public-data-only multistart fits. High photon likelihood
  can coexist with wrong planted geometry; this is an ambiguity diagnostic,
  not proof of computational hardness or a basis to reject equivalent fits.
- `fit_images.py`: 128 public images screened without true poses/classes;
  57 fit at reduced chi-square below 1.5; the remainder include local minima.
- `fit_crossmodal_pilot.py`: clusters reliable image-derived geometries, then
  fits photon kinetics from independent initializations. No truth initialization.
  It is a sequential initialization baseline, not the full joint objective.
- Post-fit diagnostic: all four mobile-domain geometries recovered within
  0.092 angstrom RMSD, but some fitted rates differ from generating values by
  about 4.4-fold on the short traces. Truth was read only after fitting.
- `check_identifiability.py`: finite-difference local curvature agrees across
  step sizes. Some rates have very large local uncertainty; short-trace
  tight parameter grading would be unjustified. The inverse-Hessian diagnostic
  is not a calibrated confidence interval or global identifiability proof.

Outstanding: full joint geometry/pose/kinetic/cooling fit, noise-replicate and
identifiability calibration, harder but still identifiable structural instances,
resource-constrained reference, complete verifier and conditional model trials.
The simple known-core hinge pilot is a numerical foundation, not a demonstrated
hard benchmark; its successful image reconstruction is not hidden or discounted.
