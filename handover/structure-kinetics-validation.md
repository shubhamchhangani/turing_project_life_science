# Structure–kinetics prototype progress

This is a research prototype, not a validated Turing task or a demonstrated
model-breaking instance. No frontier trial or ZIP has been produced.

## Completed work

- Parameterized forward model for four conformations, cryo-like projections,
  temperature-dependent detailed-balance kinetics, cooling and dye photophysics.
- Reproducible public/private data separation: 128 images and 24 traces.
  Initial pilot has 31,891 photons; longer pilot has 317,957 photons.
- Three-way photon-time likelihood validation (spectral, direct exponential,
  uniformization), including no-photon survival: max log difference 2.85e-14.
- Projection versus independent 3D numerical integration: max difference 6.94e-18.
- Input-only photon fits from three independent initializations.
- Input-only image geometry and orientation recovery, followed by crossmodal
  initialization. All four domains recovered within 0.092 angstrom RMSD.
- Broadening pose search recovered all 128 images to noise-level residuals.
  A subsequent, separate truth audit found 128/128 correct conformation
  assignments and median rotation error 0.1133 degrees.

No true poses, state labels, geometries or kinetic parameters initialized these
inverse solvers. Truth was inspected only by author diagnostics after fitting.
The known rigid core is an explicit pilot simplification; successful recovery
must not be presented as proof that the final inverse task is difficult.

## Important unresolved finding

Good likelihood does not guarantee precise kinetic parameters on these data.
Short-trace local curvature showed very weak rate constraints. Ten times more
photons reduced the weakest local uncertainty, but a roughly ten-fold local
uncertainty factor remains for some rates. This is an inverse-Hessian diagnostic
conditional on estimated geometry, not a calibrated confidence interval.
Longer-trace fits can differ from planted rates by around 6–14 fold.
It would be invalid to grade those rates against arbitrary tight truth tolerances.
The full coupled fit and data-design/identifiability checks must address this.

## In progress at report creation

A coupled kinetic objective incorporates all images, their inferred poses and
structures, cooling-dependent mixture populations and raw photon likelihoods.
Geometry and poses are fixed in this pilot; joint refinement and uncertainty
marginalization remain outstanding. Machine-readable current status is in
`structure-kinetics-progress.json`.

## Remaining release gates

Full joint reference and calibrated uncertainty; a genuinely harder but
identifiable structural instance; independent output verifier; one-CPU/2-GiB
runtime validation; schema and adversarial-output audits; then Claude first.
GPT is allowed only after a substantive Claude failure. ZIP release requires
both normal failures on identical frozen task bytes. No all-model-failure or
portal-acceptance guarantee is made.
