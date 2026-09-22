# Proposed CMB contribution

Upstream: `knygren/glmbayes`, `main` at `7959a42d32f822a4b65ddfaf94927190ad5db536`
(verified against GitHub on 2026-09-21). Local starting point: CMB v1 at
`611a88fa960d1617a1c6c2416caabaa0b0fe8920`, package version 0.9.76.9001.
Prepared for the personal fork of `jbogomolovas2`.

This branch starts directly from upstream. It does not include private
development commits, raw emu observations, study reports, machine-specific
build configuration or discarded adaptation experiments. The upstream author,
package identity and license are retained. The README identifies this as a fork.

## CMB functionality

* A finite-support likelihood with natural predictors theta = X beta and
  nu = Z gamma; both use identity links, with a joint normal coefficient prior.
* `cmb()` and `cmb_augment()` for the augmented `rglmb` matrix interface.
* `cmb_fitted()` integrates the PMF per draw for expected counts/proportions;
  `cmb_nu()` extracts nu. CMB-aware fitted/residual/summary handling is included.
* Documentation and a synthetic example with separate predictor designs.

This is a natural-parameter model, not a mean-link CMB model. Nu is unconstrained
on the real line (finite support gives a normalizable distribution). Random
effects and generic stats GLM-family compatibility are not provided.

## Shared sampling changes

* Safeguarded Newton tangency refinement with final residual reporting and a
  bounded-memory fixed-point fallback for larger problems.
* Log-space mixture weights/selection and stable extreme-tail proposal bounds.
* Bounded Normal-GLM rejection work; failure returns no partial accepted sample.
* CPU row-indexed RNG streams controlled by R seeds, including independent
  Normal-Gamma sampling, unaffected by worker scheduling.
* Fit diagnostics retained through matrix/formula interfaces and summaries.

These shared changes warrant review separately from the distribution itself.
There is no automatic grid adaptation, cost controller or new tuning layer.

## Validation

Fresh public-tree validation on 2026-09-21 passed 302 expectations, with zero
failures or test warnings. Eleven unavailable OpenCL test blocks were skipped.
This includes direct PMF/derivative/binomial-limit checks and shared regressions.
All six synthetic cases passed in both CPU modes: 12 fits, 240,000 posterior
draws, checked against numerical integration references. The synthetic example
also completed. See `inst/validation/CMB_VALIDATION.md` to reproduce the checks.

CPU validation does not certify OpenCL. Finite tests do not prove correctness
for every input or practical envelope cost at arbitrary dimension. The capped
proposal count does not cap time. Exact RNG results are not promised across
different platforms or numerical-library builds. No broad speed comparison
with glmmTMB is claimed.

## Suggested review boundaries

1. Shared numerical/sampling changes with their regression tests.
2. CMB likelihood, augmentation and post-fit methods.
3. Fit diagnostics and documentation.

The branch is prepared as one reviewable snapshot. These are proposed review
boundaries, not claims that every subset cherry-picks independently. The
maintainer's preferred interface and contribution structure should guide PRs.
