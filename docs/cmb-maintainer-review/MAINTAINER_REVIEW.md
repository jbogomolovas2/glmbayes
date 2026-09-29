# Maintainer review: CMB contribution and tangency refinement

### Review at a glance

The contribution adds CMB regression and several shared sampler improvements. Here is what changes:

1. **[Signatures](#1-changed-signatures):** `EnvelopeSort` adds `logPLSD`; the internal/native envelope builder adds `refine` and `refine_maxit`. Most public fitting arguments stay the same. The handwritten envelope wrapper needs a missed update: it passes 16 arguments to an 18-argument native function and can crash R.
2. **[Returned items](#2-additional-returned-items):** envelopes gain log weights and refinement diagnostics; fits and summaries gain sampling diagnostics. The tables explain their types, availability and ordering.
3. **[New functions/files](#3-new-functions-and-files):** CMB adds family, data-preparation and fitted-value helpers, plus likelihood kernels. Shared additions cover refinement, tail calculations, random-number streams and diagnostics.
4. **[CMB implementation](#4-cmb-implementation):** two predictors model `theta` and `nu`, using the existing normal-prior sampler. The supported entry point is the augmented matrix interface; formula support is incomplete.

In the zero-heavy CMB example, refinement allowed all three 20,000-draw runs to finish; without it, all three hit the proposal cap before accepting a draw. The other five CMB examples worked either way. [Results in brief](EVIDENCE_SUMMARY.md) and [datasets/model specifications](cases/README.md) provide the context.

Refinement is currently enabled by default in the shared sampler. Its effect on the paper's asymptotic guarantees remains unresolved.

### Versions and supporting material

This review compares contribution **`b971596511c90f33de290497b2c7e8feb1a4e580`** with base **`7959a42d32f822a4b65ddfaf94927190ad5db536`**: 51 changed paths, 2,583 insertions and 234 deletions. Source links refer to that contribution. The tests used a clean build and public or generated data.

The four sections below provide the detailed inventory. Longer signature and file tables can be expanded as needed. [Full diff](reviewed.diff) · [Function inventory](r-function-inventory.csv) · [Results](results/summary.csv) · [Publication notes](PUBLICATION_NOTES.md)

## 1. Changed signatures

### Public R functions

Only `EnvelopeSort` changes the formal argument list of an existing exported R function. It appends `logPLSD = NULL` to `(l1, l2, GIndex, G3, cbars, logU, logrt, loglt, logP, LLconst, PLSD, a1, E_draws, lg_prob_factor = NULL, UB2min = NULL)`. Existing positional/named calls remain valid; the new optional vector is reordered with `PLSD` and returned only when supplied. Both regular and independent Normal-Gamma envelope builders call it; the latter need not supply the new field. [Source](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R)

The following existing public functions have **unchanged arguments and defaults**, but changed behavior. The exact old/new formals for each are recorded in the signature appendix below and CSV.

| Functions | Callers and compatibility implications |
|---|---|
| `dNormal`, its internal `plinks` closure | Prior constructors now admit CMB/identity; other prior-family dispatch is unchanged. |
| `glmbfamfunc`, family-local `f1`, `f2`, `f3`, `f4`, `f7` | Used by mode optimization, deviance and information calculations; adds CMB implementations of the existing callback contracts. The branch-local definitions are new implementations, not changes to existing families' callback arguments. |
| `rNormal_reg` | Called through `rglmb`/prior simulation dispatch; admits CMB, skips `glmb.wfit` for it, removes `fit` by assigning NULL. Existing post-fit consumers that require an IRLS fit cannot assume it exists. |
| `rindepNormalGamma_reg` | Called by `rglmb`/`rlmb`; now returns the actual `Envelope`, previously NULL. More retained memory; callers testing `is.null(Envelope)` observe a change. |
| `rglmb`, `rlmb` | Called by users and formula interfaces; add diagnostics and structured proposal-cap errors via `.glmb_run_simulation`. Public controls still do not expose refinement. Existing parallel/serial controls and defaults are unchanged. |
| `glmb`, `lmb` | Formula interfaces copy simulation diagnostics. No complete CMB formula interface is added: the supported CMB entry is augmented matrices with `rglmb`. |
| `residuals.rglmb` | User/S3 residual calls explicitly error for CMB; per-observation deviance residuals remain unsupported. |
| `summary.rglmb` | User/S3 summary calls use CMB finite-sum fitted means; suppress GLM MLE/SE values with NA and return NULL model for CMB. Carries diagnostics. |
| `summary.glmb`, `print.summary.glmb`, `print.summary.rglmb` | Carry/print structured diagnostics when present; otherwise retain the old candidate-count print. Printed output changes without argument changes. |

The internal `.uni_lmb` also copies diagnostics without changing its formals. No other existing R function definition changes in the parsed diff. The parser records the final same-named family closure in each version; the CMB `f1`–`f7` rows therefore denote branch additions, not replacements of earlier branches. The full diff is authoritative for these branches.

### R wrapper → registration → C++ trace, including the defect

The unchanged public formals are `EnvelopeBuild(bStar, A, y, x, mu, P, alpha, wt, family="binomial", link="logit", Gridtype=2L, n=1L, n_envopt=NULL, sortgrid=FALSE, use_opencl=FALSE, verbose=FALSE)` in both revisions. The unchanged manual `.EnvelopeBuild_cpp(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)` has no defaults in either revision. Its unchanged signature is nevertheless affected by the native arity change.

- `rglmb` → `rNormal_reg` → `.rNormalGLM_cpp` → generated `_glmbayes_rNormalGLM_cpp_export` (18 arguments, unchanged) → `rNormalGLM_cpp_export` → `glmbayes::sim::rNormalGLM` → `glmbayes::env::EnvelopeBuild(..., true, 60)`. This path bypasses the defective R envelope wrapper. [R dispatcher](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simfunction.R), [native wrapper](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/export_wrappers.cpp), [sampler](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalGLM.cpp).
- Internal generated `EnvelopeBuild_cpp_export` changes from 16 arguments to the same 16 plus `refine = TRUE, refine_maxit = 60L`. Its `.Call`, generated C trampoline, native registration (16 → 18), export wrapper and C++ implementation all agree. A call to this R function using the old 16 arguments receives its R defaults and is source-compatible. [Generated R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/RcppExports.R), [generated C++](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/RcppExports.cpp).
- **Defect D1:** public `EnvelopeBuild` has unchanged formals; for non-Gaussian families it calls `.EnvelopeBuild_cpp`, whose unchanged manual wrapper passes only 16 arguments directly to the now-18-argument native symbol. C++ default arguments do not apply to an R `.Call` pointer invocation. Native registration reports 18; `length(formals(.EnvelopeBuild_cpp))` is 16. On this build the isolated reproduction **segfaulted** at this call, rather than producing a catchable R error. Gaussian public `EnvelopeBuild` instead uses the separate independent Normal-Gamma builder. [Manual wrapper](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rcpp_wrappers.R), [public wrapper](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R), [crash log](regression-export.log), [safe arity inspection](defects.log).
- The comparisons use the generated export, which has the correct argument list.

### Existing C++ interfaces

| Function / layer | Old → new arguments/defaults | Callers / implications |
|---|---|---|
| `glmbayes::env::EnvelopeBuild` | Append `bool refine=true, int refine_maxit=60` after `verbose=false`; all earlier parameters/defaults unchanged. | Native GLM sampler explicitly passes true/60; export wrapper forwards settings. Source callers including the header can omit them. C++ ABI changes; recompile all dependents. |
| `EnvelopeBuild_cpp_export` | Append `bool refine=true, int refine_maxit=60` after mandatory `verbose`; preceding 16 arguments unchanged. | Rcpp-generated R/C++ wrappers; generated native entry must receive 18. |
| `_glmbayes_EnvelopeBuild_cpp_export` | 16 SEXP parameters → 18, appended `refineSEXP`, `refine_maxitSEXP`; no C defaults. | Registration count changes accordingly. Manual wrapper defect above. |
| `rNormalGLM_worker` constructor | Add mandatory `uint64_t job_seed_` after `n_`; add `RVector<double> logPLSD_r_`, `bool use_log_sel_` after `PLSD_r_`; add `RMatrix<double> logU_r_` after `logrt_r_`. | Internal parallel/pilot construction sites updated. No R API change. `operator()` arguments stay `(size_t begin,size_t end)`; per-row streams, cap failure and CMB scoring change its body. |
| `rIndepNormalGammaReg_worker` constructor | Add mandatory `uint64_t job_seed_` after `n_`, `RMatrix<double> logU_r_` after `logrt_r_`. | Internal parallel caller updated. `operator()` signature unchanged. |
| `glmbayes::rng::rnorm_ct` | Preserve `(double lgrt,double lglt,double mu,double sigma)`; add overload with mandatory final `double log_mass`. No defaults. | Normal-GLM and independent Normal-Gamma serial/parallel paths supply `logU`; legacy four-argument callers remain possible but cannot encode both extreme-tail endpoints reliably. |
| `runif_safe()` | Signature unchanged; engine changes from OS-seeded `mt19937`/standard distribution to per-draw-seeded `mt19937_64`, explicit nonzero 53-bit uniform. | Shared envelope samplers and truncated draws; results differ from the base's RNG streams. |
| `f2_f3_non_opencl`, `rNormalGLM`, `rNormalGLM_std`, `rNormalGLM_std_parallel`, `run_rcppparallel_pilot` | Signatures/defaults unchanged. | CMB dispatch; log-weight/tail safeguards; proposal cap now aborts without force acceptance; reproducible CPU streams. Pilot checks report cap failure after worker return. |
| `rIndepNormalGammaReg`, `rIndepNormalGammaReg_std`, `rIndepNormalGammaReg_std_parallel` | Signatures/defaults unchanged. | Return envelope; use seeded row streams and five-argument truncated normal. |
| `rNormalReg` | Signature/defaults unchanged. | Adds explicit `Envelope=NULL` to direct Gaussian output; diagnostics can identify direct sampling. |

[Envelope declaration](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/Envelopefuncs.h), [Normal GLM implementation](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalGLM.cpp), [independent Normal-Gamma implementation](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rIndepNormalGammaReg.cpp), [RNG declarations](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rng_utils.h), [tail sampler](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rnorm_ct.cpp). There are no new CMB `.Call` registrations: new family kernels are reached by existing dispatch.

Exact old/new C++ formal lists (including worker constructors and unchanged native interfaces) are in [cpp-function-inventory.csv](cpp-function-inventory.csv). Declarations record defaults; definitions do not repeat header defaults.

### Exact R signature appendix

<details>
<summary>Show the complete R signature table</summary>

Rows below include all changed existing R definitions, including body-only changes. An argument without a default is required. New definitions are catalogued in section 3 and the CSV.

| Function | Old formals | New formals |
|---|---|---|
| [`EnvelopeBuild_cpp_export`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/RcppExports.R) | `(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose)` | `(bStar, A, y, x, mu, P, alpha, wt, family, link, Gridtype, n, n_envopt, sortgrid, use_opencl, verbose, refine = TRUE, refine_maxit = 60L)` |
| [`glmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/glmb.R) | `(formula, family = binomial, pfamily = dNormal(mu, Sigma, dispersion = 1), n = 1000, data, weights, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, subset, offset, na.action, Gridtype = 2, n_envopt = NULL, start = NULL, etastart, mustart, control = list(...), model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL, ...)` | unchanged (same arguments/defaults) |
| [`lmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/lmb.R) | `(formula, pfamily, n = 1000, data, subset, weights, na.action, method = "qr", model = TRUE, x = TRUE, y = TRUE, qr = TRUE, singular.ok = TRUE, contrasts = NULL, offset, Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, ...)` | unchanged (same arguments/defaults) |
| [`.uni_lmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/lmb.R) | `(formula, pfamily, n = 1000, data, subset, weights, na.action, method = "qr", model = TRUE, x = TRUE, y = TRUE, qr = TRUE, singular.ok = TRUE, contrasts = NULL, offset, Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, ...)` | unchanged (same arguments/defaults) |
| [`dNormal`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/pfamily.R) | `(mu, Sigma, dispersion = NULL)` | unchanged (same arguments/defaults) |
| [`dNormal::plinks`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/pfamily.R) | `(family)` | unchanged (same arguments/defaults) |
| [`residuals.rglmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/residuals.glmb.R) | `(object, ysim = NULL, ...)` | unchanged (same arguments/defaults) |
| [`rglmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rglmb.R) | `(n = 1, y, x, family = gaussian(), pfamily, offset = NULL, weights = 1, Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE)` | unchanged (same arguments/defaults) |
| [`rlmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rlmb.R) | `(n = 1, y, x, pfamily, offset = rep(0, nobs), weights = NULL, Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, progbar = FALSE)` | unchanged (same arguments/defaults) |
| [`rindepNormalGamma_reg`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simfunction.R) | `(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(), Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, progbar = TRUE)` | unchanged (same arguments/defaults) |
| [`rNormal_reg`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simfunction.R) | `(n, y, x, prior_list, offset = NULL, weights = 1, family = gaussian(), Gridtype = 2, n_envopt = NULL, use_parallel = TRUE, use_opencl = FALSE, verbose = FALSE, progbar = FALSE)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(family, lik_shape = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f1`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, alpha = 0, wt = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f2`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, mu, P, alpha = 0, wt = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f3`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, mu, P, alpha = 0, wt = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f4`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, alpha = 0, wt = 1, dispersion = 1)` | unchanged (same arguments/defaults) |
| [`glmbfamfunc::f7`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(b, y, x, mu, P, alpha = 0, wt = 1)` | unchanged (same arguments/defaults) |
| [`EnvelopeSort`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | `(l1, l2, GIndex, G3, cbars, logU, logrt, loglt, logP, LLconst, PLSD, a1, E_draws, lg_prob_factor = NULL, UB2min = NULL)` | `(l1, l2, GIndex, G3, cbars, logU, logrt, loglt, logP, LLconst, PLSD, a1, E_draws, lg_prob_factor = NULL, UB2min = NULL, logPLSD = NULL)` |
| [`summary.glmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.glmb.R) | `(object, ...)` | unchanged (same arguments/defaults) |
| [`print.summary.glmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.glmb.R) | `(x, digits = max(3, getOption("digits") - 3), ...)` | unchanged (same arguments/defaults) |
| [`summary.rglmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.rglmb.R) | `(object, ...)` | unchanged (same arguments/defaults) |
| [`print.summary.rglmb`](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.rglmb.R) | `(x, digits = max(3, getOption("digits") - 3), ...)` | unchanged (same arguments/defaults) |

</details>

## 2. Additional returned items

### Envelopes and sampler returns

| Field | Type / availability | Meaning and consumers |
|---|---|---|
| `Envelope$logPLSD` | Numeric vector, one entry per cell, regular `EnvelopeBuild` sorted, unsorted and sort-memory-fallback paths | Normalized log mixture weights. CPU serial/parallel Normal-GLM selection uses these when a finite log weight exponentiates to zero. Ordinary weights still use the bounded ordinary scan. `EnvelopeSort` returns it only if supplied; no blanket promise for other builders. |
| `Envelope$refinement` | List, regular builder on all return paths; empty when disabled (also potentially incomplete after unusable initial mass) | Detailed refinement state below. Fit diagnostics summarize it; public tests and this harness inspect it. |
| `Envelope` in native `rIndepNormalGammaReg` and R `rindepNormalGamma_reg` | Actual list `Env3`, previously absent/native or NULL/R | Shared diagnostics can now distinguish rejection sampling from direct sampling. It is a different envelope builder: do not assume regular-builder refinement fields exist. |
| `Envelope` in native `rNormalReg` | Explicit NULL | Conjugate Gaussian path, no envelope/refinement. |
| Existing `thetabars`, `cbars`, `LLconst`, `logP`, `PLSD` | Same names; values now reflect refined tangencies | `G3` is synchronized with refined `G4` before `setlogP_C2`, avoiding stale tangencies paired with new gradients. `logP` returned here is the unnormalized box log probability (first internal column), not the full cell envelope mass. |

The full cell log mass used in this review is `logP - LLconst + rowSums(cbars^2)/2`; total mass is its log-sum-exp. These are masses against the standardized normal base, with the likelihood's retained constants. They are comparable between on/off for a given case, not as evidence across different datasets. Acceptance prediction uses the matching integration normalization.

`refinement` fields (from [implementation](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/envelope_refine.cpp)):

| Field(s) | Type / meaning |
|---|---|
| `iterations` | Integer evaluated update passes, at most configured maxit (60 here). |
| `residual` | Numeric scalar maximum absolute `m(r)-r` over axes/cells, final state. |
| `converged` | Logical `residual < 1e-6`; efficiency diagnostic, not an MCMC convergence test or accuracy certificate. |
| `accelerated` | Logical eligibility for Newton path; does not mean every cell accepted a Newton update. |
| `log_mass_initial`, `log_mass_final` | Numeric scalars, log-sum-exp of cell envelope masses before/after. |
| `cell_log_mass` | Numeric vector of final per-cell log masses. |
| `step_fraction` | Numeric vector of final per-cell damping factors. |
| `accepted_updates`, `rejected_updates`, `newton_updates` | Numeric vectors of per-cell update counts; accepted Newton steps are a subset of accepted updates. |
| `cell_grid_index` | Numeric matrix encoding integer region codes; builder attaches original `GIndex` to identify the diagnostic cells. |

**Sorted-path caveat:** `EnvelopeSort` reorders main arrays and `logPLSD` by the same permutation. The builder attaches `refinement` afterward; its vectors remain in **original grid order**. Join them to sorted envelope rows using `cell_grid_index` versus `GridIndex`; do not zip them by position. Unsorted and memory-fallback returns carry both fields as well. [Harness check](harness-check.log) verifies normal sorted/unsorted alignment; the allocation-failure fallback is inspected in source, not forced at runtime. The fallback rebuilds an unsorted list and does not retain a `sort_ok` field.

### Fit, summary and error diagnostics

`rglmb`, `rlmb`, `glmb`, `lmb`/`.uni_lmb`, `summary.glmb`, and `summary.rglmb` add/copy `diagnostics`, a list of class `glmb_diagnostics`. The low-level `rNormal_reg` called directly does not itself invoke the new diagnostic wrapper. The two summary printers consume this object; `print.glmb_diagnostics` is a new S3 method. [Definitions](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/fit_diagnostics.R)

| Field(s) | Type / availability / meaning |
|---|---|
| `status`, `method`, `sampling_mode` | Character scalars: completed/refinement_incomplete; direct/envelope_rejection; serial/parallel/not_applicable. Parallel describes the requested CPU mode, not an independently measured device path. |
| `requested_parallel`, `requested_opencl` | Logical requested flags; no proof of actual OpenCL execution. |
| `draws`, `total_proposals`, `proposals_per_draw`, `acceptance_rate`, `max_proposals` | Numeric/integer counts or rates; proposal fields NA if no envelope or `iters` is not a positive integral count per returned draw. |
| `envelope$status`, `$cells`, `$iterations`, `$residual` | Character status and integer/numeric summaries. Status can be converged, incomplete, not_reported, not_applicable; absent refinement values become NA. |
| `build$version`, `$library`, `$R_version`, `$platform` | Character build identity; library is local installation path. |
| `build$fingerprint` | Named character MD5 vector `compiled_code`, `r_code`; NA for unavailable files. Cached once per namespace. |
| `rng_kind` | Character vector from `RNGkind()`. |

Cap errors passing through `.glmb_run_simulation` have classes `glmbayes_sampling_error`, `error`, `condition`, and fields `message` (character), `call`, `parent` (original error), `diagnostics`. Error diagnostics include `status='failed'`, `reason='proposal_limit'`, `draws_returned=0L`, numeric `proposal_limit`, `total_proposals=NA`, character `sampling_mode`, `action`, and `build`. Native low-level harness calls receive the original error, not this R wrapper class. No partial samples are returned. Other errors are rethrown unchanged.

### New CMB return contracts

`cmb()` returns a class-`family` list with `family='cmb'`, `link='identity'`, `valideta` and `validmu` functions (both return TRUE). Despite that class it is not a complete `stats::family` object. `cmb_augment` returns numeric `x` (2N × (p+q)), `y`, `weights`, `offset` (each length 2N); details in section 4. `cmb_fitted` and `cmb_nu` each return draw × N numeric matrices. `glmbfamfunc(cmb())` returns the existing callback names, implemented for CMB. Native `f2_f3_cmb` returns `qf` (one numeric value per parameter column) and `grad` (grid-points × parameters matrix), consumed by `EnvelopeEval` as `NegLL` and `cbars`.

For CMB `summary.rglmb` additionally changes existing values: `fitted.values` is draw × N, `linear.predictors` still contains 2N augmented entries; `family` is the minimal CMB family; `model` is NULL; MLE coefficient/SE comparison columns are NA.

## 3. New functions and files

### New functions and dependencies

| Group | New names, purpose and visibility |
|---|---|
| CMB public R | `cmb(link='identity')`, `cmb_augment(X,yc,m,Z=NULL)`, `cmb_fitted(object,scale=c('proportion','count'))`, `cmb_nu(object)` are exported. Base R finite sums and matrix operations; sampler requires existing normal prior and standardization. |
| CMB internal R | `.is_cmb`, `.cmb_nobs`, `.cmb_m`; `glmbfamfunc` closures `.cmb_N`, `.cmb_tab`, `.cmb_mom`, `.cmb_dat`, `.cmb_eta`, `.cmb_sat` and its `EYf`; CMB-specific `f1/f2/f3/f4/f7`. Layout checks, sufficient-statistic moments, likelihood/prior derivatives, saturated-likelihood solve. |
| CMB C++ | `f1_cmb`, `f2_cmb`, `f3_cmb`, `f2_f3_cmb`, `f2_cmb_rmat`; private `CMBData`, `cmb_prep_t`, `cmb_moments`, `cmb_N`. Internal namespace kernels, declared in `famfuncs.h`; Rcpp/Armadillo/Rmath/RcppParallel container dependencies. No new exported Rcpp entry points. |
| Shared refinement | `refine_tangency`; private bounds/target/residual helpers and lambdas. Internal C++ API; depends on `EnvelopeEval`, Rmath tails and Armadillo linear solves. |
| Numerical stability | `mills_ratio`, `trunc_norm_mean`, `log_box_prob`, `logaddexp_safe`, `select_cell_log`; private `R_pos`, `logaddexp`, `cell_bounds`; new `rnorm_ct` overload. Shared C++ utilities, not R exports. |
| Shared RNG/cap | `seed_from_R`, `seed_draw`, private `check_proposal_limit` and cap constant. Seed extraction only on R thread; per-row C++ streams. |
| Diagnostics | `.glmb_build_info`, `.glmb_fit_diagnostics`, `.glmb_run_simulation`, cache environment; exported S3 `print.glmb_diagnostics`. Base R build hashes and count checking; no new package dependency. |
| Validation | `quadratic_envelope` test helper; matrix script `lse`, `quad`, `importance` and local callbacks; six new test files and Python subprocess runner. Tests use testthat; importance fallback also requires **mvtnorm**, used via `::` but not declared in DESCRIPTION. |

CMB's finite-support likelihood, two-predictor augmentation, dispatch and mean conversion are its requirements. Refinement, log mixture weights, tail stability, row-seeded RNG, shared proposal failure handling and diagnostics are **broader proposed improvements**. The comparisons below test refinement while keeping the other changes in place.

### Complete revision file ledger

<details>
<summary>Show all 51 changed paths</summary>

All 51 paths are listed below. `A/M/D` means added/modified/deleted; generated files are labeled. [Download the file list](diff-inventory.tsv).

| Status | File | Area / visibility and dependencies |
|---|---|---|
| M | [.Rbuildignore](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/.Rbuildignore) | **validation / packaging** — Exclude local analyses, backups and generated machine Makevars from builds; R build tooling. |
| M | [.gitignore](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/.gitignore) | **validation / packaging** — Ignore local analyses, editor/build files; Git only. |
| A | [CMB_CONTRIBUTION.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/CMB_CONTRIBUTION.md) | **validation** — Contribution scope and historical public results; documentation only. |
| M | [DESCRIPTION](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/DESCRIPTION) | **CMB / packaging** — Development version 0.9.76.9001 and X-CMB-Release v1; no added dependencies. |
| M | [NAMESPACE](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/NAMESPACE) | **CMB / diagnostics** — Generated-style namespace file with CMB exports and appended diagnostic S3 registration. |
| M | [R/RcppExports.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/RcppExports.R) | **shared sampling** — GENERATED Rcpp R wrapper; refinement arguments/defaults; Rcpp native symbol. |
| A | [R/fit_diagnostics.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/fit_diagnostics.R) | **diagnostics** — New fit/build/error diagnostics and S3 printer; base R hashes, RNGkind. |
| M | [R/glmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/glmb.R) | **diagnostics** — Formula GLM diagnostic propagation; rglmb. |
| M | [R/lmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/lmb.R) | **diagnostics** — Formula linear/univariate helper diagnostic propagation; rlmb. |
| M | [R/pfamily.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/pfamily.R) | **CMB** — Normal-prior dispatch, public augmentation/fitted helpers and family; finite sums. |
| M | [R/residuals.glmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/residuals.glmb.R) | **CMB** — Explicit unsupported CMB residual error; S3 methods. |
| M | [R/rglmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rglmb.R) | **diagnostics / shared sampling** — Diagnostic wrapper and reproducibility docs; simfunction. |
| M | [R/rlmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/rlmb.R) | **diagnostics / shared sampling** — Diagnostic wrapper; simfunction. |
| M | [R/simfunction.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simfunction.R) | **CMB / diagnostics** — CMB normal-prior dispatch, IRLS bypass, ING envelope retention. |
| M | [R/simulationpipeline.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R) | **CMB / numerical stability** — CMB callback closures, sorted log weights; base R, native builders. |
| M | [R/summary.glmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.glmb.R) | **diagnostics** — Summary carries diagnostics; printer uses new S3 method. |
| M | [R/summary.rglmb.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/summary.rglmb.R) | **CMB / diagnostics** — CMB fitted means/MLE limitations and diagnostics; cmb_fitted. |
| M | [README.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/README.md) | **CMB / validation** — Public usage/scope and validation pointers; documentation. |
| A | [inst/examples/Ex_CMB.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/examples/Ex_CMB.R) | **CMB / validation** — Generated-data public example; augmented matrices and dNormal. |
| A | [inst/validation/CMB_V1.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/CMB_V1.md) | **validation** — Release scope; documentation. |
| A | [inst/validation/CMB_VALIDATION.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/CMB_VALIDATION.md) | **validation** — Public suite and integration protocol; testthat/mvtnorm/Python. |
| A | [inst/validation/RNG_REPRODUCIBILITY.md](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/RNG_REPRODUCIBILITY.md) | **shared sampling / validation** — CPU RNG contract and limitations; documentation. |
| A | [inst/validation/cmb-matrix.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/cmb-matrix.R) | **CMB / validation** — Six generated cases, quadrature/importance references and fits; mvtnorm fallback. |
| A | [inst/validation/run-cmb-matrix.py](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/inst/validation/run-cmb-matrix.py) | **validation** — Sequential case timeout/CSV runner; Python standard library and Rscript. |
| A | [man/cmb.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/cmb.Rd) | **CMB** — GENERATED roxygen help for four CMB public helpers. |
| A | [man/print.glmb_diagnostics.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/print.glmb_diagnostics.Rd) | **diagnostics** — GENERATED roxygen S3 printer help. |
| M | [man/rglmb.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/rglmb.Rd) | **diagnostics / shared sampling** — GENERATED roxygen return/RNG help. |
| M | [man/rlmb.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/rlmb.Rd) | **diagnostics** — GENERATED roxygen diagnostic return help. |
| M | [man/simfuncs.Rd](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/man/simfuncs.Rd) | **diagnostics** — GENERATED roxygen ING envelope return help. |
| M | [src/EnvelopeBuild.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/EnvelopeBuild.cpp) | **shared sampling / numerical stability** — Refine fixed tangencies, synchronize G3, add log weights/diagnostics; EnvelopeEval/refine/EnvelopeSort. |
| M | [src/EnvelopeEval.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/EnvelopeEval.cpp) | **CMB** — CPU CMB joint likelihood/gradient dispatch; famfuncs_cmb. |
| M | [src/Envelopefuncs.h](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/Envelopefuncs.h) | **shared sampling** — Append defaulted refinement arguments to builder declaration. |
| D | [src/Makevars](https://github.com/jbogomolovas2/glmbayes/commit/b971596511c90f33de290497b2c7e8feb1a4e580#diff-57dba0e6774ac3caeb9123b3600bddbcfe0e095d1f5b2f2d7115e5e739b0ad2e) | **validation / packaging** — DELETED machine-specific generated build configuration; configure creates environment-appropriate file. |
| M | [src/RcppExports.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/RcppExports.cpp) | **shared sampling** — GENERATED Rcpp trampoline/registration arity 18; Rcpp/export_wrappers. |
| A | [src/envelope_refine.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/envelope_refine.cpp) | **shared sampling / numerical stability** — New refinement, stable means/mass and log categorical selection; Rmath, Armadillo, EnvelopeEval. |
| A | [src/envelope_refine.h](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/envelope_refine.h) | **shared sampling / numerical stability** — New internal declarations, defaults and rationale. |
| M | [src/export_wrappers.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/export_wrappers.cpp) | **shared sampling** — Handwritten Rcpp export wrapper forwards refinement arguments. |
| M | [src/famfuncs.h](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/famfuncs.h) | **CMB** — New internal family kernel declarations; Rcpp/Armadillo/RcppParallel. |
| A | [src/famfuncs_cmb.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/famfuncs_cmb.cpp) | **CMB** — New finite-support normalizer and likelihood/gradient kernels; Rmath/Armadillo. |
| M | [src/rIndepNormalGammaReg.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rIndepNormalGammaReg.cpp) | **shared sampling / diagnostics** — Row streams, interval-mass-aware normal draws, returned envelope. |
| M | [src/rNormalGLM.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalGLM.cpp) | **CMB / shared sampling** — CMB accept score, refinement, cap failure, log selection, row streams and logU. |
| M | [src/rNormalReg.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalReg.cpp) | **diagnostics** — Direct Gaussian sampler returns Envelope=NULL. |
| M | [src/rng_utils.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rng_utils.cpp) | **shared sampling** — Job seed from R, per-row stream assignment, explicit uniform conversion. |
| M | [src/rng_utils.h](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rng_utils.h) | **shared sampling / numerical stability** — Seed function and log_mass overload declarations. |
| M | [src/rnorm_ct.cpp](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rnorm_ct.cpp) | **numerical stability** — Endpoint-preserving tail inverse-CDF draw; shared RNG/Rmath. |
| A | [tests/testthat/test-cmb-public.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-cmb-public.R) | **validation** — Finite-sum likelihood/gradient/binomial identity and fitted means. Depends on testthat and built package. |
| A | [tests/testthat/test-envelope-refinement.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-envelope-refinement.R) | **validation** — Stiff unsplit, already-fixed and large-p fallback examples. Depends on testthat and built package. |
| A | [tests/testthat/test-fit-diagnostics.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-fit-diagnostics.R) | **validation** — Fit/build/summary diagnostics and failure conditions. Depends on testthat and built package. |
| A | [tests/testthat/test-rejection-cap.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-rejection-cap.R) | **validation** — No force acceptance or partial samples on cap. Depends on testthat and built package. |
| A | [tests/testthat/test-rng-reproducibility.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-rng-reproducibility.R) | **validation** — CPU row-stream repeatability and serial/parallel equivalence. Depends on testthat and built package. |
| A | [tests/testthat/test-truncated-normal-tails.R](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/tests/testthat/test-truncated-normal-tails.R) | **validation** — Extreme-tail endpoint and interval behavior. Depends on testthat and built package. |

</details>

## 4. CMB implementation

### Likelihood, parameterization, augmentation and dispatch

For count `yc_i` in `{0,…,m_i}`, the implementation uses

`log p(yc_i | theta_i, nu_i) = theta_i*yc_i + nu_i*log choose(m_i,yc_i) - log sum_{j=0}^{m_i} exp(theta_i*j + nu_i*log choose(m_i,j))`.

Here `theta = X beta + offset_theta`, `nu = Z gamma + offset_nu`, both with identity links. `nu` is unrestricted (including negative values) on finite support. At `nu=1`, `theta` is the binomial logit; for other `nu`, applying `plogis(theta)` is not the mean. The normalizer is a stabilized finite log-sum-exp. The log likelihood is concave in the natural parameters because the normalizer is convex; linear predictors preserve this. A proper positive-definite normal prior supplies posterior curvature, including when the likelihood is weakly identified. All `m=1` inputs are rejected because the second sufficient statistic vanishes. [R family/helpers](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/pfamily.R), [native kernels](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/famfuncs_cmb.cpp).

`cmb_augment` constructs `x = rbind(cbind(X,0),cbind(0,Z))`, defaults `Z` to an intercept, sets `y=c(yc/m,0_N)`, `weights=c(m,1_N)`, and zero offsets of length 2N. The final N response/weight entries are placeholders, not observations. Both predictor functionals are rows of the design; they survive the sampler's right-multiplication/rotation of `x`. Treating `nu` as a single fixed coordinate after rotation would be wrong. Explicit offsets can be supplied to `rglmb` after augmentation; the fitted helper defect below affects their reporting.

For `T=(Y,log choose(m,Y))`, the score is the observed statistic minus its finite-sum expectation. The negative posterior gradient is `-x'*(observed-expected)+P*(b-mu)`. R `f7` constructs information from `Var(Y)`, `Var(log choose(m,Y))` and their covariance, with cross-block terms; it is likelihood information without the prior. `f2/f3` drive BFGS mode finding and `optim`'s numerical Hessian drives standardization. Native `f2_f3_cmb` jointly computes likelihood and gradient; native acceptance scoring uses `f2_cmb` or `f2_cmb_rmat`. There is no inner normalizing solve; work scales with the sum of the finite support sizes, number of evaluated parameter points and repeated evaluations. Tables are built per call; a cheap closed-form binomial kernel is not substituted.

`f4` computes total deviance using a saturated likelihood at the draw's own `nu`; interior observations use bracketed Newton/bisection with tolerance 1e-10 and at most 200 iterations, endpoints use limiting saturated log likelihood zero. CMB dispatch is CPU-only here: `EnvelopeEval` adds CMB to `f2_f3_non_opencl`, and serial/parallel Normal-GLM acceptance adds CMB scoring. `rNormalGLM` forces dispersion to 1 for CMB so trial-count weights are not rescaled. [Callbacks](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/R/simulationpipeline.R), [CPU dispatch](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/EnvelopeEval.cpp), [sampler](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/rNormalGLM.cpp).

Limitations: no complete GLM family apparatus, CMB formula interface, IRLS fit, ordinary inverse link, generic influence/residual method, or CMB OpenCL certification. Use `cmb_fitted`/`cmb_nu` for zero-offset fitted summaries; generic inherited methods need individual assessment. Validation covers small 2–4 coefficient CMB examples with moderate trial sizes, not high-dimensional or very large-support performance.

### Known defects

1. **D1 — native wrapper mismatch/crash:** described in section 1; public non-Gaussian `EnvelopeBuild` is broken. The regression tests use the generated export and missed this path.
2. **D2 — offsets omitted by fitted helpers:** `cmb_fitted` and `cmb_nu` use only `x %*% coefficients`, ignoring fitted `offset`/`offset2`. The supplied synthetic check with theta offset 1 and nu offset .5 returns count mean 2 instead of 2.710088, and nu 1 instead of 1.5. Affects CMB summaries using those helpers. [Reproducer](check-defects.R), [result](defects.log).
3. **D3 — insufficient count validation:** `cmb_augment` accepts fractional counts/trial sizes (e.g. 1.2 of 4.4); downstream R/C++ code rounds implied counts and m, silently changing the stated data. It also does not comprehensively reject zero/nonfinite trial sizes at the public boundary. The finite-support model requires integer m≥1 and integer counts in range. Generated validation inputs satisfy this.
4. **D4 — validation dependency metadata:** `cmb-matrix.R` calls `mvtnorm` for difficult references, but DESCRIPTION does not declare it. It was present here (version recorded in reference metadata); a reproduction needs it explicitly.
5. **Documentation qualifications:** comments saying failed refinement retains the original tangency are too broad: a rejected update retains the last accepted point; the algorithm does not roll the whole cell back to initialization. Nested diagnostic order differs from sorted envelope order. Header claims that refinement can “never” affect correctness assume exact valid evaluation; floating-point/inverse-CDF limitations remain. The unchanged A08 vignette transcription omits square roots in the omega formula and presents Theorem 3 as an inequality; the paper's displayed Theorem 3 is a limit. This is inherited documentation, not a new contribution defect. The C++ initialization matches the paper's square-root formula.

### Original initialization and refinement algorithm

For standardized axis i with curvature `a_i=diag(A)_i`, initialization is

`omega_i = (sqrt(2)-exp(-1.20491-.7321*sqrt(.5+a_i)))/sqrt(1+a_i)`.

Tangencies start at `bStar_i + {-omega_i,0,omega_i}`, and cutpoints at `bStar_i ± omega_i/2`. Gridtype 1 selects one/three points via a curvature threshold, 2 uses `EnvelopeOpt`'s cost criterion with `n_envopt`, 3 always splits, 4 leaves axes unsplit. Products form cells; unsplit axes have region code 4. Refinement leaves these cutpoints, grid size and standardization fixed. [Builder](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/EnvelopeBuild.cpp)

Writing `F(r)` for the negative residual log likelihood including the shifted prior quadratic and `c=grad F(r)`, the restricted proposal is `N(-c,I)` within the fixed cell. Its unnormalized log envelope mass is

`log W_A(r) = -F(r) + r'c + .5*c'c + log Pr[N(-c,I) in A]`.

The target equation is `r=m(r)`, the restricted proposal mean. The code's differentiation gives `grad log W_A = Hessian(log likelihood)*(m-r)`. Under its concavity assumptions, the raw `m-r` direction is a descent direction. Singular curvature means stationarity need not imply this particular fixed point; neither this calculation nor the implementation establishes global uniqueness/convergence.

Each pass computes target means using continued-fraction Mills ratios in positive tails and reflection on the left. For at most 16 parameters and `p*p*cells <= 2,000,000`, it forms forward finite-difference Jacobians of `r-m(r)` (step `1e-5*max(1,abs(r_axis))`) and solves per-cell Newton systems. Failed/nonfinite solves retain the fixed-point proposal. Larger problems use damped fixed-point steps without Jacobian storage. Candidates are projected into their original cell; nonfinite coordinates stay at their current value.

A candidate must have finite evaluated likelihood, gradient and box mass. Its cell log mass must not increase, except a narrowly defined Newton root correction: old residual <1e-3, new residual less than half, and mass increase within `16*epsilon*(1+abs(F)+abs(linear)+abs(quadratic)+abs(box))`. Rejection halves that cell's damping; accepted Newton steps can double it back toward 1. Last accepted values are retained. An unusable initial total mass (NaN/+Inf) returns immediately without moves or full diagnostics; it does not guarantee a usable envelope afterward.

Stopping: maximum residual <1e-6, no representable move, or 60 evaluated passes in this builder. The helper has a `min_gain` default 1e-3, but **the builder passes 0**, disabling the nonnegative-gain early-exit rule (`gain < min_gain`). Cost includes initial evaluation, up to p extra full-grid gradient evaluations per Newton pass, one trial evaluation per pass, per-cell dense solves and bounded Jacobian storage; backtracking can consume further passes. Construction and total time must therefore accompany acceptance rates. [Implementation](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/src/envelope_refine.cpp)

### Relation to the paper

The fixed-cell mean/tangency condition in [Nygren & Nygren (2006)](https://doi.org/10.1198/016214506000000357), section 3.1, motivates the refinement target. The paper does not analyze this particular Newton implementation.

Under the subgradient assumptions, moving a tangency preserves the supporting bound in exact arithmetic. That is separate from how efficiently the sampler runs and how the construction behaves asymptotically. Sections 3.2–3.3 establish limits for the stated normal-model constructions and discuss extension through normal approximation; preservation of those guarantees under this refinement has not been established.

The comparisons below measure acceptance, construction cost and sampling time on specific examples. They also provide data for examining non-normality and correlation, but do not yet isolate which causes the difficult case.

### Controlled comparison protocol

[compare.R](compare.R) toggles `refine=TRUE/FALSE` through the generated envelope builder, then uses the same CPU sampler. The tail, log-weight and RNG safeguards stay enabled. An intercept case was checked against public `rglmb`: proposal counts were identical and coefficient differences were at most 4.89e-12. [Check results](harness-check.log)

Each of ten cases has one warm-up and three measured runs per setting: **80 fits**, requesting **20,000 draws** each. Runs use paired seeds (9100 for warm-up, 9101–9103 for measurements), alternating setting order, fresh R processes, serial CPU sampling and single-threaded numerical libraries. The limits are 200,000 proposals per draw and 300 seconds per fit. `Gridtype=2` and `n_envopt=1000` are fixed, except the single-cell stiff-axis test uses `Gridtype=4`.

The [datasets and model specifications](cases/README.md) describe all six CMB cases. Other inputs are the public Menarche binomial, Dobson Poisson and plant-weight Gaussian examples. The Gaussian comparison uses a fixed plug-in dispersion and forces the envelope path; ordinary Gaussian fits use direct sampling. The stiff-axis example deliberately displaces the starting tangency, so it measures sensitivity to poor placement.

Construction time includes mode finding, standardization and envelope building. Total time also includes sampling, coefficient transformation and saving intermediate results; it excludes R startup, data setup and reference integration. Disabled-refinement residuals are calculated afterward with zero refinement passes, checking that tangencies and weights stay unchanged.

### Results and tradeoffs

Values are **median [minimum, maximum]** over the three measured runs. All completed runs returned 20,000 draws; none timed out. Capped runs report time to failure, not time to complete the requested sample. Their errors establish zero accepted draws from 200,000 proposals; inaccessible native counts remain NA in the raw table. [Failure counts](results/failure-counts.csv) · [Full results](results/summary.csv)

| Case | Refinement | Complete/cap/timeout | Construction s | Sampling s | Total s | Proposals | log mass | Residual |
|---|---|---|---|---|---|---|---|---|
| intercepts | TRUE | 3/0/0 | 0.031 [0.03, 0.032] | 0.306 [0.306, 0.309] | 0.353 [0.351, 0.356] | 25513 [25395, 25540] | -128.55 [-128.55, -128.55] | 4.63989e-07 [4.63989e-07, 4.63989e-07] |
| intercepts | FALSE | 3/0/0 | 0.029 [0.029, 0.029] | 0.308 [0.307, 0.308] | 0.352 [0.351, 0.353] | 25768 [25710, 25835] | -128.541 [-128.541, -128.541] | 0.00614947 [0.00614947, 0.00614947] |
| theta_slope | TRUE | 3/0/0 | 0.035 [0.035, 0.036] | 0.349 [0.34, 0.354] | 0.406 [0.396, 0.414] | 28838 [28702, 28848] | -128.253 [-128.253, -128.253] | 2.31869e-10 [2.31869e-10, 2.31869e-10] |
| theta_slope | FALSE | 3/0/0 | 0.033 [0.032, 0.034] | 0.35 [0.343, 0.351] | 0.405 [0.397, 0.407] | 29252 [29230, 29439] | -128.238 [-128.238, -128.238] | 0.011722 [0.011722, 0.011722] |
| nu_slope | TRUE | 3/0/0 | 0.034 [0.033, 0.035] | 0.338 [0.338, 0.343] | 0.393 [0.391, 0.398] | 28877 [28737, 29037] | -133.582 [-133.582, -133.582] | 3.46712e-10 [3.46712e-10, 3.46712e-10] |
| nu_slope | FALSE | 3/0/0 | 0.033 [0.031, 0.033] | 0.342 [0.34, 0.342] | 0.396 [0.392, 0.396] | 29270 [29250, 29582] | -133.566 [-133.566, -133.566] | 0.00992152 [0.00992152, 0.00992152] |
| strong_prior | TRUE | 3/0/0 | 0.04 [0.04, 0.042] | 0.369 [0.369, 0.374] | 0.436 [0.436, 0.444] | 31951 [31789, 32129] | -124.475 [-124.475, -124.475] | 1.2255e-09 [1.2255e-09, 1.2255e-09] |
| strong_prior | FALSE | 3/0/0 | 0.034 [0.033, 0.034] | 0.374 [0.374, 0.383] | 0.436 [0.436, 0.443] | 32539 [32266, 32659] | -124.457 [-124.457, -124.457] | 0.118226 [0.118226, 0.118226] |
| weak_prior | TRUE | 3/0/0 | 0.042 [0.041, 0.045] | 0.387 [0.38, 0.454] | 0.459 [0.451, 0.524] | 33208 [33086, 33319] | -136.903 [-136.903, -136.903] | 7.71018e-08 [7.71018e-08, 7.71018e-08] |
| weak_prior | FALSE | 3/0/0 | 0.036 [0.035, 0.037] | 0.39 [0.386, 0.392] | 0.455 [0.449, 0.457] | 34130 [34070, 34163] | -136.871 [-136.871, -136.871] | 0.00641066 [0.00641066, 0.00641066] |
| zero_heavy | TRUE | 3/0/0 | 0.129 [0.129, 0.13] | 0.437 [0.436, 0.468] | 0.594 [0.593, 0.628] | 45569 [45388, 45644] | -49.7923 [-49.7923, -49.7923] | 4.45723e-08 [4.45723e-08, 4.45723e-08] |
| zero_heavy | FALSE | 0/3/0 | 0.036 [0.036, 0.037] | 1.217 [1.215, 1.222] | 1.275 [1.275, 1.28] | NA | 751131 [751131, 751131] | 929.771 [929.771, 929.771] |
| stiff_axis | TRUE | 3/0/0 | 0.003 [0.002, 0.003] | 4.642 [4.626, 4.819] | 4.671 [4.651, 4.844] | 2.10052e+06 [2.09023e+06, 2.10521e+06] | -0.918939 [-0.918939, -0.918939] | 1.04592e-08 [1.04592e-08, 1.04592e-08] |
| stiff_axis | FALSE | 0/3/0 | 0.003 [0.002, 0.003] | 0.431 [0.43, 0.438] | 0.437 [0.435, 0.443] | NA | 500554 [500554, 500554] | 1001 [1001, 1001] |
| gaussian | TRUE | 3/0/0 | 0.003 [0.003, 0.003] | 0.213 [0.211, 0.214] | 0.234 [0.233, 0.235] | 25441 [25366, 25513] | -25.2209 [-25.2209, -25.2209] | 3.25469e-08 [3.25469e-08, 3.25469e-08] |
| gaussian | FALSE | 3/0/0 | 0.003 [0.002, 0.003] | 0.214 [0.212, 0.214] | 0.234 [0.234, 0.235] | 25441 [25367, 25514] | -25.2209 [-25.2209, -25.2209] | 0.000113381 [0.000113381, 0.000113381] |
| binomial | TRUE | 3/0/0 | 0.003 [0.003, 0.003] | 0.226 [0.224, 0.229] | 0.246 [0.246, 0.251] | 25473 [25416, 25539] | -62.6083 [-62.6083, -62.6083] | 5.44887e-09 [5.44887e-09, 5.44887e-09] |
| binomial | FALSE | 3/0/0 | 0.003 [0.002, 0.003] | 0.23 [0.229, 0.234] | 0.251 [0.251, 0.256] | 25491 [25408, 25539] | -62.6077 [-62.6077, -62.6077] | 0.00268295 [0.00268295, 0.00268295] |
| poisson | TRUE | 3/0/0 | 0.017 [0.017, 0.019] | 0.29 [0.289, 0.296] | 0.344 [0.343, 0.352] | 36843 [36396, 36855] | -36.0385 [-36.0385, -36.0385] | 8.8669e-11 [8.8669e-11, 8.8669e-11] |
| poisson | FALSE | 3/0/0 | 0.003 [0.003, 0.004] | 0.29 [0.287, 0.291] | 0.33 [0.327, 0.331] | 36941 [36770, 37087] | -36.0326 [-36.0326, -36.0326] | 0.0186753 [0.0186753, 0.0186753] |


The **zero-heavy CMB case** finishes only with refinement under these settings: median 0.594 seconds and 45,569 proposals, or about 44% acceptance. Refinement takes 31 passes over 81 cells. All three disabled runs reach the cap before their first accepted draw.

The **stiff-axis stress test** also finishes only with refinement, but still needs about 105 proposals per draw. Improving the starting tangency helps; the single unsplit cell remains inefficient.

The other five CMB cases finish with either setting. Their 0.9–3.2% expected proposal reductions produce little change in total time. Gaussian is essentially unchanged, and Menarche has a small timing difference. For Poisson, extra construction work outweighs the acceptance gain: median total time increases from 0.330 to 0.344 seconds. Three short repetitions are insufficient for a general speed claim.

### Regression, references, environment and reproduction

The public suite recorded **302 passing assertions, no failures or test errors, and 11 OpenCL skips** across 68 tests. The wrapper crash was reproduced separately and remains a known defect. [Test log](regression.log) · [Table](regression.csv) · [R results](regression.rds)

For five CMB cases, two scaled Gauss-Hermite integrations agree within 1e-5 in log mass, means and covariance. The zero-heavy case uses the existing Student-t importance-sampling fallback: two seeds, effective sample sizes around 579,000, and agreement within 1.147 combined standard errors. All six references resolved. [Reference status](references/status.csv)

All **33 completed measured CMB runs** pass the six-standard-error checks for means, covariance and proposal rates, including reference uncertainty. The three capped zero-heavy runs are failures. Integration shares the compiled likelihood; separate finite-sum and derivative tests check the likelihood itself. [Per-run checks](results/reference-checks.csv)

The recorded environment is R 4.5.2 on arm64 macOS Ventura, glmbayes 0.9.76.9001, Apple clang 14.0.3 and Accelerate BLAS. [Full environment](environment.txt) · [Reproduction instructions](README.md#inspect-or-reproduce) · [Checksums](SHA256SUMS)

The [later upstream changes](UPSTREAM_CHANGES.md) were not included in these measurements and need an integration review. The reviewed contribution still enables refinement by default; this documentation branch does not change package behavior.
