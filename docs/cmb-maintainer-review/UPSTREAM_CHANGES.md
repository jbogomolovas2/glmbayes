# Subsequent upstream changes (outside this review binary)

Compared base `7959a42d32f822a4b65ddfaf94927190ad5db536` to observed upstream HEAD `dc2d905d2244350c76068609abe031bb01421132`. GitHub API reports ahead_by=8, behind_by=0, merge base equal to the recorded base. Local repository is shallow; raw local reachability counts include imported older history and are not the comparison commit count. Net snapshot diff: 39 files, 2,058 insertions, 1,123 deletions.

These are revision metadata, not contribution source links. No rebasing or integration was performed.

| Commit | Subject |
|---|---|
| `b24df6e57dd7a00cf564542f3c4d8d305fc71eb3` | Delegate get_opencl_core_count to opencltools::get_opencl_core_count() |
| `2e66bc150f3876673a6076da59703be875ee601c` | Add ING prior guard |
| `5477dc3efff8f3c86159844faadf995381ba5a69` | Gate cat()/print() output per glmbayesCore  CRAN reviewer feedback. Also return values for methods. |
| `9e977d7c3b6d566745ba97bd066893bc7980fff3` | Update progress bar utilities in line with glmbayeCore |
| `b636dd86369b3f22dfbe500abb809838c2b9e091` | Update dic_info, residuals, and summary functions to make them consistent with glmbayescore |
| `5fda1d3e821db825246dece1f881866eb85ccae8` | Add glmbayesCore to Imports and Prior_Setup to importFrom |
| `603357bb981ad36d407b41cc7f2a799af776b58b` | Merge origin/main; reconcile NEWS and Stage 0 migration docs |
| `dc2d905d2244350c76068609abe031bb01421132` | UPdate glmbayes |

| Status | Path |
|---|---|
| modified | `DESCRIPTION` |
| modified | `NAMESPACE` |
| modified | `NEWS.md` |
| added | `R/dic_info.R` |
| modified | `R/glmb.R` |
| modified | `R/gpu_diagnostics.R` |
| added | `R/ing_prior_guard.R` |
| modified | `R/prior.R` |
| modified | `R/reexports.R` |
| modified | `R/residuals.glmb.R` |
| modified | `R/simfunction.R` |
| modified | `R/summary.rglmb.R` |
| modified | `data-raw/CORE_MIGRATION_DIFF.md` |
| added | `data-raw/CPP_R_CALLBACK_INVENTORY.md` |
| added | `data-raw/_diff_paths.txt` |
| added | `data-raw/_gb_only_paths.txt` |
| added | `data-raw/_gc_only_paths.txt` |
| added | `data-raw/_inventory_core_migration.R` |
| added | `data-raw/cpp_r_callback_inventory.R` |
| modified | `man/Prior_Setup.Rd` |
| modified | `man/gpu_diagnostics.Rd` |
| modified | `man/residuals.glmb.Rd` |
| modified | `man/summary.rglmb.Rd` |
| modified | `src/EnvelopeBuild.cpp` |
| modified | `src/EnvelopeBuild_Ind_Normal_Gamma.cpp` |
| modified | `src/EnvelopeOrchestrator.cpp` |
| modified | `src/EnvelopeSize.cpp` |
| modified | `src/Envelopefuncs.h` |
| modified | `src/R_interface.h` |
| modified | `src/kernel_loader.cpp` |
| modified | `src/kernel_wrappers.cpp` |
| modified | `src/openclPort.h` |
| added | `src/package_ns.h` |
| modified | `src/progress_utils.cpp` |
| modified | `src/progress_utils.h` |
| modified | `src/rGammaGamma.cpp` |
| modified | `src/rGammaGaussian.cpp` |
| modified | `src/rIndepNormalGammaReg.cpp` |
| modified | `src/rNormalGammaReg.cpp` |
