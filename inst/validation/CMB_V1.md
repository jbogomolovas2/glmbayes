# CMB v1 scope

This contribution derives from local CMB v1, version 0.9.76.9001, based on
upstream GLMBayes 0.9.76. It is not an official upstream release.

Retained: natural-parameter CMB regression, mean conversion, bounded rejection,
safeguarded refinement, extreme-tail correction, reproducible CPU RNG and fit
diagnostics. CMB uses the documented augmented matrix interface.

Excluded: experimental adaptive splitting, cost controllers, raw study data,
local analysis outputs and backups. No random-effects functionality was added.
See `CMB_VALIDATION.md` and the repository's `CMB_CONTRIBUTION.md`.
