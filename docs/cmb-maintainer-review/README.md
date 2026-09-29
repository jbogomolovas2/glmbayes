# CMB contribution: maintainer review

This is the review requested in [glmbayes issue #61](https://github.com/knygren/glmbayes/issues/61). It describes contribution **`b971596511c90f33de290497b2c7e8feb1a4e580`** against upstream base **`7959a42d32f822a4b65ddfaf94927190ad5db536`**. This branch adds review material only; the reviewed package code and sampler defaults are unchanged.

## Start here

1. **[Four-part overview and detailed inventory](MAINTAINER_REVIEW.md#review-at-a-glance):** changed signatures, additional returned items, new functions/files, and the CMB implementation.
2. **[Evidence summary](EVIDENCE_SUMMARY.md):** controlled refinement comparisons, regression results, defects, and limits of the evidence.
3. **[Timing and sampling table](results/summary.md):** medians and ranges for the ten comparison cases. [CSV](results/summary.csv), [individual runs](results/runs.csv), and [reference checks](results/reference-checks.csv) support the summaries.
4. **[Datasets and exact model specifications](cases/README.md):** six downloadable synthetic datasets, priors, predictor formulas, and the paired zero-heavy settings.
5. **[Later upstream changes](UPSTREAM_CHANGES.md):** changes observed on September 28, 2026, outside the measured binary.

The strongest CMB result is the generated zero-heavy case: all three measured 20,000-draw runs completed with refinement (median 0.594 seconds); all three disabled runs hit the 200,000-proposal cap before accepting their first draw. The other five CMB cases completed either way, with modest proposal improvements. Poisson was slightly slower with refinement. The stiff-axis example is an artificial placement stress test.

All 33 completed measured CMB fits passed the existing integration tolerances. The public regression suite recorded 302 passing assertions and 11 OpenCL skips. Capped runs are failures. These CPU results do not establish universal speed improvements, OpenCL behavior, or preservation of the paper's asymptotic guarantees.

**Known defects remain in the reviewed code:** the public envelope wrapper has a 16/18-argument mismatch, CMB fitted helpers omit offsets, count validation permits fractional inputs, and the validation dependency `mvtnorm` is undeclared. They are documented here for subsequent correction, not silently fixed in this snapshot.

## Inspect or reproduce

Git, Python 3 (standard library only), R, and the compiler/toolchain needed to install this R package are required. Install the pinned package's `Depends`, `Imports`, and `LinkingTo` dependencies, plus `testthat` and `mvtnorm`, in your normal R library first. The [pinned DESCRIPTION](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/DESCRIPTION) lists package/build dependencies; [environment.txt](environment.txt) and the reference logs record the versions used for the saved results. Optional test dependencies and OpenCL availability can change skips on another machine.

```sh
git clone --branch maintainer-review https://github.com/jbogomolovas2/glmbayes.git
cd glmbayes/docs/cmb-maintainer-review
python3 verify-package.py
sh reproduction.sh --check
sh reproduction.sh
```

The check command verifies dependencies and the pinned Git revision without running fits. The final command creates a fresh temporary run directory, prints its location, archives the pinned source, builds a separate library, and runs the public suite, 80 serial comparisons, six reference integrations, and the summaries. Alternatively, pass a **new, nonexistent output directory** as the single argument. Saved evidence in this repository is never overwritten. The runner stops on build/test/reference errors; the two documented capped cases remain explicitly recorded comparison outcomes. Each fit has a 300-second timeout; each reference has a separate 30-minute timeout. Reference integration can take substantially longer than the fits.

The optional `reproduce-wrapper-defect.R` deliberately reproduces an R-process crash and is excluded from the normal runner. Individual `.rds` files may be downloaded from GitHub and opened with `readRDS()` for reanalysis.

[Publication notes](PUBLICATION_NOTES.md) explain the publication edits and path redactions. [SHA256SUMS](SHA256SUMS) covers every published file except the manifest itself. The complete evidence is available as individual files; no duplicate ZIP, private study data, supplied paper, or installed build library is committed here.
