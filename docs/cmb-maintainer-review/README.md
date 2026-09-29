# CMB contribution: maintainer review

This review answers Kjell's questions in [issue #61](https://github.com/knygren/glmbayes/issues/61).

- **[Four-part review](MAINTAINER_REVIEW.md#review-at-a-glance):** changed arguments, returned fields, new functions/files, and the CMB implementation.
- **[Results in brief](EVIDENCE_SUMMARY.md):** where refinement helped, its cost, and the remaining problems.
- **[Datasets and models](cases/README.md):** all six synthetic CMB datasets, with predictors, priors and comparison settings.
- **[Timing table](results/summary.md):** results for each case, with [individual runs](results/runs.csv) and [reference checks](results/reference-checks.csv).

The main finding is that the zero-heavy CMB example completed with refinement but stalled without it. The other CMB examples worked either way. The review also identifies several bugs to fix before integration. Whether refinement preserves the paper's asymptotic guarantees remains open.

The reviewed contribution is `b971596511c90f33de290497b2c7e8feb1a4e580`, based on `7959a42d32f822a4b65ddfaf94927190ad5db536`. This branch adds documentation and evidence; package code is unchanged. [Later upstream changes](UPSTREAM_CHANGES.md) are listed separately.

## Inspect or reproduce

You need Git, Python 3, R and the tools to compile the package. Install the dependencies in the [reviewed DESCRIPTION](https://github.com/jbogomolovas2/glmbayes/blob/b971596511c90f33de290497b2c7e8feb1a4e580/DESCRIPTION), including `LinkingTo`, plus `testthat` and `mvtnorm`. [environment.txt](environment.txt) records the versions used here.

```sh
git clone --branch maintainer-review https://github.com/jbogomolovas2/glmbayes.git
cd glmbayes/docs/cmb-maintainer-review
python3 verify-package.py
sh reproduction.sh --check
sh reproduction.sh
```

The first two checks verify saved files and required dependencies. The final command builds the reviewed version in a fresh temporary directory and runs the tests, comparisons and reference calculations. It prints the output location. You can instead pass a new output-directory path as its argument; existing directories are refused.

Each fit has a five-minute timeout and each reference calculation has a 30-minute timeout. Reference calculations can take longer than the fits. Skips may differ with installed packages and OpenCL availability.

Saved draws, envelopes and references are provided as `.rds` files for use with `readRDS()`. The optional `reproduce-wrapper-defect.R` deliberately reproduces an R crash and is excluded from the normal runner.

[Publication notes](PUBLICATION_NOTES.md) · [File checksums](SHA256SUMS)
