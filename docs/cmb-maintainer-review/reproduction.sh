#!/bin/sh
set -eu
# Run the immutable contribution in a new directory, preserving saved evidence.
PIN=b971596511c90f33de290497b2c7e8feb1a4e580
review_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(git -C "$review_dir" rev-parse --show-toplevel)
git -C "$repo_dir" cat-file -e "$PIN^{commit}"
case "${1-}" in
  --check) check_only=yes ;;
  --*) echo 'Usage: sh reproduction.sh [--check | NEW_OUTPUT_DIRECTORY]' >&2; exit 2 ;;
  *) check_only=no ;;
esac
if [ "$#" -gt 1 ]; then
  echo 'Usage: sh reproduction.sh [--check | NEW_OUTPUT_DIRECTORY]' >&2
  exit 2
fi
command -v python3 >/dev/null
command -v R >/dev/null
command -v Rscript >/dev/null
git -C "$repo_dir" show "$PIN:DESCRIPTION" | Rscript --vanilla -e '
d <- read.dcf(file("stdin"))
fields <- intersect(c("Depends", "Imports", "LinkingTo"), colnames(d))
packages <- trimws(gsub("[[:space:]]*[(][^)]*[)]", "", unlist(strsplit(paste(d[1, fields], collapse=","), ","))))
packages <- setdiff(unique(c(packages, "testthat", "mvtnorm")), c("R", ""))
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing)) stop("Install required R packages first: ", paste(missing, collapse=", "))
cat("Pinned revision and required R packages are available.\n")
'
if [ "$check_only" = yes ]; then exit 0; fi
if [ "$#" -eq 0 ]; then
  run_dir=$(mktemp -d "${TMPDIR:-/tmp}/cmb-review.XXXXXX")
else
  if [ -e "$1" ]; then echo 'Output directory already exists; choose a new directory.' >&2; exit 2; fi
  mkdir -p -- "$1"
  run_dir=$(CDPATH= cd -- "$1" && pwd)
fi
echo "Reproduction directory: $run_dir"
for script in compare.R run-comparisons.py references.R run-references.py validate-results.R summarize.py regression.R reproduce-wrapper-defect.R; do
  cp "$review_dir/$script" "$run_dir/$script"
done
cd "$run_dir"
mkdir source library
git -C "$repo_dir" archive --output="$run_dir/source.tar" "$PIN"
tar -xf source.tar -C source
R CMD INSTALL --preclean --no-multiarch -l library source > install.log 2>&1
export R_LIBS="$run_dir/library${R_LIBS:+:$R_LIBS}" OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 RCPP_PARALLEL_NUM_THREADS=1
Rscript --vanilla regression.R > regression.log 2>&1
python3 run-comparisons.py > comparisons.log 2>&1
python3 run-references.py > references.log 2>&1
Rscript --vanilla validate-results.R > validation.log 2>&1
python3 summarize.py
Rscript --vanilla -e '
r <- read.csv("results/runs.csv")
refs <- read.csv("references/status.csv")
checks <- read.csv("results/reference-checks.csv")
stopifnot(nrow(r)==80L, nrow(refs)==6L, all(refs$status=="resolved"))
expected_cap <- !r$refine & r$case %in% c("zero_heavy", "stiff_axis")
stopifnot(all(r$status[!expected_cap]=="completed"), all(r$status[expected_cap] %in% c("completed", "proposal_cap")))
completed_cmb <- sum(r$status=="completed" & r$repetition>0 & r$case %in% refs$case)
stopifnot(nrow(checks)==completed_cmb, all(checks$passed))
cat("Completed CMB comparisons pass reference checks; capped comparisons remain failures.\n")
'
echo "Completed. Results and logs are in $run_dir"
