\name{print.glmb_diagnostics}
\alias{print.glmb_diagnostics}
\alias{glmbayes_sampling_error}
\title{Print sampling diagnostics}
\usage{\method{print}{glmb_diagnostics}(x, ...)}
\arguments{
\item{x}{A fit's diagnostics object.}
\item{...}{Unused additional arguments.}
}
\description{Fits returned by \code{rglmb}, \code{rlmb}, and their formula interfaces retain sampling information in \code{fit$diagnostics}. Summaries print this information when present.}
\details{
The status is \code{completed} or \code{refinement_incomplete}. Envelope status is \code{converged}, \code{incomplete}, \code{not_reported}, or \code{not_applicable}. Refinement convergence describes tangency optimization, not MCMC convergence or proof of posterior accuracy. A converged envelope can still have poor acceptance.

Proposal totals, proposals per draw, maximum proposals for a draw, and acceptance rate are calculated from the returned iteration counts. They exclude setup, pilot, and calibration work. Missing or invalid counts are \code{NA}; direct sampling has no envelope acceptance rate. Envelope cell count, refinement passes, and final residual are retained when reported by the backend.

Requested parallel and OpenCL flags record the request, not proof that a GPU was used. Build information includes package version, loaded library path, R version, platform, and MD5 fingerprints of the loaded compiled library and installed R database. Fingerprints identify build files and are not security guarantees. The R RNG kind is recorded without advancing its state.

A rejection proposal limit raises a \code{glmbayes_sampling_error} condition with a \code{diagnostics} field containing the limit, failure reason, and suggested action. No partial draws are returned. Other errors retain their original classes. This condition does not bound elapsed execution time.
}
\value{The diagnostics object, invisibly.}
\examples{
\dontrun{
print(fit$diagnostics)
tryCatch(rglmb(...), glmbayes_sampling_error = function(e) e$diagnostics)
}
}
