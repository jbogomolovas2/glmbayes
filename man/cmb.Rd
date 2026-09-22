\name{cmb}
\alias{cmb}
\alias{cmb_augment}
\alias{cmb_fitted}
\alias{cmb_nu}
\title{Conway-Maxwell-binomial regression helpers}
\usage{
cmb(link = "identity")
cmb_augment(X, yc, m, Z = NULL)
cmb_fitted(object, scale = c("proportion", "count"))
cmb_nu(object)
}
\arguments{
\item{link}{The natural-parameter identity link.}
\item{X}{N by p numeric design matrix for theta.}
\item{yc}{N integer counts between zero and m.}
\item{m}{Positive integer trial counts, either length N or a scalar.}
\item{Z}{N by q numeric design matrix for nu; NULL gives an intercept.}
\item{object}{A fitted CMB model produced using the augmented design.}
\item{scale}{Return expected proportions or counts.}
}
\description{Construct and summarize CMB models using the normal-prior
\code{rglmb} interface. Probability is proportional to
\eqn{\exp(\theta y + \nu\log {m\choose y})} for y from zero to m.}
\details{Theta and nu are X beta and Z gamma. Nu has an identity link and is
not constrained positive: finite support permits finite real natural parameters.
Nu equal to one gives a binomial with probability \code{plogis(theta)}. In
general that transformation is not the expected proportion; use \code{cmb_fitted}.

The family is a specialized marker, not a complete stats GLM family. Pass the
outputs of \code{cmb_augment} to \code{rglmb}. The design stacks N theta rows
above N nu rows. The trailing response/weight values are placeholders. Generic
GLM methods assuming a scalar mean link may not apply. Random effects are not
supported.}
\value{\code{cmb} returns a family marker. \code{cmb_augment} returns x (2N by
p+q), y, weights and offset. \code{cmb_fitted} and \code{cmb_nu} return matrices
with one row per posterior draw and one column per original observation.
Use column means or quantiles to summarize these draws.}
\examples{
X <- matrix(1, 4, 1)
a <- cmb_augment(X, c(0, 1, 3, 4), m = 5)
dim(a$x)
}
