% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CBPSMain.R
\name{summary.CBPS}
\alias{summary.CBPS}
\title{Summarizing Covariate Balancing Propensity Score Estimation}
\usage{
\method{summary}{CBPS}(object, ...)
}
\arguments{
\item{object}{an object of class \dQuote{CBPS}, usually, a result of a call
to CBPS.}

\item{...}{Additional arguments to be passed to summary.}
}
\value{
\item{call}{The matched call.} \item{deviance.residuals}{The five
number summary and the mean of the deviance residuals.}
\item{coefficients}{A table including the estimate for the each coefficient
and the standard error, z-value, and two-sided p-value for these estimates.}
\item{J}{Hansen's J-Statistic for the fitted model.}
\item{Log-Likelihood}{The log-likelihood of the fitted model.}
}
\description{
Prints a summary of a fitted CBPS object.
}
\details{
Prints a summmary of a CBPS object, in a format similar to glm.  The
variance matrix is calculated from the numerical Hessian at convergence of
CBPS.
}
\seealso{
\link{CBPS}, \link{summary}
}
\author{
Christian Fong, Marc Ratkovic, and Kosuke Imai.
}
