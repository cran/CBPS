\name{summary.CBPSlogit}
\alias{summary.CBPSlogit}
\title{Summarizing Covariate Balancing Propensity Score for Logistic Regression}
\description{Prints a summary of a fitted CBPSlogit object.}
\usage{\method{summary}{CBPSlogit}(object,...)}
\arguments{
  \item{object}{an object of class \dQuote{CBPSlogit}, usually, a result of a call to CBPSlogit.}
  \item{...}{Additional arguments to be passed to summary.}
}
\details{Prints a summmary of a CBPSlogit object, in a format similar to glm.  The variance matrix is calculated from the numerical Hessian at convergence of CBPSlogit.}
\value{
  \item{call}{The matched call.}
  \item{deviance.residuals}{The five number summary and the mean of the deviance residuals.}
  \item{coefficients}{A table including the estimate for the each coefficient and the standard error, z-value, and two-sided p-value for these estimates.}
  \item{J}{Hansen's J-Statistic for the fitted model.}
}

\author{Marc Ratkovic}

\seealso{\link{CBPSlogit}, \link{summary}}