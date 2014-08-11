\name{plot.CBPS}
\alias{plot.CBPS}
\title{Plotting Covariate Balancing Propensity Score Estimation}
\description{Plots the absolute difference in standardized means before and after weighting.}
\usage{\method{plot}{CBPS}(x, covars = NULL, silent = TRUE,...)}
\arguments{
  \item{x}{an object of class \dQuote{CBPS}, usually, a result of a call to CBPS.}
  \item{covars}{Indices of the covariates to be plotted (excluding the intercept).  For example, if only the first two covariates from \code{balance} are desired, set \code{covars} to 1:2.  The default is \code{NULL}, which plots all covariates.}
  \item{silent}{If set to \code{FALSE}, returns the absolute imbalance for each contrast pair before and after weighting.  This helps the user to create his or her own customized plot. Default is \code{TRUE}, which returns nothing.}
  \item{...}{Additional arguments to be passed to plot.}
}
\details{The "Before Weighting" plot gives the balance before weighting, and the "After Weighting" plot gives the balance after weighting.}
\value{
For binary and multi-valued treatments, plots the absolute difference in standardized means by contrast for all covariates before and after weighting.  For continuous treatments, plots the absolute Pearson correlation between the treatment and each covariate.
}

\author{Christian Fong, Marc Ratkovic, and Kosuke Imai.}

\seealso{\link{CBPS}, \link{plot}}