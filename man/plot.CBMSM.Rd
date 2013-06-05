\name{plot.CBMSM}
\alias{plot.CBMSM}
\title{Plotting Covariate Balancing Propensity Score Estimation}
\description{Plots the absolute difference in standardized means before and after weighting.}
\usage{\method{plot}{CBMSM}(x,stabilized=TRUE,...)}
\arguments{
  \item{x}{an object of class \dQuote{CBMSM}.}
  \item{stabilized}{Set to \code{TRUE} to use stabilized weights, and \code{FALSE} to use unstabilized weights.  Default is \code{TRUE}.}
  \item{...}{Additional arguments to be passed to plot.}
}
\details{Covariate balance is improved if the plot's points are below the plotted line of y=x..}
\value{
For each contrast and all covariates, plots the absolute mean difference before weighting on the x-axis and after weighting on the y-axis.
}

\author{Marc Ratkovic}

\seealso{\link{CBPS}, \link{plot}}