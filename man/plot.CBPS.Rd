\name{plot.CBPS}
\alias{plot.CBPS}
\title{Plotting Covariate Balancing Propensity Score Estimation}
\description{Plots the absolute difference in standardized means before and after weighting.}
\usage{\method{plot}{CBPS}(x,...)}
\arguments{
  \item{x}{an object of class \dQuote{CBPS}, usually, a result of a call to CBPS.}
  \item{...}{Additional arguments to be passed to plot.}
}
\details{The "Before Weighting" plot gives the absolute difference in standardized means before weighting, and the "After Weighting" plot gives the after.}
\value{
Plots the absolute difference in standardized means by contrast for all covariates, both before and after weighting.  Prints a mapping of the contrast names given in the axes to the names of the treatments to the console.
}

\author{Marc Ratkovic}

\seealso{\link{CBPS}, \link{plot}}