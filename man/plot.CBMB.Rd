\name{plot.CBMB}
\alias{plot.CBMB}
\title{Plotting Covariate Balancing Propensity Score Estimation for Multiple Simultaneously Administered Binary Treatments}
\description{Plots the absolute difference in standardized means before and after weighting.}
\usage{\method{plot}{CBMB}(x,stabilized=TRUE,...)}
\arguments{
  \item{x}{an object of class \dQuote{CBMB}.}
  \item{stabilized}{Set to \code{TRUE} to use stabilized weights, and \code{FALSE} to use unstabilized weights.  Default is \code{TRUE}.}
  \item{...}{Additional arguments to be passed to plot.}
}
\details{The "Before Weighting" plot gives the absolute difference in standardized means before weighting, and the "After Weighting" plot gives the after. By default, uses stabilized weights.}
\value{
Plots the absolute difference in standardized means by contrast for all covariates, both before and after weighting.  Prints a mapping of the contrast names given in the axes to the names of the treatments to the console.
}

\author{Marc Ratkovic}

\seealso{\link{CBPS}, \link{plot}}