\name{balance}
\alias{balance}
\title{Optimal Covariate Balance}
\description{Returns the mean and standardized mean associated with each treatment group, before and after weighting.}
\usage{balance(object, ...)}
\arguments{
  \item{object}{A CBPS, CBMB, or CBMSM object.}
  \item{...}{Additional arguments to be passed to balance.}
}
\details{Each of the matrices' rows are the covariates and whose columns are the weighted mean, and standardized mean associated with each treatment group.  The standardized mean is the weighted mean divided by the standard deviation of the covariate for the whole population.}
\value{
  Returns a list of two matrices, "original" (before weighting) and "balanced" (after weighting).  
}
\author{Marc Ratkovic, Kosuke Imai, and Christian Fong.}

\seealso{\link{CBPS}}