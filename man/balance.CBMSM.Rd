\name{balance.CBMSM}
\alias{balance.CBMSM}
\title{Optimal Covariate Balance for Marginal Structural Models}
\description{Returns the mean and standardized mean associated with each treatment group, before and after weighting.}
\usage{\method{balance}{CBMSM}(object, stabilized, ...)}
\arguments{
  \item{object}{A CBMSM object.}
  \item{stabilized}{Set to \code{TRUE} to use stabilized weights, and \code{FALSE} to use unstabilized weights.  Default is \code{TRUE}.}
  \item{...}{Additional arguments to be passed to balance.}
}
\details{Each of the matrices' rows are the covariates and whose columns are the weighted mean, and standardized mean associated with each treatment group.  The standardized mean is the weighted mean divided by the standard deviation of the covariate for the whole population.}
\value{
  Returns a list whose elements are lists associated with each time step in the marginal structural model.   Each of the sublists has two matrices, "original" (before weighting) and "balanced" (after weighting).  
}
\author{Marc Ratkovic, Christian Fong, and Kosuke Imai.}

\examples{

\dontrun{
###
### Example: Optimal Covariate Balance
###
data1<-MSMdata(200)
attach(data1)
formulas.msm<-list(c(treat.1~X1, treat.2~X2, treat.3~X3))
msm1<-CBPS(formulas.msm,method="over",type="MSM")
balance(msm1)
detach(data1)
}
}