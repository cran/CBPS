\name{balance.CBMB}
\alias{balance.CBMB}
\title{Optimal Covariate Balance for Multiple Simulatenously Administered Binary Treatments}
\description{Returns the mean and standardized mean associated with each treatment group, before and after weighting.}
\usage{\method{balance}{CBMB}(object, stabilized, ...)}
\arguments{
  \item{object}{A CBMB object.}
   \item{stabilized}{Set to \code{TRUE} to use stabilized weights, and \code{FALSE} to use unstabilized weights.  Default is \code{TRUE}.}
   \item{...}{Additional arguments to be passed to balance.}
}
\details{Each of the matrices' rows are the covariates and whose columns are the weighted mean, and standardized mean associated with each treatment group.  The standardized mean is the weighted mean divided by the standard deviation of the covariate for the whole population.}
\value{
  Returns a list of two matrices, "original" (before weighting) and "balanced" (after weighting).  
}
\author{Christian Fong, Marc Ratkovic, and Kosuke Imai.}

\examples{

\dontrun{
###
### Example: Optimal Covariate Balance
###
n<-200
k<-4
set.seed(1040)
X1<-cbind(1,matrix(rnorm(n*k),ncol=k))

betas.1<-betas.2<-betas.3<-c(2,4,4,-4,3)/5


probs.1<-probs.2<-probs.3<-(1+exp(-X1 \%*\% betas.1))^-1

treat.1<-rbinom(n=length(probs.1),size=1,probs.1)
treat.2<-rbinom(n=length(probs.2),size=1,probs.2)
treat.3<-rbinom(n=length(probs.3),size=1,probs.3)

formulas.multi<-list(c(treat.1~X1, treat.2~X1, treat.3~X1))

multibin1<-CBPS(formulas.multi,method="over",type="MultiBin")
balance(multibin1)
}
}