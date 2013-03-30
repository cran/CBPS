\name{MSMdata}
\alias{MSMdata}
\title{Marginal Structural Model Data Generator}
\description{
  \code{MSMdata} simulates data to be fit with a marginal structural model
}
\usage{
	 MSMdata(n)
}
\arguments{
  \item{n}{The desired sample size.}
}
\details{\code{MSMdata} simulates data to be fit with a marginal structural model}
\value{
  Returns an object with the following:
  \item{treat.1, treat.2, treat.3}{The observed binary treatment at times 1, 2, and 3.}
  \item{X1, X2, X3}{Covariate matrices for times 1, 2, and 3.  The last column is the previous observed treatment vector.}
  \item{probs.joint, probs.uncond}{The true joint probabilities and unconditional for each treatment sequence, of length n.  The true (stabilized) marginal structural weights are the unconditional probabilities divided by the joint probabilities. }
  \item{y}{The outcome, of length n.}
  \item{n.t}{Number of unique time values}
}
\references{Imai, Kosuke and Marc Ratkovic}

\author{Marc Ratkovic and Kosuke Imai.}

\seealso{\link{CBPS}}

\examples{

##Generates a dataset with 100 observations.
data1<-MSMdata(100)
}