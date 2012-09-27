\name{IPW}
\alias{IPW}
\title{Inverse Probability Weighting Estimators}
\description{
  \code{IPW} calculates the inverse probability weighting estimators from Lunceford and Davidian (2004).
}
\usage{
	  IPW(outcome, treat, data, pscore, k)
}
\arguments{
  \item{outcome}{A vector of the outcome variable.}
  \item{treat}{A vector of treatment assignments, where 1 indicates assignment to the treatment group and 0 indicates assignment to the control group.}
  \item{data}{an optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{CBPS} is called.}
  \item{pscore}{A vector of propensity scores.}
  \item{k}{The number of covariates considered in the generation of the propensity scores.}
}
\details{Calculates the inverse probability weighting estimators from Lunceford and Davidian (2004).  For use with CBPS, run only with ATE.}
\value{
  Returns a data frame whose columns are the point estimates and standard errors for these inverse probability weighting estimators, and whose rows correspond to IPW1, IPW2, and IPW3.
}
\references{Lunceford, Jared and Marie Davidian.  ``Stratification and Weighting Via the Propensity Score in Estimation of Causal Treatment Effects: A Comparative Study.'' Statistics in Medicine 23, 15 October 2004.}

\author{Marc Ratkovic, Kosuke Imai, and Christian Fong.}

\seealso{\link{DR}}

\examples{
###
### Example: propensity score matching
###

##Load the LaLonde data
data(LaLonde)
## Estimate CBPS via logistic regression for ATE (Run with ATE only)
fit <- CBPS(treat ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0), data = LaLonde, ATT = FALSE)
##Calculate the IPW estimators
inverse.prob.wt <- IPW(outcome = re78, treat = treat, data = LaLonde, pscore = fitted(fit), k = length(coef(fit)))

}