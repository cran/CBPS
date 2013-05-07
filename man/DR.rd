\name{DR}
\alias{DR}
\title{Doubly Robust Estimator}
\description{
  \code{DR} calculates the doubly-robust estimator as as described in Lunceford and Davidian (2004) and originally detailed in Robins, Rotnitzky, and Zhao (1994).
}
\usage{
	  DR(formula,model,data,treat,pscore)
}
\arguments{
  \item{formula}{An object of class \code{formula} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{model}{Set to \dQuote{lm} to use a linear model, set to \dQuote{glm} for the generalized linear model.}
  \item{data}{An optional data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model. If not found in data, the variables are taken from \code{environment(formula)}, typically the environment from which \code{CBPS} is called.}
  \item{treat}{A vector of treatment assignments, where 1 indicates assignment to the treatment group and 0 indicates assignment to the control group.}
  \item{pscore}{A vector of propensity scores.}
}
\details{Fits covariate balancing propensity scores.}
\value{
  Returns a vector with the point estimate of the doubly-robust estimator as well as its standard error.  For use with CBPS, run only with the ATE and binary treatment.
}
\references{Lunceford, Jared and Marie Davidian.  ``Stratification and Weighting Via the Propensity Score in Estimation of Causal Treatment Effects: A Comparative Study.'' Statistics in Medicine 23, 15 October 2004.}

\author{Marc Ratkovic, Kosuke Imai, and Christian Fong.}

\seealso{\link{IPW}}

\examples{

\dontrun{
###
### Example: Doubly robust estimator
###

##Load the LaLonde data
data(LaLonde)
## Estimate CBPS via logistic regression for ATE.  Run only with ATE.
fit <- CBPS(treat ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0), data = LaLonde, ATT = FALSE)
## Find doubly robust estimator under GLM
doubly.robust <- DR(re78 ~ age + educ + re75 + re74 + I(re75==0) + I(re74==0), model="glm", data = LaLonde, treat = treat, pscore = fitted(fit))
}
}