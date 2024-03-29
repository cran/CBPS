% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CBPSMain.R
\name{plot.CBPS}
\alias{plot.CBPS}
\title{Plotting Covariate Balancing Propensity Score Estimation}
\usage{
\method{plot}{CBPS}(x, covars = NULL, silent = TRUE, boxplot = FALSE, ...)
}
\arguments{
\item{x}{an object of class \dQuote{CBPS} or \dQuote{npCBPS}, usually, a
result of a call to \code{CBPS} or \code{npCBPS}.}

\item{covars}{Indices of the covariates to be plotted (excluding the
intercept).  For example, if only the first two covariates from
\code{balance} are desired, set \code{covars} to 1:2.  The default is
\code{NULL}, which plots all covariates.}

\item{silent}{If set to \code{FALSE}, returns the imbalances used to
construct the plot.  Default is \code{TRUE}, which returns nothing.}

\item{boxplot}{If set to \code{TRUE}, returns a boxplot summarizing the
imbalance on the covariates instead of a point for each covariate.  Useful
if there are many covariates.}

\item{...}{Additional arguments to be passed to plot.}
}
\value{
For binary and multi-valued treatments, plots the absolute
difference in standardized means by contrast for all covariates before and
after weighting.  This quantity for a single covariate and a given pair of
treatment conditions is given by \eqn{\frac{\sum_{i=1}^{n} w_i * (T_i == 1)
* X_i}{\sum_{i=1}^{n} (T_i == 1) * w_i} - \frac{\sum_{i=1}^{n} w_i * (T_i ==
0) * X_i}{\sum_{i=1}^{n} (T_i == 0) * w_i}}{[\sum w_i * (T_i == 1) *
X_i]/[\sum w_i * (T_i == 1)] - [\sum w_i * (T_i == 0) * X_i]/[\sum w_i *
(T_i == 0)]}.  For continuous treatments, plots the weighted absolute
Pearson correlation between the treatment and each covariate.  See
\url{https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Weighted_correlation_coefficient.
}
}
\description{
This function plots the absolute difference in standardized means before and after
weighting.  To access more sophisticated graphics for assessing covariate balance,
consider using Noah Greifer's \code{cobalt} package.
}
\details{
The "Before Weighting" plot gives the balance before weighting, and the
"After Weighting" plot gives the balance after weighting.

### @aliases plot.CBPS plot.npCBPS
}
\seealso{
\link{CBPS}, \link{plot}
}
\author{
Christian Fong, Marc Ratkovic, and Kosuke Imai.
}
