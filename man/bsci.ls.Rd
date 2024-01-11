\name{bsci.ls}
\alias{bsci.ls}
\title{Calculating bootstrap confidence interval for modified least-squares regression based on the quasi-score statistic}
\usage{
bsci.ls(formula, data, x.name=NULL, B=1000, cl=0.95, C0=10^-5,
 digits=4, seed=527916)
}
\arguments{
\item{formula}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}

\item{data}{A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

\item{x.name}{The variable name that the confidence interval is calculated for the regression coefficient; should be involved in \code{formula} as an explanatory variable. Specify as a character object.}

\item{B}{The number of bootstrap resampling (default: 1000)}

\item{cl}{Confidence level for calculating confidence intervals (default: 0.95)}

\item{C0}{A tuning parameter to control the precisions of numerical computations of confidence limits (default: 10^-5).}

\item{digits}{Number of decimal places in the output (default: 4).}

\item{seed}{Seed to generate random numbers (default: 527916).}
}
\value{
Results of the modified least-squares analyses are presented. Three objects are provided: Results of the modified least-squares regression with the Wald-type approximation by \code{rqlm}, the bootstrap-based confidence interval for the corresponding covariate, and P-value for the bootstrap test of \code{RD=0}.
}
\description{
Recent studies revealed the robust standard error estimates of the modified least-squares regression analysis are generally biased under small or moderate sample settings. To adjust the bias and to provide more accurate confidence intervals, confidence interval and P-value of the test for risk difference by modified least-squares regression are calculated based on the bootstrap approach of Noma et al. (2024).}
\references{
Noma, H. et al. (2024+). Bootstrap confidence intervals based on quasi-score statistics for the modified Poisson and least-squares regressions to estimate risk ratio and risk difference. Forthcoming.
}
\examples{
data(exdata01)

bsci.ls(y ~ x1 + x2 + x3 + x4, data=exdata01, "x3", B=10)
# For illustration. B should be >= 1000 (the number of bootstrap resampling).
}
