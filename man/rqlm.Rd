\name{rqlm}
\alias{rqlm}
\title{Modified Poisson and least-squares regression analyses for binary outcomes}
\usage{
rqlm(formula, data, family=poisson, eform=FALSE, cl=0.95, digits=4)
}
\arguments{
\item{formula}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}

\item{data}{A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

\item{family}{A description of the error distribution and link function to be used in the model. \code{gaussian}: Modified least-squares regression. \code{poisson}: Modified Poisson regression.}

\item{eform}{A logical value that specify whether the outcome should be transformed by exponential function (default: \code{FALSE})}

\item{cl}{Confidence level for calculating confidence intervals (default: 0.95)}

\item{digits}{Number of decimal places in the output (default: 4).}
}
\value{
Results of the modified Poisson and least-squares regression analyses.
\itemize{
\item \code{coef}: Coefficient estimates; transformed to the exponential scale if \code{eform=TRUE}.
\item \code{SE}: Robust standard error estimates for \code{coef}.
\item \code{CL}: Lower limits of confidence intervals.
\item \code{CU}: Upper limits of confidence intervals.
\item \code{P-value}: P-values for the coefficient tests.
}
}
\description{
Modified Poisson and least-squares regression analyses for binary outcomes are performed. This function is handled by a similar way with \code{lm} or \code{glm}. The model fitting to the binary data can be specified by \code{family}. Also, the resultant coefficients and confidence limits can be transformed to exponential scales by specifying \code{eform}. The standard error estimates are calculated using the standard robust variance estimator by \code{sandwich} package.}
\references{
Cheung, Y. B. (2007). A modified least-squares regression approach to the estimation of risk difference. \emph{American Journal of Epidemiology} \strong{166}, 1337-1344.

Noma, H. and Gosho, M. (2024). Bootstrap confidence intervals based on quasi-likelihood estimating functions for the modified Poisson and least-squares regressions for binary outcomes. Forthcoming.

White, H. (1982). Maximum likelihood estimation of misspecified models. \emph{Econometrica}, \strong{50}, 1-25.

Zou, G. (2004). A modified poisson regression approach to prospective studies with binary data. \emph{American Journal of Epidemiology} \strong{159}, 702-706.
}
\examples{
data(exdata02)

rqlm(y ~ x1 + x2 + x3 + x4, data=exdata02, family=poisson, eform=TRUE)
# Modifed Poisson regression analysis
# Coefficient estimates are translated to risk ratio scales

rqlm(y ~ x1 + x2 + x3 + x4, data=exdata02, family=gaussian)
# Modifed least-squares regression analysis

rqlm(y ~ x1 + x2 + x3 + x4, data=exdata02, family=gaussian, digits=3)
# Modifed least-squares regression analysis
# Number of decimal places can be changed by specifying "digits"
}
