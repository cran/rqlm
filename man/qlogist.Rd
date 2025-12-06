\name{qlogist}
\alias{qlogist}
\title{Augmented (modified) logistic regression analyses for estimating risk ratio}
\usage{
qlogist(formula, data, eform=TRUE, cl=0.95, digits=4, var.method="MBN")
}
\arguments{
\item{formula}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}

\item{data}{A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

\item{eform}{A logical value that specify whether the outcome should be transformed by exponential function (default: \code{TRUE})}

\item{cl}{Confidence level for calculating confidence intervals (default: 0.95)}

\item{digits}{Number of decimal places in the output (default: 4).}

\item{var.method}{Method for estimating standard errors. Standard robust variance estimator (\code{standard}), Morel-Bokossa-Neerchaal-type corrected estimator (\code{MBN}), Gosho-Sato-Takeuchi-type corrected estimator (\code{GST}), and Wang-Long-type corrected estimator (\code{WL}) are available (default: \code{MBN}).}
}
\value{
Results of the augmented (modified) logistic regression analysis.
\itemize{
\item \code{coef}: Coefficient estimates; transformed to the exponential scale if \code{eform=TRUE}.
\item \code{SE}: Robust standard error estimates for \code{coef}.
\item \code{CL}: Lower limits of confidence intervals.
\item \code{CU}: Upper limits of confidence intervals.
\item \code{P-value}: P-values for the coefficient tests.
}
}
\description{
Logistic regression with augmented pseudo-observations for estimating risk ratios is performed. This function is handled by a similar way with \code{lm} or \code{glm}. Also, the resultant coefficients and confidence limits can be transformed to exponential scales by specifying \code{eform}. The Morel-Bokossa-Neerchaal-type small-sample corrected estimator is adopted for standard error estimation as the default method.}
\references{
Diaz-Quijano, F. A. (2012). A simple method for estimating relative risk using logistic regression. \emph{BMC Medical Research Methodology} \strong{12}, 14.

Gosho, M., Sato, Y., and Takeuchi, H. (2014). Robust covariance estimator for small-sample adjustment in the generalized estimating equations: a simulation study. \emph{Science Journal of Applied Mathematics and Statistics} \strong{2}, 20-25.

Morel, J. G., Bokossa, M., and Neerchal, N. (2003). Small sample correction for the variance of GEE estimators. \emph{Biometrical Journal} \strong{45}, 395-409.

Noma, H. (2025). Robust variance estimators for risk ratio estimators from logistic regression in cohort and case-cohort studies. Forthcoming.

Noma, H., and Gosho, M. (2025). Logistic mixed-effects model analysis with pseudo-observations for estimating risk ratios in clustered binary data analysis. \emph{Statistics in Medicine} \strong{44}, e70280.

Schouten, E. G., Dekker, J. M., Kok, F. J., et al. (1993). Risk ratio and rate ratio estimation in case-cohort designs: hypertension and cardiovascular mortality. \emph{Statistics in Medicine} \strong{12}, 1733-1745.

Shiiba, H., and Noma, H. (2025). Confidence intervals of risk ratios for the augmented logistic regression with pseudo-observations. \emph{Stats} \strong{8}, 83.

Wang, M., and Long, Q. (2011). Modified robust variance estimator for generalized estimating equations with improved small-sample performance. \emph{Statistics in Medicine} \strong{30}, 1278-1291.
}
\examples{
data(exdata02)

qlogist(y ~ x1 + x2 + x3 + x4, data=exdata02)
# Augmented logistic regression analysis
# Coefficient estimates are translated to risk ratio scales
# MBN robust variance estimator is adopted.

qlogist(y ~ x1 + x2 + x3 + x4, data=exdata02, var.method="GST")
# GST robust variance estimator is adopted.

qlogist(y ~ x1 + x2 + x3 + x4, data=exdata02, var.method="WL")
# WL robust variance estimator is adopted.
}
