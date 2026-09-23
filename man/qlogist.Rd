\name{qlogist}
\alias{qlogist}
\title{Augmented (modified) logistic regression analyses for estimating risk ratio}
\usage{
qlogist(formula, data, eform=TRUE, cl=0.95, digits=4, var.method="MBN", id=NULL)
}
\arguments{
\item{formula}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}

\item{data}{A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

\item{eform}{A logical value that specify whether the outcome should be transformed by exponential function (default: \code{TRUE})}

\item{cl}{Confidence level for calculating confidence intervals (default: 0.95)}

\item{digits}{Number of decimal places in the output (default: 4).}

\item{var.method}{Method for estimating standard errors. Standard robust variance estimator (\code{standard}), Morel-Bokossa-Neerchaal-type corrected estimator (\code{MBN}), Gosho-Sato-Takeuchi-type corrected estimator (\code{GST}), and Wang-Long-type corrected estimator (\code{WL}) are available (default: \code{MBN}).}
\item{id}{Optional column name in \code{data}, supplied without quotes or as one character string, identifying independent clusters. With \code{NULL}, each original observation and its pseudo-observation form one cluster. With an external ID, all original and pseudo-observations with that ID are grouped together. Only \code{standard} and \code{MBN} support an external ID. Missing IDs in the analysis sample are not allowed.}
}
\value{
An object of class \code{"rqlm"} containing the following components.
\itemize{
\item \code{call}: The matched function call.
\item \code{formula}: The original model formula.
\item \code{coefficients}: Coefficient estimates on the log-risk scale, irrespective of \code{eform}.
\item \code{se}: Robust standard errors on the coefficient scale.
\item \code{cl}, \code{cu}: Lower and upper confidence limits on the coefficient scale.
\item \code{z}, \code{p}: Wald statistics and two-sided P-values.
\item \code{eform}, \code{cl.level}, \code{digits}, \code{var.method}: The corresponding analysis and display settings.
\item \code{vcov}: The full estimated covariance matrix, with coefficient names.
\item \code{model}: The fitted augmented-data \code{glm} object.
\item \code{n}: The number of original observations in the analysis sample.
\item \code{n.clusters}: The number of independent clusters.
}
}
\description{
Logistic regression with augmented pseudo-observations for estimating risk ratios is performed. This function is handled by a similar way with \code{lm} or \code{glm}. Also, the resultant coefficients and confidence limits can be transformed to exponential scales by specifying \code{eform}. The Morel-Bokossa-Neerchaal-type small-sample corrected estimator is adopted for standard error estimation as the default method.}
\details{
The model frame and design matrix are constructed from the original
observations before augmentation. This keeps transformations and matrix
predictors aligned with their original observations and avoids overwriting
user columns named \code{id} or \code{d}. Missing model values are handled
by \code{model.frame} using the current \code{na.action} option; IDs are
then aligned with the retained rows. Adaptive transformations such as
\code{poly} or splines are therefore defined before, rather than after,
augmentation.

At least two independent clusters are needed for \code{standard}.
For \code{MBN}, their number must exceed the number of model parameters.
The existing MBN formula, including the cluster-count adjustment in
\code{sandwich::vcovCL} before the additional MBN factors, is retained.
The \code{GST} and \code{WL} calculations are restricted to
\code{id = NULL}.

The methods \code{coef}, \code{vcov}, and \code{family} are available;
see \code{\link{rqlm-methods}}. The family returned by \code{family}
is the binomial-logit family of the augmented-data working model, not a
logit model for the original outcome. The coefficients represent log
risks in the original cohort. When using \code{dlnm::crosspred}, specify
\code{model.link = "log"} explicitly and use the outer object to retain
the robust covariance. The response predictions from \code{model}
refer to the augmented outcome; they are not the original risks.
}
\references{
Diaz-Quijano, F. A. (2012). A simple method for estimating relative risk using logistic regression. \emph{BMC Medical Research Methodology} \strong{12}, 14.

Gosho, M., Sato, Y., and Takeuchi, H. (2014). Robust covariance estimator for small-sample adjustment in the generalized estimating equations: a simulation study. \emph{Science Journal of Applied Mathematics and Statistics} \strong{2}, 20-25.

Morel, J. G., Bokossa, M., and Neerchal, N. (2003). Small sample correction for the variance of GEE estimators. \emph{Biometrical Journal} \strong{45}, 395-409.

Noma, H. (2026). Robust variance estimators for risk ratio estimators from logistic regression in cohort and case-cohort studies. \emph{Statistics and Probability Letters} \strong{234}, 110698.

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
