\name{rqlm}
\alias{rqlm}
\title{Modified Poisson and Least-Squares Regression Analyses for Binary Outcomes}
\usage{
rqlm(formula, data, family = poisson, eform = FALSE, cl = 0.95,
     digits = 4, var.method = "MBN")
}
\arguments{
\item{formula}{
An object of class \code{"formula"} (or one that can be coerced to that
class): a symbolic description of the model to be fitted.
}

\item{data}{
A data frame, list, or environment (or an object coercible by
\code{as.data.frame}) containing the variables in the model.
}

\item{family}{
A description of the error distribution and link function to be used in
the model. Specify \code{gaussian} for modified least-squares regression
or \code{poisson} for modified Poisson regression.
}

\item{eform}{
A logical value indicating whether the coefficient estimates and
confidence limits should be exponentiated (default: \code{FALSE}).
}

\item{cl}{
Confidence level used to calculate confidence intervals
(default: \code{0.95}).
}

\item{digits}{
Number of decimal places displayed in the printed output
(default: \code{4}).
}

\item{var.method}{
Method used to estimate the covariance matrix and standard errors.
Available options are:
\itemize{
  \item \code{"standard"}: the standard heteroskedasticity-robust
  sandwich covariance estimator;
  \item \code{"MBN"}: the Morel--Bokossa--Neerchal-type corrected
  estimator;
  \item \code{"GST"}: the Gosho--Sato--Takeuchi-type corrected
  estimator;
  \item \code{"WL"}: the Wang--Long-type corrected estimator;
  \item \code{"HAD"}: the Hadamard finite-sample unbiased covariance
  estimator;
  \item \code{"HAD0"}: the Hadamard estimator with negative
  observation-level variance estimates truncated at zero;
  \item \code{"HAD025"}: the Hadamard estimator with observation-level
  variance estimates truncated to the interval \eqn{[0, 1/4]}.
}
The \code{"HAD"}, \code{"HAD0"}, and \code{"HAD025"} options are
available only with a Gaussian family and identity link. The default is
\code{"MBN"}.
}
}
\value{
An object of class \code{"rqlm"}, containing the following components:
\item{call}{The matched function call.}
\item{formula}{The model formula.}
\item{coefficients}{Coefficient estimates on the model scale.}
\item{se}{Estimated robust standard errors.}
\item{cl}{Lower confidence limits on the model scale.}
\item{cu}{Upper confidence limits on the model scale.}
\item{z}{Wald \eqn{z} statistics.}
\item{p}{Two-sided Wald P-values.}
\item{eform}{The value supplied to \code{eform}.}
\item{cl.level}{The confidence level.}
\item{digits}{The number of displayed decimal places.}
\item{var.method}{The selected covariance estimator.}
\item{vcov}{The estimated covariance matrix.}
\item{model}{The fitted \code{glm} object.}
\item{n}{The number of observations in the analysis sample after
missing-value handling.}
\item{vhat.had}{
For a Hadamard method, the observation-level variance estimates after
any requested truncation; otherwise \code{NULL}.
}
\item{vhat.had.raw}{
For a Hadamard method, the untruncated observation-level variance
estimates; otherwise \code{NULL}.
}
\item{kappa.had}{
For a Hadamard method, the estimated condition number of the
Hadamard-squared residual-maker matrix; otherwise \code{NULL}.
}
\item{n.neg.had}{
For a Hadamard method, the number of negative untruncated
observation-level variance estimates.
}
\item{n.above025.had}{
For a Hadamard method, the number of untruncated observation-level
variance estimates exceeding \eqn{1/4}.
}
\item{min.vhat.had, max.vhat.had}{
The minimum and maximum of the untruncated observation-level variance
estimates.
}
}
\description{
Performs modified Poisson or modified least-squares regression for a
binary outcome. The function is used in a manner similar to \code{lm}
or \code{glm}. The regression family determines the fitted mean model,
and setting \code{eform = TRUE} exponentiates the coefficient estimates
and confidence limits. The Morel--Bokossa--Neerchal-type small-sample
corrected estimator is used by default.
}
\details{
The Hadamard methods are implemented for modified least-squares
regression fitted using a Gaussian family with an identity link.

Let \eqn{X} denote the model matrix,
\eqn{H = X(X^\top X)^{-1}X^\top} the hat matrix, and
\eqn{M = I - H} the residual-maker matrix. If \eqn{r} is the vector of
least-squares residuals, the \code{"HAD"} option estimates the
observation-specific error variances by
\deqn{
\widehat{\boldsymbol{\sigma}}^2 =
(M \mathbin{\circ} M)^{-1}(r \mathbin{\circ} r),
}
where \eqn{\mathbin{\circ}} denotes the elementwise Hadamard product.
The covariance matrix is then estimated as
\deqn{
(X^\top X)^{-1}
X^\top \mathrm{diag}(\widehat{\boldsymbol{\sigma}}^2)X
(X^\top X)^{-1}.
}
Under independent observations, a correctly specified linear
conditional-mean model, and invertibility of
\eqn{M \mathbin{\circ} M}, this covariance estimator is conditionally
finite-sample unbiased.

Unbiased observation-level variance estimates can be negative in a
particular sample. The \code{"HAD0"} option replaces negative values by
zero. For a Bernoulli outcome, whose conditional variance cannot exceed
\eqn{1/4}, \code{"HAD025"} truncates the estimates to
\eqn{[0, 1/4]}. These truncations can improve numerical stability but
do not preserve exact finite-sample unbiasedness.

The matrix \eqn{M \mathbin{\circ} M} is an \eqn{n \times n} matrix.
Consequently, the Hadamard methods can require substantial memory and
computation for large datasets. An error is returned if this matrix is
singular or numerically non-invertible, and a warning is returned when
it is severely ill-conditioned.
}
\references{
Cheung, Y. B. (2007).
A modified least-squares regression approach to the estimation of risk
difference.
\emph{American Journal of Epidemiology} \strong{166}, 1337--1344.

Gosho, M., Ishii, R., Noma, H., and Maruo, K. (2023).
A comparison of bias-adjusted generalized estimating equations for
sparse binary data in small-sample longitudinal studies.
\emph{Statistics in Medicine} \strong{42}, 2711--2727.

Gosho, M., Sato, Y., and Takeuchi, H. (2014).
Robust covariance estimator for small-sample adjustment in the
generalized estimating equations: a simulation study.
\emph{Science Journal of Applied Mathematics and Statistics}
\strong{2}, 20--25.

Morel, J. G., Bokossa, M., and Neerchal, N. (2003).
Small sample correction for the variance of GEE estimators.
\emph{Biometrical Journal} \strong{45}, 395--409.

Noma, H. and Gosho, M. (2025).
Finite-sample improved confidence intervals based on the estimating
equation theory for the modified Poisson and least-squares regressions.
\emph{Epidemiologic Methods} \strong{14}, 20240030.

Noma, H. and Kitano, T. (2026).
Finite-sample unbiased covariance estimation for modified least-squares regression with binary outcomes.
In Preparation.

Noma, H., Sunada, H., and Gosho, M. (2025).
Quasi-likelihood ratio tests and the Bartlett-type correction for
improved inferences of the modified Poisson and least-squares
regressions for binary outcomes.
\emph{Statistica Neerlandica} \strong{79}, e70012.

Wang, M. and Long, Q. (2011).
Modified robust variance estimator for generalized estimating
equations with improved small-sample performance.
\emph{Statistics in Medicine} \strong{30}, 1278--1291.

White, H. (1982).
Maximum likelihood estimation of misspecified models.
\emph{Econometrica} \strong{50}, 1--25.

Zou, G. (2004).
A modified Poisson regression approach to prospective studies with
binary data.
\emph{American Journal of Epidemiology} \strong{159}, 702--706.
}
\examples{
data(exdata02)

rqlm(
  y ~ x1 + x2 + x3 + x4,
  data = exdata02,
  family = poisson,
  eform = TRUE
)
# Modified Poisson regression analysis.
# Coefficients and confidence limits are displayed as risk ratios.
# The MBN robust covariance estimator is used by default.

rqlm(
  y ~ x1 + x2 + x3 + x4,
  data = exdata02,
  family = gaussian
)
# Modified least-squares regression analysis.
# Coefficients are risk differences.

rqlm(
  y ~ x1 + x2 + x3 + x4,
  data = exdata02,
  family = gaussian,
  var.method = "HAD"
)
# Modified least-squares regression with the finite-sample unbiased
# Hadamard covariance estimator.

rqlm(
  y ~ x1 + x2 + x3 + x4,
  data = exdata02,
  family = gaussian,
  var.method = "HAD0"
)
# Negative observation-level variance estimates are truncated at zero.

rqlm(
  y ~ x1 + x2 + x3 + x4,
  data = exdata02,
  family = gaussian,
  var.method = "HAD025"
)
# Observation-level variance estimates are truncated to [0, 1/4].

rqlm(
  y ~ x1 + x2 + x3 + x4,
  data = exdata02,
  family = poisson,
  eform = TRUE,
  var.method = "GST"
)

rqlm(
  y ~ x1 + x2 + x3 + x4,
  data = exdata02,
  family = poisson,
  eform = TRUE,
  var.method = "WL"
)
}
