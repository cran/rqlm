\name{SumStat}
\alias{SumStat}
\title{Creating summary table for IPTW analysis using stabilized weights}
\usage{
SumStat(formula, data, trunc=c(0.01,0.99), digits=3)
}
\arguments{
\item{formula}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the logistic regression model to be fitted.}

\item{data}{A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

\item{trunc}{Quantiles to be truncated the weights (default: 0.01, 0.99).}

\item{digits}{Number of decimal places in the output (default: 4).}
}
\value{
Summary statistics for unadjusted and IPTW-adjusted populations are provided. For each covariate, the mean and standard deviation, as well as the standardized mean difference (SMD), were calculated. Truncated stabilized weights were used for weighting. Setting the \code{trunc} argument to c(0, 1) allows weighting with untruncated stabilized weights.
\itemize{
\item \code{mean0}: Mean for group 0.
\item \code{mean1}: Mean for group 1.
\item \code{sd0}: SD for group 0.
\item \code{sd1}: SD for group 1.
\item \code{SMD}: SMD for the original population.
\item \code{wmean0}: IPTW-weighted mean for group 0.
\item \code{wmean1}: IPTW-weighted mean for group 1.
\item \code{wsd0}: IPTW-weighted SD for group 0.
\item \code{wsd1}: IPTW-weighted SD for group 1.
\item \code{wSMD}: IPTW-weighted SMD for the adjusted pseudo-population.
\item \code{type}: Type of the variable.
}
}
\description{
Summary statistics are computed before and after weighting for IPTW analyses. For all covariates included in the weight model, the mean, standard deviation, and standardized mean difference (SMD) are calculated. Weighting can be performed using either stabilized or unstabilized weights, with optional truncation. For binary covariates, summary statistics are computed using unbiased estimators of the binomial mean and variance.
}
\references{
Austin, P. C. (2009). Balance diagnostics for comparing the distribution of baseline covariates between treatment groups in propensity-score matched samples. \emph{Statistics in Medicine} \strong{28}, 3083-3107.
}
\examples{
data(exdata02)

SumStat(x2 ~ x1 + x3 + x4, data=exdata02)
}
