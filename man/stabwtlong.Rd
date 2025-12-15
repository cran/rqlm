\name{stabwtlong}
\alias{stabwtlong}
\title{Calculating stabilized weights for IPW analysis: Longitudinal data}
\usage{
stabwtlong(formula_denom, formula_num, data, trunc=c(0.01,0.99), digits=4)
}
\arguments{
\item{formula_denom}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the denominator logistic regression model to be fitted.}

\item{formula_num}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the numerator logistic regression model to be fitted.}

\item{data}{A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

\item{trunc}{Quantiles to be truncated the weights (default: 0.01, 0.99).}

\item{digits}{Number of decimal places in the output (default: 4).}
}
\value{
Truncated and untruncated stabilized weights are calculated.
\itemize{
\item \code{sw1}: Untruncated stabilized weights.
\item \code{sw2}: Truncated stabilized weights.
}
}
\description{
This function calculates stabilized weights for IPW analysis of longitudinal data using logistic regressions. Both of untruncated and utruncated weights are provided.
}
\references{
Cole, S. R., and Hernan, M. A. (2008). Constructing inverse probability weights for marginal structural models. \emph{American Journal of Epidemiology} \strong{168}, 656-664.

Westreich, D., Edwards, J. K., Lesko, C. R., Stuart, E., and Cole, S. R. (2017). Transportability of trial results using inverse odds of sampling weights. \emph{American Journal of Epidemiology} \strong{186}, 1010-1014.
}
\examples{
data(exdata04)

stabwtlong(formula_denom = A ~ L1 + L2 + L3, formula_num = A ~ L1, data = exdata04)
}
