\name{coeff}
\alias{coeff}
\title{Computation of the ordinary confidence intervals and P-values using the model variance estimator}
\usage{
coeff(gm, eform=FALSE, cl=0.95, digits=4)
}
\arguments{
\item{gm}{An output object of \code{lm}, \code{glm}, \code{lmer}, or \code{glmer}.}

\item{eform}{A logical value that specify whether the outcome should be transformed by exponential function (default: \code{FALSE})}

\item{cl}{Confidence level for calculating confidence intervals (default: 0.95)}

\item{digits}{Number of decimal places in the output (default: 4).}
}
\value{
Results of inferences of the regression coefficients using the ordinary model variance estimators.
\itemize{
\item \code{coef}: Coefficient estimates; transformed to the exponential scale if \code{eform=TRUE}.
\item \code{SE}: Robust standard error estimates for \code{coef}.
\item \code{CL}: Lower limits of confidence intervals.
\item \code{CU}: Upper limits of confidence intervals.
\item \code{P-value}: P-values for the coefficient tests.
}
}
\description{
Confidence intervals and P-values for the generalized linear model and generalized linear-mixed-effects model can be calculated using the ordinary model variance estimators. Through simply entering the output objects of \code{lm}, \code{glm}, \code{lmer}, or \code{glmer}, the inference results are fastly computed. For the linear regression model, the exact confidence intervals and P-values based on the t-distribution are calculated. Also, for the generalized linear model, the Wald-type confidence intervals and P-values based on the asymptotic normal approximation are computed. The resultant coefficients and confidence limits can be transformed to exponential scales by specifying \code{eform}.
}
\examples{
data(exdata02)

gm1 <- glm(y ~ x1 + x2 + x3 + x4, data=exdata02, family=binomial)
coeff(gm1,eform=TRUE)
# Logistic regression analysis
# Coefficient estimates are translated to odds ratio scales

lm1 <- lm(x1 ~ x2 + x3 + x4, data=exdata02)
coeff(lm1)
# Linear regression analysis

data(mch)

if(requireNamespace("lme4", quietly = TRUE)) {
    lmr1 <- lme4::lmer(Y ~ x + (1|SOUM), data=mch)
    coeff(lmr1)
}
# Linear mixed-effects model analysis

if(requireNamespace("lme4", quietly = TRUE)) {
	gmr1 <- lme4::glmer(y ~ x + (1|SOUM), nAGQ=25, family=binomial, data=mch)
	coeff(gmr1, eform=TRUE)
}
# Logistic mixed-effects model analysis
# Coefficient estimates are translated to odds ratio scales
}
