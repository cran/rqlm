\name{ttemsm}
\alias{ttemsm}
\title{Pooled logistic regression for target trial emulation}
\usage{
ttemsm(formula, data, id, weight, family=quasibinomial(link="cloglog"), 
eform=TRUE, cl=0.95, digits=4, var.method="MBN")
}
\arguments{
\item{formula}{An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.}

\item{data}{A data frame, list or environment (or object coercible by as.data.frame to a data frame) containing the variables in the model.}

\item{id}{id variable for individual participants containing in the \code{data} object.}

\item{weight}{weight variable for individual participants containing in the \code{data} object.}

\item{family}{A description of the error distribution and link function to be used in the model. When the goal is to estimate hazard ratios using a discrete-time Cox model, the complementary log–log (cloglog) link should be used. When the objective is to estimate risk differences or cumulative incidence (e.g., via the g-formula or standardization), the logit link should be used.}

\item{eform}{A logical value that specify whether the outcome should be transformed by exponential function (default: \code{TRUE})}

\item{cl}{Confidence level for calculating confidence intervals (default: 0.95)}

\item{digits}{Number of decimal places in the output (default: 4).}

\item{var.method}{Method for estimating standard errors. Standard cluster-robust variance estimator (\code{standard}) and Morel-Bokossa-Neerchaal-type corrected cluster-robust estimator (\code{MBN}) are available (default: \code{MBN}).}
}
\value{
Results of the pooled logistic regression analysis.
\itemize{
\item \code{coef}: Coefficient estimates; transformed to the exponential scale if \code{eform=TRUE}.
\item \code{SE}: Cluster-robust standard error estimates for \code{coef}.
\item \code{CL}: Lower limits of confidence intervals.
\item \code{CU}: Upper limits of confidence intervals.
\item \code{P-value}: P-values for the coefficient tests.
}
}
\description{
This function implements pooled logistic regression for use in target trial emulation. Before running the function, the user must prepare an analysis dataset by stacking the sequential‐trial datasets in long format. If inverse probability weights (such as IPCW) are required, the corresponding \code{weight} variable should be included in the dataset in advance. The regression model is specified through the \code{formula} argument, and the individual identifier must be provided via the \code{id} argument, so that cluster‐robust standard errors are computed across repeated trials for each individual. When a \code{weight} variable is specified through weights, the function performs inverse probability weighting. The output includes point estimates of the hazard ratio, corresponding confidence intervals, and P-values.
}
\references{
Gosho, M., Ishii, R., Noma, H., and Maruo, K. (2023). A comparison of bias-adjusted generalized estimating equations for sparse binary data in small-sample longitudinal studies. \emph{Statistics in Medicine} \strong{42}, 2711-2727.

Hernan, M. A., Alonso, A., Logan, R., et al. (2008). Observational studies analyzed like randomized experiments: an application to postmenopausal hormone therapy and coronary heart disease. \emph{Epidemiology} \strong{19}, 766-779.

Hernan, M. A., Dahabreh, I. J., Dickerman, B. A., and Swanson, S. A. (2025). The Target Trial Framework for Causal Inference From Observational Data: Why and When Is It Helpful? \emph{Annals of Internal Medicine} \strong{178}, 402-407.

Hernan, M. A., Wang, W., and Leaf, D. E. (2022). Target Trial Emulation: A Framework for Causal Inference From Observational Data. \emph{JAMA} \strong{328}, 2446-2447.

Morel, J. G., Bokossa, M., and Neerchal, N. (2003). Small sample correction for the variance of GEE estimators. \emph{Biometrical Journal} \strong{45}, 395-409.
}
\examples{
data(exdata04)

ttemsm( Y ~ A + L1 + L2 + L3 + time + I(time^2) + trial,
  data    = exdata04, id = ID, weight = w_pp,
  family = quasibinomial(link="cloglog"),  
  eform   = TRUE, cl = 0.95, var.method="standard")
# Pooled logistic regression for target trial emulation with cloglog link
# For estimating hazard ratios using a discrete-time Cox model

ttemsm( Y ~ A + L1 + L2 + L3 + time + I(time^2) + trial,
  data    = exdata04, id = ID, weight = w_pp,
  family = quasibinomial(link="logit"),
  eform   = TRUE, cl = 0.95, var.method="standard")
# Pooled logistic regression for target trial emulation with logit link
# For estimating RDs or cumulative incidence (e.g., via the g-formula)

ttemsm( Y ~ A + L1 + L2 + L3 + time + I(time^2) + trial,
  data    = exdata04, id = ID, weight = w_pp,
  eform   = TRUE, cl = 0.95, var.method="MBN")
# Pooled logistic regression for target trial emulation with cloglog link
# Morel-Bokossa-Neerchaal-type corrected SE estimator is used.
}
