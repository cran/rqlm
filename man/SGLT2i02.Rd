\name{SGLT2i02}
\alias{SGLT2i02}
\docType{data}
\title{Simulated example data based on a target trial emulation for SGLT2i: risk difference estimation}
\usage{data(SGLT2i02)}
\description{
Simulated example data representing an active-comparator new-user
comparison of SGLT2i and DPP-4i initiation in older adults with type 2
diabetes. The outcome is binary at a fixed 36-month horizon. These are
not patient records and provide no clinical evidence about either drug.

\itemize{
  \item \code{id}: Unique synthetic individual identifier; one row per individual.
  \item \code{A}: Treatment indicator: 0 for DPP-4i and 1 for SGLT2i.
  \item \code{treatment}: Treatment factor with levels DPP-4i and SGLT2i.
  \item \code{death36}: Simulated all-cause death by 36 months, coded 0/1.
  \item \code{age}: Baseline age in years.
  \item \code{female}: Female indicator, coded 0/1.
  \item \code{bmi}: Body mass index in kg/m^2.
  \item \code{hba1c}: Glycated haemoglobin in percent.
  \item \code{egfr}: Estimated glomerular filtration rate in mL/min/1.73 m^2.
  \item \code{proteinuria}: Proteinuria indicator, coded 0/1.
  \item \code{prior_heart_failure}: History of heart failure, coded 0/1.
  \item \code{prior_stroke}: History of stroke, coded 0/1.
  \item \code{recent_hospitalization}: Recent hospitalization indicator, coded 0/1.
  \item \code{trial_period}: Factor identifying the synthetic treatment-initiation period.
  \item \code{risk_score}: Fixed baseline-covariate score used in the simulation; not a validated clinical score.
  \item \code{followup_months}: Common ascertainment horizon, equal to 36 for every individual; not an event time.
}
}

\format{
A data frame with 700 observations and 16 variables.
}
\details{
The datasets were independently simulated and include no participant records
from the motivating study. They are intended solely for education, software
testing, and reproducible examples; numerical results obtained from them have
no clinical interpretation.
}
\references{
Noma, H., Goto, A., Sugimoto, T., Sunada, H., Oda, F., Maeda, M., and
Fukuda, H. (2026). Real-world effectiveness of SGLT2 inhibitors in adults aged
75 years or older: a target trial emulation. \emph{Age and Ageing} \strong{55}, afag246.
}
\seealso{\code{\link{rqlm}}, \code{\link{qlogist}}, \code{\link{SGLT2i01}}}
\examples{
data(SGLT2i02)
rqlm(death36 ~ A + I(exp(2.4 * risk_score)),
      data=SGLT2i02, family=gaussian, var.method="standard")
\donttest{
rqlm(death36 ~ A + I(exp(2.4 * risk_score)),
     data=SGLT2i02, family=gaussian, var.method="HAD")
}
}
\keyword{datasets}
