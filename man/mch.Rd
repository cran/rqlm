\docType{data}
\name{mch}
\alias{mch}
\title{A cluster-randomised trial dataset for the maternal and child health handbook}
\format{A data frame with 500 participants with 18 soums.}
\usage{
data(mch)
}
\description{
A cluster-randomised trial dataset with binomial outcome.
\itemize{
  \item \code{ID}: ID variable of participants.
  \item \code{SOUM}: ID variable of soums (involving 18 soums).
  \item \code{x}: Binary variable specifying intervention groups (1=Intervention, 0=Control).
  \item \code{mage}: Mother's age.
  \item \code{medu}: Mother's education (1=uneducated, 2=elementary, 3=incomplete secondary, 4=complete secondary, 5=incomplete high, 6=high (completed collage or university)).
  \item \code{mmarry}: Mother's marital status (1=single, 2=married/cohabitating, 3=separated/divorce, 4=windowed/other).
  \item \code{mprig1}: First pregnancy (1=Yes, 2=No).
  \item \code{height}: Mother's height.
  \item \code{weight}: Mother's weight.
  \item \code{time}: Travel time from mother's home to antenatal care clinic.
  \item \code{Y}: Outcome variable: Number of antenatal visits.
  \item \code{y}: Outcome variable: Whether the number of antenatal visits is >= 6 (0 or 1).
  \item \code{ses}: Quintile groups by the social-economic index (= 1, 2, 3, 4, 5).
}
}
\references{
Mori, R., Yonemoto, N., Noma, H., et al. (2015). 
The Maternal and Child Health (MCH) handbook in Mongolia: a cluster-randomized, controlled trial. 
\emph{PloS One} 
\strong{10}: e0119772.
}
\keyword{datasets}
