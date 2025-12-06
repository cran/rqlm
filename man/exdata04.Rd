\docType{data}
\name{exdata04}
\alias{exdata04}
\title{A simulated example dataset for target trial emulation}
\format{A simulated data for sequentially nested emulated trials (n=600).}
\usage{
data(exdata04)
}
\description{
A simulated data for sequentially nested emulated trials.
\itemize{
  \item \code{ID}: The id variable for individual participants.
  \item \code{trial}: The trial id variable.
  \item \code{time}: The discrete time index within each emulated trial, representing the interval since trial baseline.
  \item \code{A}: Treatment variable.
  \item \code{L1}: Baseline confounding variable.
  \item \code{Y}: Dichotomous outcome variable.
  \item \code{w_pp}: The weight variable.
}
}
\keyword{datasets}
