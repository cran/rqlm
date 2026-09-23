\name{rqlm-methods}
\alias{rqlm-methods}
\alias{coef.rqlm}
\alias{vcov.rqlm}
\alias{family.rqlm}
\alias{coef.ttemsm}
\alias{vcov.ttemsm}
\alias{family.ttemsm}
\title{Extract Coefficients, Covariance Matrices, and Model Families}
\usage{
\method{coef}{rqlm}(object, ...)
\method{vcov}{rqlm}(object, ...)
\method{family}{rqlm}(object, ...)
\method{coef}{ttemsm}(object, ...)
\method{vcov}{ttemsm}(object, ...)
\method{family}{ttemsm}(object, ...)
}
\arguments{
\item{object}{An object returned by \code{rqlm}, \code{qlogist}, or \code{ttemsm}.}
\item{...}{Additional arguments; currently unused.}
}
\description{
Extract coefficient estimates, the full estimated robust covariance
matrix, or the family of the underlying fitted working model.
}
\details{
The \code{qlogist} function returns an object of class \code{"rqlm"}
and uses the same methods. Coefficients and their covariance matrix
remain on the model scale, even when \code{eform = TRUE}. The covariance
matrix corresponds to the selected \code{var.method}.

For \code{qlogist}, \code{family} returns the binomial-logit family of
the augmented-data working model. Its coefficients parameterize log
risks for the original outcome, not ordinary log odds for that outcome.
For \code{dlnm::crosspred}, use \code{model.link = "log"} explicitly
for modified Poisson or augmented logistic risk-ratio models, and pass
the outer object to use its robust covariance. Supplying the inner
\code{glm} object instead uses that object's model-based covariance.
A \code{family} method alone does not ensure automatic link detection
by every other package.

Objects saved with an earlier version may not contain \code{vcov} or
\code{model}. They must be refitted before the corresponding extraction
method can be used. A full covariance matrix cannot be reconstructed
from marginal standard errors alone.
}
\value{
For \code{coef}, a named numeric coefficient vector.
For \code{vcov}, the estimated covariance matrix.
For \code{family}, the family object of the underlying fitted model.
}
\seealso{\code{\link{rqlm}}, \code{\link{qlogist}}, \code{\link{ttemsm}}}
\examples{
data(SGLT2i01)
fit <- rqlm(death36 ~ A + risk_score, data=SGLT2i01,
            family=poisson, eform=TRUE)
coef(fit)
vcov(fit)
family(fit)
}
\keyword{internal}