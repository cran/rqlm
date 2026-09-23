coef.rqlm <- function(object, ...) {
  object$coefficients
}

vcov.rqlm <- function(object, ...) {
  if (is.null(object$vcov))
    stop("This object does not contain a covariance matrix; refit it with the updated package.")
  object$vcov
}

family.rqlm <- function(object, ...) {
  if (is.null(object$model))
    stop("This object does not contain a fitted model; refit it with the updated package.")
  stats::family(object$model)
}

coef.ttemsm <- coef.rqlm
vcov.ttemsm <- vcov.rqlm
family.ttemsm <- family.rqlm

.rqlm_id <- function(data, expr, omit=NULL, n=NULL) {
  if (is.null(expr)) return(NULL)
  if (is.symbol(expr)) expr <- as.character(expr)
  if (!is.character(expr) || length(expr) != 1L ||
      is.na(expr) || !(expr %in% names(data)))
    stop("id must be NULL or a column name in data.")
  z <- data[[expr]]
  if (!is.atomic(z) || !is.null(dim(z)))
    stop("The id column must be a vector.")
  if (!is.null(omit)) z <- z[-as.integer(omit)]
  if (!is.null(n) && length(z) != n)
    stop("The id column could not be aligned with the analysis sample.")
  if (!is.null(n) && anyNA(z))
    stop("The id column contains missing values in the analysis sample.")
  if (is.factor(z)) z <- droplevels(z)
  z
}
