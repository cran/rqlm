rqlm <- function(formula, data, family = poisson, eform = FALSE,
                 cl = 0.95, digits = 4, var.method = "MBN") {

  call <- match.call()

  var.method <- match.arg(
    var.method,
    choices = c("standard", "MBN", "GST", "WL",
                "HAD", "HAD0", "HAD025")
  )

  gm1 <- glm(
    formula,
    data = data,
    family = family,
    x = TRUE,
    y = TRUE
  )

  ## Use the actual analysis sample after handling missing values.
  n <- nrow(gm1$x)

  cc <- 1 - 0.5 * (1 - cl)

  coef1 <- gm1$coefficients
  Vout <- NULL

  ## Hadamard-estimator diagnostics.
  vhat.had <- NULL
  vhat.had.raw <- NULL
  kappa.had <- NULL
  n.neg.had <- NA_integer_
  n.above025.had <- NA_integer_
  min.vhat.had <- NA_real_
  max.vhat.had <- NA_real_

  if (var.method == "standard") {

    V1 <- sandwich::sandwich(gm1)
    Vout <- V1
    se1 <- sqrt(diag(V1))

  }

  if (var.method == "MBN") {

    V1 <- sandwich::sandwich(gm1)
    Ainv <- vcov(gm1)
    A <- solve(Ainv)

    p1 <- ncol(V1)

    if (n <= p1) {
      stop("The number of observations must be larger than the number of model parameters for var.method = \"MBN\".")
    }

    Q1 <- (n - 1) / (n - p1)
    Q2 <- n / (n - 1)

    Q3 <- sum(diag(V1 %*% A)) / p1

    delta <- min(0.5, p1 / (n - p1))
    gamma <- max(1, Q3)

    V2 <- Q1 * Q2 * V1 + delta * gamma * Ainv
    Vout <- V2
    se1 <- sqrt(diag(V2))

  }

  if (var.method == "GST") {

    Ainv <- vcov(gm1)

    p1 <- ncol(Ainv)

    X <- gm1$x
    pi <- gm1$fitted.values
    y <- gm1$y

    BG1 <- BG2 <- matrix(numeric(p1 * p1), p1)

    for (i in seq_len(n)) {

      xi <- matrix(X[i, ], ncol = 1)

      Ai <- pi[i] * (1 - pi[i]) * (xi %*% t(xi))
      Ui <- xi * (y[i] - pi[i])

      Aiih <- mat_inv_sqrt(Ai)
      BG1 <- BG1 + Aiih %*% Ui %*% t(Ui) %*% Aiih

    }

    BG1 <- BG1 / (n - p1)

    for (i in seq_len(n)) {

      xi <- matrix(X[i, ], ncol = 1)

      Ai <- pi[i] * (1 - pi[i]) * (xi %*% t(xi))
      Aih <- mat_sqrt(Ai)
      BG2 <- BG2 + Aih %*% BG1 %*% Aih

    }

    V2 <- Ainv %*% BG2 %*% Ainv
    V2 <- (V2 + t(V2)) / 2

    Vout <- V2
    se1 <- sqrt(diag(V2))

  }

  if (var.method == "WL") {

    Ainv <- vcov(gm1)

    p1 <- ncol(Ainv)

    X <- gm1$x
    pi <- gm1$fitted.values
    y <- gm1$y

    BW1 <- BW2 <- matrix(numeric(p1 * p1), p1)

    for (i in seq_len(n)) {

      xi <- matrix(X[i, ], ncol = 1)

      Ai <- pi[i] * (1 - pi[i]) * (xi %*% t(xi))
      Ui <- xi * (y[i] - pi[i])

      Aih <- mat_sqrt(Ai)
      Aiih <- mat_inv_sqrt(Ai)
      Hi <- Aih %*% Ainv %*% Aih
      Fi <- solve(diag(p1) - Hi)

      BW1 <- BW1 + Aiih %*% Fi %*% Ui %*% t(Ui) %*% Fi %*% Aiih

    }

    BW1 <- BW1 / n

    for (i in seq_len(n)) {

      xi <- matrix(X[i, ], ncol = 1)

      Ai <- pi[i] * (1 - pi[i]) * (xi %*% t(xi))
      Aih <- mat_sqrt(Ai)
      BW2 <- BW2 + Aih %*% BW1 %*% Aih

    }

    V2 <- Ainv %*% BW2 %*% Ainv
    V2 <- (V2 + t(V2)) / 2

    Vout <- V2
    se1 <- sqrt(diag(V2))

  }

  if (var.method %in% c("HAD", "HAD0", "HAD025")) {

    ## These estimators are defined here only for unweighted
    ## ordinary least-squares regression with an identity link.
    if (!identical(gm1$family$family, "gaussian") ||
        !identical(gm1$family$link, "identity")) {
      stop(
        paste0(
          "var.method = \"", var.method,
          "\" is available only for modified least-squares regression ",
          "with family = gaussian or gaussian(link = \"identity\")."
        )
      )
    }

    X <- gm1$x
    y <- as.vector(gm1$y)
    p1 <- ncol(X)

    if (qr(X)$rank < p1 || anyNA(coef1)) {
      stop("The design matrix is rank deficient; the Hadamard estimator cannot be computed.")
    }

    if (!all(y %in% c(0, 1))) {
      warning(
        paste0(
          "The response is not coded entirely as 0/1. ",
          "HAD and HAD0 remain heteroskedastic OLS covariance estimators, ",
          "but the upper truncation at 1/4 used by HAD025 is justified ",
          "specifically for Bernoulli outcomes."
        ),
        call. = FALSE
      )
    }

    ## OLS bread matrix.
    XtX.inv <- tryCatch(
      solve(crossprod(X)),
      error = function(e) {
        stop(
          "X'X could not be inverted; the Hadamard estimator cannot be computed.",
          call. = FALSE
        )
      }
    )

    ## Residual-maker matrix M = I - H.
    H <- X %*% XtX.inv %*% t(X)
    M <- diag(n) - H

    ## Elementwise (Hadamard) square of M.
    M2 <- M * M

    ## E(r^2 | X) = (M o M) sigma^2.
    r2 <- as.vector(residuals(gm1, type = "response"))^2

    vhat.had.raw <- tryCatch(
      as.vector(solve(M2, r2)),
      error = function(e) {
        stop(
          paste0(
            "M o M is singular or numerically non-invertible; ",
            "the Hadamard estimator cannot be computed for this design."
          ),
          call. = FALSE
        )
      }
    )

    kappa.had <- tryCatch(
      kappa(M2, exact = FALSE),
      error = function(e) Inf
    )

    if (!is.finite(kappa.had) || kappa.had > 1e10) {
      warning(
        paste0(
          "M o M is ill-conditioned. ",
          "The Hadamard variance estimate may be numerically unstable."
        ),
        call. = FALSE
      )
    }

    n.neg.had <- sum(vhat.had.raw < 0)
    n.above025.had <- sum(vhat.had.raw > 0.25)
    min.vhat.had <- min(vhat.had.raw)
    max.vhat.had <- max(vhat.had.raw)

    vhat.had <- vhat.had.raw

    if (var.method == "HAD" && n.neg.had > 0) {
      warning(
        paste0(
          n.neg.had,
          " element(s) of the unbiased observation-level variance estimate ",
          "are negative. This can occur with HAD. ",
          "HAD0 or HAD025 may be used as practical sensitivity analyses."
        ),
        call. = FALSE
      )
    }

    if (var.method == "HAD0") {
      vhat.had <- pmax(vhat.had.raw, 0)
    }

    if (var.method == "HAD025") {
      vhat.had <- pmin(pmax(vhat.had.raw, 0), 0.25)
    }

    ## X' diag(vhat) X.
    Meat <- crossprod(X, sweep(X, 1L, vhat.had, `*`))

    V2 <- XtX.inv %*% Meat %*% XtX.inv
    V2 <- (V2 + t(V2)) / 2

    dimnames(V2) <- list(names(coef1), names(coef1))
    Vout <- V2

    diag.V2 <- diag(V2)
    diag.tol <- 100 * .Machine$double.eps *
      max(1, max(abs(diag.V2), na.rm = TRUE))

    if (any(diag.V2 < -diag.tol, na.rm = TRUE)) {
      warning(
        paste0(
          "The estimated covariance matrix has a negative diagonal element. ",
          "The corresponding standard error, confidence interval, and P-value ",
          "are reported as NA."
        ),
        call. = FALSE
      )
    }

    se1 <- rep(NA_real_, length(diag.V2))
    names(se1) <- names(coef1)

    nonnegative <- diag.V2 >= -diag.tol
    se1[nonnegative] <- sqrt(pmax(diag.V2[nonnegative], 0))

  }

  cl1 <- coef1 - qnorm(cc) * se1
  cu1 <- coef1 + qnorm(cc) * se1

  Z <- coef1 / se1
  P <- 2 * pnorm(-abs(Z))

  res <- list(
    call = call,
    formula = formula,
    coefficients = coef1,
    se = se1,
    cl = cl1,
    cu = cu1,
    z = Z,
    p = P,
    eform = eform,
    cl.level = cl,
    digits = digits,
    var.method = var.method,
    vcov = Vout,
    model = gm1,
    n = n,
    vhat.had = vhat.had,
    vhat.had.raw = vhat.had.raw,
    kappa.had = kappa.had,
    n.neg.had = n.neg.had,
    n.above025.had = n.above025.had,
    min.vhat.had = min.vhat.had,
    max.vhat.had = max.vhat.had
  )

  class(res) <- "rqlm"
  res

}


print.rqlm <- function(x, digits = x$digits, ...) {

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat(
    "Coefficient estimates and CIs with robust SE estimator (",
    x$var.method,
    "):\n",
    sep = ""
  )

  est <- x$coefficients
  se <- x$se
  Lower <- x$cl
  Upper <- x$cu
  zval <- x$z
  pval <- x$p

  if (!isTRUE(x$eform)) {

    tab <- cbind(
      "Estimate" = est,
      "Robust SE" = se,
      "Lower" = Lower,
      "Upper" = Upper,
      "z value" = zval,
      "Pr(>|z|)" = pval
    )

  } else {

    tab <- cbind(
      "Estimate" = est,
      "Robust SE" = se,
      "exp(coef)" = exp(est),
      "Lower" = exp(Lower),
      "Upper" = exp(Upper),
      "z value" = zval,
      "Pr(>|z|)" = pval
    )

  }

  print(round(tab, digits))

  if (x$var.method %in% c("HAD", "HAD0", "HAD025")) {

    cat("\nHadamard-estimator diagnostics:\n")
    cat("  Condition number of M o M: ",
        signif(x$kappa.had, digits), "\n", sep = "")
    cat("  Negative raw variance estimates: ",
        x$n.neg.had, "\n", sep = "")
    cat("  Raw variance estimates above 1/4: ",
        x$n.above025.had, "\n", sep = "")
    cat("  Range of raw variance estimates: [",
        signif(x$min.vhat.had, digits), ", ",
        signif(x$max.vhat.had, digits), "]\n", sep = "")

    if (x$var.method == "HAD") {
      cat("  HAD uses the untruncated finite-sample unbiased estimator.\n")
    } else {
      cat(
        "  The truncation used by ",
        x$var.method,
        " does not preserve exact finite-sample unbiasedness.\n",
        sep = ""
      )
    }

  }

  invisible(x)

}


rqlm0 <- function(formula, data, family = poisson, eform = FALSE,
                  cl = 0.95, digits = 4) {

  gm1 <- glm(formula, data = data, family = family, x = TRUE)

  cc <- 1 - 0.5 * (1 - cl)

  coef1 <- gm1$coefficients
  V1 <- sandwich::sandwich(gm1)
  se1 <- sqrt(diag(V1))
  cl1 <- coef1 - qnorm(cc) * se1
  cu1 <- coef1 + qnorm(cc) * se1

  Z <- coef1 / se1
  P <- 2 * pnorm(-abs(Z))

  out <- data.frame(coef1, se1, cl1, cu1, P)
  colnames(out) <- c("coef", "SE", "CL", "CU", "P-value")

  if (isTRUE(eform)) {

    coef1 <- exp(coef1)
    cl1 <- exp(cl1)
    cu1 <- exp(cu1)

    out <- data.frame(coef1, se1, cl1, cu1, P)
    colnames(out) <- c("exp(coef)", "SE", "CL", "CU", "P-value")

  }

  round(out, digits)

}
