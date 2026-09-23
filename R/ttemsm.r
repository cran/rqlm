ttemsm <- function(formula, data, id, weight=NULL, family=quasibinomial(link="cloglog"), eform=TRUE, cl=0.95, digits=4, var.method="MBN"){

  data <- as.data.frame(data)
  
  call <- match.call()
  var.method <- match.arg(var.method, c("standard", "MBN"))

  yname <- all.vars(formula)[1]
  if (!(yname %in% names(data)))
    stop("Outcome variable not found in data.")

  id_vec <- .rqlm_id(data, substitute(id))
  if (is.null(id_vec)) stop("id must specify a column in data.")
  wt_name <- substitute(weight)
  if (is.null(wt_name)) {
    data$w_vec <- rep(1, nrow(data))
  } else {
    if (is.symbol(wt_name)) wt_name <- as.character(wt_name)
    if (!is.character(wt_name) || length(wt_name) != 1L ||
        !(wt_name %in% names(data)))
      stop("weight must be NULL or a column name in data.")
    data$w_vec <- data[[wt_name]]
  }
  
  gm1 <- glm(formula, data = data, family = family, weights=w_vec)		# cloglogリンクで、離散Cox回帰（HRを推定したいという場合）。リスク差・累積発生率（g-formula / 標準化）の推定が目的なら、logitリンクに変える。

  L <- length(gm1$y)
  id_vec <- .rqlm_id(data, substitute(id), gm1$na.action, L)
  K <- length(unique(id_vec))
  if (K < 2L) stop("At least two independent clusters are required.")
  Vout <- NULL

	cc <- 1 - 0.5*(1 - cl)
	
	coef1 <- gm1$coefficients

	if(var.method=="standard"){

		V1 <- vcovCL(gm1, cluster = id_vec)
		Vout <- V1
		se1 <- sqrt(diag(V1))

	}
	
	if(var.method=="MBN"){
	
		L <- length(gm1$y)
		K <- length(unique(id_vec))
	
		V1 <- vcovCL(gm1, cluster = id_vec)
		Ainv <- vcov(gm1)
		A <- solve(Ainv)

		p1 <- dim(V1)[1]

		if (K <= p1 || L <= p1)
		  stop("The number of independent clusters must be larger than the number of model parameters for var.method = \"MBN\".")

		Q1 <- (L - 1)/(L - p1)
		Q2 <- K / (K - 1)
		
		Q3 <- sum(diag(V1%*%A))/p1
		
		delta <- min(0.5,p1/(K-p1))
		gamma <- max(1,Q3)
		
		V2 <- Q1*Q2*V1 + delta*gamma*Ainv
		Vout <- V2
		se1 <- sqrt(diag(V2))
	
	}
	
	cl1 <- coef1 - qnorm(cc)*se1
	cu1 <- coef1 + qnorm(cc)*se1

	Z <- coef1/se1
	P <- 2*pnorm(-abs(Z))
	

  ## オブジェクトとしてまとめて返す
  res <- list(
    call       = call,
    formula    = formula,
    coefficients = coef1,
    se          = se1,
    cl          = cl1,
    cu          = cu1,
    z           = Z,
    p           = P,
    eform       = eform,
    cl.level    = cl,
    digits      = digits,
    var.method  = var.method,
    vcov        = Vout,
    model       = gm1,
    n           = L,
    n.clusters  = K
  )
  class(res) <- "ttemsm"
  return(res)
  
}



if (getRversion() >= "2.15.1") {
  utils::globalVariables(c("w_vec"))
}



print.ttemsm <- function(x, digits = x$digits, ...) {

  cat("Call:\n")
  print(x$call)
  cat("\n")

  ## ラベル（CRとかMBNとか）を表示
  cat("Coefficient estimates and CIs with cluster-robust SE estimator (", x$var.method, "):\n", sep = "")

  ## ログスケールの推定値・SE・z・p
  est  <- x$coefficients
  se   <- x$se
  Lower <- x$cl
  Upper <- x$cu
  
  zval <- x$z
  pval <- x$p

 if(x$eform==FALSE){

  tab <- cbind(
    "Estimate" = est,
    "Robust SE" = se,
	"Lower" = Lower,
	"Upper" = Upper,
    "z value"   = zval,
    "Pr(>|z|)"  = pval
  )

 }

 if(x$eform==TRUE){

  tab <- cbind(
    "Estimate" = est,
    "Robust SE" = se,
    "exp(coef)"   = exp(est),
    "Lower"   = exp(Lower),
    "Upper"   = exp(Upper),
    "z value"   = zval,
    "Pr(>|z|)"  = pval
  )

 }

  tab <- round(tab, digits)
  print(tab)
  invisible(x)
  
}
