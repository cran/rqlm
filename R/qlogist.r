qlogist <- function(formula, data, eform=TRUE, cl=0.95, digits=4, var.method="MBN"){

	call <- match.call()

  yname <- all.vars(formula)[1]
  if (!(yname %in% names(data)))
    stop("Outcome variable not found in data.")

  data$id <- factor(seq_len(nrow(data)))

  cases     <- data[data[[yname]] == 1, , drop = FALSE]
  subcohort <- data                           # full cohort

  cases$d    <- 1
  subcohort$d <- 0
  
  n <- dim(data)[1]
  n1 <- dim(cases)[1]

  mdata <- rbind(cases, subcohort)

  mformula <- update(formula, d ~ .)
  gm1 <- glm(mformula, data = mdata, family = binomial("logit"), x=TRUE)

	cc <- 1 - 0.5*(1 - cl)
	
	coef1 <- gm1$coefficients

	if(var.method=="standard"){

		V1 <- vcovCL(gm1, cluster = mdata$id)
		se1 <- sqrt(diag(V1))
	
	}
	
	if(var.method=="MBN"){
	
		V1 <- vcovCL(gm1, cluster = mdata$id)
		Ainv <- vcov(gm1)
		A <- solve(Ainv)

		p1 <- dim(V1)[1]

		Q1 <- (n + n1 - 1)/(n + n1 - p1)
		Q2 <- n / (n - 1)
		
		Q3 <- sum(diag(V1%*%A))/p1
		
		delta <- min(0.5,p1/(n-p1))
		gamma <- max(1,Q3)
		
		V2 <- Q1*Q2*V1 + delta*gamma*Ainv
		se1 <- sqrt(diag(V2))
	
	}
	
	cl1 <- coef1 - qnorm(cc)*se1
	cu1 <- coef1 + qnorm(cc)*se1

	Z <- coef1/se1
	P <- 2*pnorm(-abs(Z))
	
	out <- data.frame(coef1,se1,cl1,cu1,P)
	colnames(out) <- c("coef","SE","CL","CU","P-value")

	if(eform==TRUE){
		
		coef1 <- exp(coef1)
		cl1 <- exp(cl1)
		cu1 <- exp(cu1)
		out <- data.frame(coef1,se1,cl1,cu1,P)
		colnames(out) <- c("exp(coef)","SE","CL","CU","P-value")
		
	}

	# out <- round(out, digits)
	
	# return(out)

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
    var.method  = var.method
  )
  class(res) <- "rqlm"
  return(res)
  
}
