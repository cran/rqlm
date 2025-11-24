rqlm <- function(formula, data, family=poisson, eform=FALSE, cl=0.95, digits=4, var.method="MBN"){

	call <- match.call()
	
	n <- dim(data)[1]

	gm1 <- glm(formula, data=data, family=family, x=TRUE, y=TRUE)
	
	cc <- 1 - 0.5*(1 - cl)
	
	coef1 <- gm1$coefficients

	if(var.method=="standard"){

		V1 <- sandwich(gm1)
		se1 <- sqrt(diag(V1))
	
	}
	
	if(var.method=="MBN"){
	
		V1 <- sandwich(gm1)
		Ainv <- vcov(gm1)
		A <- solve(Ainv)

		p1 <- dim(V1)[1]

		Q1 <- (n - 1)/(n - p1)
		Q2 <- n / (n - 1)
		
		Q3 <- sum(diag(V1%*%A))/p1
		
		delta <- min(0.5,p1/(n-p1))
		gamma <- max(1,Q3)
		
		V2 <- Q1*Q2*V1 + delta*gamma*Ainv
		se1 <- sqrt(diag(V2))
	
	}

	if(var.method=="GST"){
	
		Ainv <- vcov(gm1)

		p1 <- dim(Ainv)[1]
		
	    X  <- gm1$x
		pi <- gm1$fitted.values
		y  <- gm1$y
		w  <- pi * (1 - pi)
		
		BG1 <- BG2 <- matrix(numeric(p1*p1),p1)
		
		for(i in 1:n){
		
			Ai <- pi[i]*(1-pi[i])*(X[i,]%*%t(X[i,]))
			Ui <- X[i,]*(y[i]-pi[i])
				
			Aiih <- mat_inv_sqrt(Ai)
			BG1 <- BG1 + Aiih%*%Ui%*%t(Ui)%*%Aiih
			
		}
		
		BG1 <- BG1/(n-p1)
		
		for(i in 1:n){
		
			Ai <- pi[i]*(1-pi[i])*(X[i,]%*%t(X[i,]))
			Aih <- mat_sqrt(Ai)
			BG2 <- BG2 + Aih%*%BG1%*%Aih
		
		}		
		
		V2 <- Ainv%*%BG2%*%Ainv
		se1 <- sqrt(diag(V2))
	
	}

	if(var.method=="WL"){
	
		Ainv <- vcov(gm1)

		p1 <- dim(Ainv)[1]
		
	    X  <- gm1$x
		pi <- gm1$fitted.values
		y  <- gm1$y
		w  <- pi * (1 - pi)
		
		BW1 <- BW2 <- matrix(numeric(p1*p1),p1)
		
		for(i in 1:n){
		
			Ai <- pi[i]*(1-pi[i])*(X[i,]%*%t(X[i,]))
			Ui <- X[i,]*(y[i]-pi[i])
				
			Aih <- mat_sqrt(Ai)
			Aiih <- mat_inv_sqrt(Ai)
			Hi <- Aih%*%Ainv%*%Aih
			Fi <- solve(diag(p1) - Hi)
				
			BW1 <- BW1 + Aiih%*%Fi%*%Ui%*%t(Ui)%*%Fi%*%Aiih
		
		}
		
		BW1 <- BW1/n
		
		for(i in 1:n){
		
			Ai <- pi[i]*(1-pi[i])*(X[i,]%*%t(X[i,]))
			Aih <- mat_sqrt(Ai)
			BW2 <- BW2 + Aih%*%BW1%*%Aih
		
		}		
		
		V2 <- Ainv%*%BW2%*%Ainv
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
    var.method  = var.method
  )
  class(res) <- "rqlm"
  return(res)

}

print.rqlm <- function(x, digits = x$digits, ...) {

  cat("Call:\n")
  print(x$call)
  cat("\n")

  ## ラベル（CRとかMBNとか）を表示
  cat("Coefficients with robust variance (", x$var.method, "):\n", sep = "")

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

rqlm0 <- function(formula, data, family=poisson, eform=FALSE, cl=0.95, digits=4){

	gm1 <- glm(formula, data=data, family=family, x=TRUE)
	
	cc <- 1 - 0.5*(1 - cl)
	
	coef1 <- gm1$coefficients
	V1 <- sandwich(gm1)
	se1 <- sqrt(diag(V1))
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

	out <- round(out, digits)
	
	return(out)

}
