rqlm <- function(formula, data, family=poisson, eform=FALSE, cl=0.95, digits=4){

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
