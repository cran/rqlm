mi_glm <- function(ice, formula, family=gaussian, offset=NULL, eform=FALSE, cl=0.95, digits=4){

	B <- ice$m

	cc <- 1 - 0.5*(1 - cl)
	
	data1 <- complete(ice,1)
	gm1 <- glm(formula, data=data1, family=family, x=TRUE)
	coef1 <- gm1$coefficients

	n <- dim(data1)[1]

	p <- length(coef1)
	
	beta.b <- matrix(numeric(p*B),B)
	V.b <- matrix(numeric(p*p),p)

	for(b in 1:B){
	
		data.b <- complete(ice,b)

		if(is.null(offset)) gm.b <- glm(formula, data=data.b, family=family, x=TRUE)
		if(is.null(offset)==FALSE) gm.b <- glm(formula, data=data.b, family=family, offset=offset, x=TRUE)
	
		beta.b[b,] <- gm.b$coefficients
		V.b <- V.b + vcov(gm.b)

	}

	coef1 <- apply(beta.b,2,mean)

	W1 <- V.b/B
	B1 <- cov(beta.b)

	V1 <- W1 + (1+1/B)*B1

	##
	
	r <- rep(NA,times=p)
	for(j in 1:p)	r[j] <- (1+1/B)*B1[j,j]/W1[j,j]
	dof <- (B-1)*((1+1/r)^2)

	##

	se1 <- sqrt(diag(V1))
	cl1 <- coef1 - qt(cc,df=dof)*se1
	cu1 <- coef1 + qt(cc,df=dof)*se1

	Z <- coef1/se1
	P <- 2*pt(-abs(Z),df=dof)
	
	out <- data.frame(coef1,se1,cl1,cu1,dof,P)
	colnames(out) <- c("coef","SE","CL","CU","df","P-value")

	if(eform==TRUE){
		
		coef1 <- exp(coef1)
		cl1 <- exp(cl1)
		cu1 <- exp(cu1)
		out <- data.frame(coef1,se1,cl1,cu1,dof,P)
		colnames(out) <- c("exp(coef)","SE","CL","CU","df","P-value")
		
	}

	out <- round(out, digits)
	
	return(out)

}
