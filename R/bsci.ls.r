bsci.ls <- function(formula, data, x.name=NULL, B=1000, cl=0.95, C0=10^-5, digits=4, seed=527916){

	set.seed(seed)

	a0 <- 1 - cl		# significant level

	qm1 <- rqlm(formula, data=data, family=gaussian, cl=cl, digits=digits)
	qm2 <- rqlm(formula, data=data, family=gaussian, cl=cl, digits=999999)

	wx <- which(rownames(qm2)==x.name)

	e1 <- qm2[wx,1]
	se1 <- qm2[wx,2]
	cl1 <- qm2[wx,3]
	cu1 <- qm2[wx,4]

	# Exploring lower limits
	
	P0 <- bs.ls(formula, data=data, x.name, beta0=e1, edfplot=FALSE, B=B)
	P1 <- bs.ls(formula, data=data, x.name, beta0=cl1, edfplot=FALSE, B=B)

	U1	<-	e1		# 上限の初期値

	if(P1 < a0)	L1 <- cl1		# 下限の初期値

	if(P1 >= a0){
	
		a1 <- 1
	
		while(P1 >= a0){
		
			a2 <- cl1 - a1*se1
			
			P1 <- bs.ls(formula, data=data, x.name, beta0=a2, edfplot=FALSE, B=B)

			a1 <- a1 + 1

			if(a1 >= 20) return("Error: Singular problems occur, please contact the authors.")
			
		}
		
		L1 <- a2

	}

	D1 <- U1 - L1

	message("The computationl process 1/4 was completed.")

	while(D1 > C0){

		M1 <- (L1 + U1)/2

		P2 <- bs.ls(formula, data=data, x.name, beta0=M1, edfplot=FALSE, B=B)
		
		if(P2 >= a0) U1 <- M1
		if(P2 < a0) L1 <- M1
		
		D1 <- U1 - L1

	}
	
	CL <- M1		# 下側の信頼限界を確定

	message("The computationl process 2/4 was completed.")

	# Exploring upper limits
	
	P0 <- bs.ls(formula, data=data, x.name, beta0=e1, edfplot=FALSE, B=B)
	P1 <- bs.ls(formula, data=data, x.name, beta0=cu1, edfplot=FALSE, B=B)

	L1	<-	e1		# 下限の初期値

	if(P1 < a0)	U1 <- cu1		# 上限の初期値

	if(P1 >= a0){
	
		a1 <- 1
	
		while(P1 >= a0){
		
			a2 <- cu1 + a1*se1
			
			P1 <- bs.ls(formula, data=data, x.name, beta0=a2, edfplot=FALSE, B=B)

			a1 <- a1 + 1

			if(a1 >= 20) return("Error: Singular problems occur, please contact the authors.")
		
		}
		
		U1 <- a2

	}

	D1 <- U1 - L1

	message("The computationl process 3/4 was completed.")

	while(D1 > C0){

		M1 <- (L1 + U1)/2

		P2 <- bs.ls(formula, data=data, x.name, beta0=M1, edfplot=FALSE, B=B)
		
		if(P2 >= a0) L1 <- M1
		if(P2 < a0) U1 <- M1
		
		D1 <- U1 - L1
		
	}
	
	CU <- M1		# 上側の信頼限界を確定

	message("The computationl process 4/4 was completed.")
	message(" ")

	out2 <- round(c(CL,CU),digits)
	
	P2 <- round(bs.ls(formula, data=data, x.name, beta0=0, edfplot=FALSE, B=B),digits)

	R1 <- list("Modified least-squares regression with the Wald-type approximation"=qm1, "Bootstrap confidence interval for the corresponding covariate"=out2,
			"Bootstrap P-value"=P2)

	return(R1)

}


