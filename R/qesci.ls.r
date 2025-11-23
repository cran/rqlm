qes.ls <- function(formula, data, x.name=NULL, beta0=0){

	gm1 <- lm(formula, data=data, y=TRUE, x=TRUE)
	
	X <- gm1$x
	y <- gm1$y

	cn <- colnames(X)

	w0 <- which(cn==x.name)

	X1 <- X[,-c(1,w0)]
	x0 <- X[,w0]
	
	dt1 <- data.frame(y,X1)
	
	z1 <- beta0*x0
	
	gm2 <- lm(y ~ . , data=dt1, offset=z1, x=TRUE, y=TRUE)

	mu3 <- predict(gm2)

	U3 <- (y - mu3) %*% X

	mu3[mu3<0] <- 0
	mu3[mu3>1] <- 1
	
	dg3 <- diag(mu3*(1 - mu3))
	V3 <- t(X) %*% dg3 %*% X

	T <- as.numeric(U3%*%ginv(V3)%*%t(U3))			# quasi-score test statistic

	P <- as.numeric(1 - pchisq(T,df=1))				# P-value

	R1 <- list(T=T,P=P,y=y,X=X,X1=X1,dt1=dt1,z1=z1,mu=mu3)

	return(R1)

}


bs.ls <- function(formula, data, x.name=NULL, beta0=0, edfplot=TRUE, B=1000){

	data <- data.frame(data)
	
	n <- dim(data)[1]

	Q1 <- qes.ls(formula, data, x.name=x.name, beta0)

	T <- Q1$T

	X1 <- Q1$X1
	X <- Q1$X

	Sr <- numeric(B)

	for(r in 1:B){
	
		yr <- rbinom(n, size=1, prob=Q1$mu)

		dt1 <- data.frame(yr,X1)
	
		lmr <- lm(yr ~ . , data=dt1, offset=Q1$z1)

		mur <- predict(lmr)		# きちんとオフセットを足しあげた上での結果変数の予測値であることを確認済み。
	
		Ur <- (yr - mur) %*% X		# quasi-score; 検算済

		mur[mur<0] <- 0
		mur[mur>1] <- 1
	
		dgr <- diag(mur*(1 - mur))
		Vr <- t(X) %*% dgr %*% X

		Sr[r] <- Ur%*%ginv(Vr)%*%t(Ur)		# quasi-score test statistic
	
	}
	
	ecdf1 <- ecdf(Sr)
	
	P <- 1 - ecdf1(T)

	if(edfplot==TRUE){
	
		plot(ecdf1,main="Bootstrap distribution of the quasi-score statistic")
		abline(0.95, 0, lty=2, col="red")
		Q95 <- round(quantile(ecdf1,0.95),4)
		title(sub=paste0("95th percentile: ",Q95))
	
	}

	return(P)

}

qesci.ls <- function(formula, data, x.name=NULL, cl=0.95, C0=10^-5, digits=4){

	a0 <- 1 - cl		# significant level

	qm1 <- rqlm0(formula, data=data, family=gaussian, cl=cl, digits=digits)
	qm2 <- rqlm0(formula, data=data, family=gaussian, cl=cl, digits=999999)

	wx <- which(rownames(qm2)==x.name)

	e1 <- qm2[wx,1]
	se1 <- qm2[wx,2]
	cl1 <- qm2[wx,3]
	cu1 <- qm2[wx,4]

	# Exploring lower limits
	
	P0 <- qes.ls(formula, data=data, x.name, beta0=e1)$P
	P1 <- qes.ls(formula, data=data, x.name, beta0=cl1)$P

	U1	<-	e1		# 上限の初期値

	if(P1 < a0)	L1 <- cl1		# 下限の初期値

	if(P1 >= a0){
	
		a1 <- 1
	
		while(P1 >= a0){
		
			a2 <- cl1 - a1*se1
			
			P1 <- qes.ls(formula, data=data, x.name, beta0=a2)$P

			a1 <- a1 + 1

			if(a1 >= 20) return("Error: Singular problems occur, possibly by the separation issue. Please try bsci.ls.")
			
		}
		
		L1 <- a2

	}

	D1 <- U1 - L1

	while(D1 > C0){

		M1 <- (L1 + U1)/2

		P2 <- qes.ls(formula, data=data, x.name, beta0=M1)$P
		
		if(P2 >= a0) U1 <- M1
		if(P2 < a0) L1 <- M1
		
		D1 <- U1 - L1

	}
	
	CL <- M1		# 下側の信頼限界を確定

	# Exploring upper limits
	
	P0 <- qes.ls(formula, data=data, x.name, beta0=e1)$P
	P1 <- qes.ls(formula, data=data, x.name, beta0=cu1)$P

	L1	<-	e1		# 下限の初期値

	if(P1 < a0)	U1 <- cu1		# 上限の初期値

	if(P1 >= a0){
	
		a1 <- 1
	
		while(P1 >= a0){
		
			a2 <- cu1 + a1*se1
			
			P1 <- qes.ls(formula, data=data, x.name, beta0=a2)$P

			a1 <- a1 + 1

			if(a1 >= 20) return("Error: Singular problems occur, possibly by the separation issue. Please try bsci.ls.")
			
		}
		
		U1 <- a2

	}

	D1 <- U1 - L1

	while(D1 > C0){

		M1 <- (L1 + U1)/2

		P2 <- qes.ls(formula, data=data, x.name, beta0=M1)$P
		
		if(P2 >= a0) L1 <- M1
		if(P2 < a0) U1 <- M1
		
		D1 <- U1 - L1
		
	}
	
	CU <- M1		# 上側の信頼限界を確定

	out2 <- round(c(CL,CU),digits)

	P2 <- round(qes.ls(formula, data=data, x.name, beta0=0)$P,digits)

	R1 <- list("Modified least-squares regression with the Wald-type approximation"=qm1, "Quasi-score confidence interval for the corresponding covariate"=out2,
			"P-value for the quasi-score test"=P2)

	return(R1)

}

