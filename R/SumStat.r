SumStat <- function(formula, data, trunc=c(0.01,0.99), digits=3){

  data <- as.data.frame(data)
  
  call <- match.call()

  yname <- all.vars(formula)[1]

  ps.mod <- glm(formula, family = binomial, data = data, y=TRUE, x=TRUE)
  ps <- predict(ps.mod, type = "response")
	
  pA  <- mean(ps.mod$y == 1)
  pC  <- 1 - pA
	
  sw1 <- ifelse(ps.mod$y == 1,
                 pA / ps,
                 pC / (1 - ps))

  Q1 <- quantile(sw1,trunc[1])
  Q2 <- quantile(sw1,trunc[2])
  
  sw2 <- sw1
  sw2[sw2>Q2] <- Q2
  sw2[sw2<Q1] <- Q1
  
  ##
  
  treat <- as.numeric(factor(ps.mod$y)) - 1
  X <- ps.mod$x[,-1]
  
  X0 <- X[treat==0,]
  X1 <- X[treat==1,]
  p <- dim(X)[2]
  
  n0 <- dim(X0)[1]
  n1 <- dim(X1)[1]
  
  w0 <- sw2[treat==0]
  w1 <- sw2[treat==1]

  ##
    
  R <- NULL
  
  for(k in 1:p){

	x0 <- X0[,k]
	x1 <- X1[,k]
		
	l0 <- length(unique(x0))
	l1 <- length(unique(x1))
		
	if((l0!=2)||(l1!=2)){		# 2値変数でない（連続変数である）場合

		m0 <- mean(X0[,k],na.rm=TRUE)
		m1 <- mean(X1[,k],na.rm=TRUE)
		s0 <- sd(X0[,k],na.rm=TRUE)
		s1 <- sd(X1[,k],na.rm=TRUE)
		cs0 <- sqrt( (s0*s0 + s1*s1)/2 )		# 2で割る。サンプルサイズで重みづけしない。Austinらもそう言っているとのこと。
		smd0 <- abs(m1 - m0)/cs0
	
		wm0 <- sum(w0*X0[,k])/sum(w0)
		wm1 <- sum(w1*X1[,k])/sum(w1)
		ws0 <- sqrt( sum( w0*( (X0[,k] - wm0)^2 ) )/sum(w0) )
		ws1 <- sqrt( sum( w1*( (X1[,k] - wm1)^2 ) )/sum(w1) )
		smd1 <- abs(wm1 - wm0)/cs0		# SMDの分母は、未重み付けSDにするらしい
		
		type <- c("conti")
		
		Rk <- data.frame(m0,m1,s0,s1,smd0,wm0,wm1,ws0,ws1,smd1,type)
		R <- rbind(R,Rk)
	
	}
  
	if((l0==2)&&(l1==2)){		# 2値変数である場合

		x0 <- as.numeric(factor(x0)) - 1
		x1 <- as.numeric(factor(x1)) - 1

		m0 <- mean(x0,na.rm=TRUE)
		m1 <- mean(x1,na.rm=TRUE)
		s0 <- sqrt( m0*(1-m0) )
		s1 <- sqrt( m1*(1-m1) )
		cs0 <- sqrt( (s0*s0 + s1*s1)/2 )
		smd0 <- abs(m1 - m0)/cs0
	
		wm0 <- sum(w0*x0)/sum(w0)
		wm1 <- sum(w1*x1)/sum(w1)
		ws0 <- sqrt( wm0*(1-wm0) )
		ws1 <- sqrt( wm1*(1-wm1) )
		smd1 <- abs(wm1 - wm0)/cs0

		type <- c("binary")
		
		Rk <- data.frame(m0,m1,s0,s1,smd0,wm0,wm1,ws0,ws1,smd1,type)
		R <- rbind(R,Rk)
	
	}
  
  }
  
  R[,1:10] <- round(R[,1:10],digits)

  rownames(R) <- colnames(X)
  colnames(R) <- c("mean0","mean1","sd0","sd1","SMD","wmean0","wmean1","wsd0","wsd1","wSMD","type")
  
  return(R)
  
}



