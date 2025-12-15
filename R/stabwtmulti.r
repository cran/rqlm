stabwtmulti <- function(formula, data, trunc=c(0.01,0.99), digits=4){

  data <- as.data.frame(data)
  
  call <- match.call()

  yname <- all.vars(formula)[1]
  
  ps.mod <- multinom(formula, data = data, trace=FALSE)
  ps <- predict(ps.mod, type = "probs")

  ps.mean <- apply(ps,2,mean)

  Q <- dim(ps)[2]
  N <- dim(ps)[1]

  treat <- factor(data[,yname])
  z <- levels(treat)

  sw0 <- ps

  for(k in 1:Q)  sw0[,k] <- ps.mean[k] / ps[,k]

  sw1 <- rep(NA,times=N)
  
  for(k in 1:Q)	sw1[treat==z[k]] <- sw0[treat==z[k],k]

  Q1 <- quantile(sw1,trunc[1])
  Q2 <- quantile(sw1,trunc[2])
  
  sw2 <- sw1
  sw2[sw2>Q2] <- Q2
  sw2[sw2<Q1] <- Q1


  ## オブジェクトとしてまとめて返す
  res <- list(
    call       = call,
    formula    = formula,
    sw1 		= sw1,
    sw2  		= sw2,
	Q1			= as.numeric(Q1),
	Q2			= as.numeric(Q2),
    digits      = digits,
	trunc		= trunc
  )
  class(res) <- "stabwt"
  return(res)
  
}


