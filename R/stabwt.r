stabwt <- function(formula, data, trunc=c(0.01,0.99), digits=4){

  data <- as.data.frame(data)
  
  call <- match.call()

  yname <- all.vars(formula)[1]

  ps.mod <- glm(formula, family = binomial, data = data)
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





print.stabwt <- function(x, digits = x$digits, ...) {

  cat("Call:\n")
  print(x$call)
  cat("\n")

  sum1 <- round(summary(x$sw1),digits)
  sum2 <- round(summary(x$sw2),digits)
  Q <- round(c(x$Q1,x$Q2),digits)

  ## ラベル（CRとかMBNとか）を表示
  cat("Summary of the stabilized weights (untruncated):\n", sep = "")
  print(sum1)
  cat("\n")

  cat("Summary of the stabilized weights (truncated):\n", sep = "")
  print(sum2)
  cat("\n")

  cat("Truncated quantiles: ", x$trunc[1], " and ", x$trunc[2], "\n", sep = "")
  
  invisible(x)

}


