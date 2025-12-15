stabwtlong <- function(formula_denom, formula_num, data, trunc=c(0.01,0.99), digits=4){

	data <- as.data.frame(data)
  
	call <- match.call()

	# --- denominator model ---
	fit.denom <- glm(formula_denom, data = data, family = binomial())

	# --- numerator model ---
	fit.num <- glm(formula_num, data = data, family = binomial())

	# --- predicted probabilities ---
	p_denom <- predict(fit.denom, type="response")
	p_num   <- predict(fit.num, type="response")

	# --- stabilized weight ---
	sw1 <- ((1 - p_num)/(1 - p_denom))^((1 - fit.num$y)) * (p_num/p_denom)^(fit.num$y)

	Q1 <- quantile(sw1,trunc[1])
	Q2 <- quantile(sw1,trunc[2])
  
	sw2 <- sw1
	sw2[sw2>Q2] <- Q2
	sw2[sw2<Q1] <- Q1


  ## オブジェクトとしてまとめて返す
  res <- list(
    call       = call,
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


