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
	
	if(var.method=="GST"){
	
		Ainv <- vcov(gm1)

		p1 <- dim(Ainv)[1]
		
	    X  <- gm1$x
		pi <- gm1$fitted.values
		y  <- mdata$d
		w  <- pi * (1 - pi)
		
		did <- as.numeric(mdata$id)
		uid <- as.numeric(unique(mdata$id))
		
		BG1 <- BG2 <- matrix(numeric(p1*p1),p1)
		
		for(i in 1:n){
		
			wi <- which(did==uid[i])
			
			if(length(wi)==1){
			
				Ai <- pi[wi]*(1-pi[wi])*(X[wi,]%*%t(X[wi,]))
				Ui <- X[wi,]*(y[wi]-pi[wi])
				
				Aiih <- mat_inv_sqrt(Ai)
				BG1 <- BG1 + Aiih%*%Ui%*%t(Ui)%*%Aiih
			
			}

			if(length(wi)==2){
			
				w1 <- wi[1]
				w2 <- wi[2]
				Ai <- 2*pi[w1]*(1-pi[w1])*(X[w1,]%*%t(X[w1,]))
				Ui <- X[w1,]*(y[w1]-pi[w1]) + X[w2,]*(y[w2]-pi[w2])

				Aiih <- mat_inv_sqrt(Ai)
				BG1 <- BG1 + Aiih%*%Ui%*%t(Ui)%*%Aiih
			
			}
		
		}
		
		BG1 <- BG1/(n-p1)
		
		for(i in 1:n){
		
			wi <- which(did==uid[i])
			
			if(length(wi)==1){
			
				Ai <- pi[wi]*(1-pi[wi])*(X[wi,]%*%t(X[wi,]))
				Aih <- mat_sqrt(Ai)
				BG2 <- BG2 + Aih%*%BG1%*%Aih
			
			}

			if(length(wi)==2){
			
				w1 <- wi[1]
				w2 <- wi[2]
				Ai <- 2*pi[w1]*(1-pi[w1])*(X[w1,]%*%t(X[w1,]))
				Aih <- mat_sqrt(Ai)

				BG2 <- BG2 + Aih%*%BG1%*%Aih
			
			}
		
		}		
		
		V2 <- Ainv%*%BG2%*%Ainv
		se1 <- sqrt(diag(V2))
	
	}

	if(var.method=="WL"){
	
		Ainv <- vcov(gm1)

		p1 <- dim(Ainv)[1]
		
	    X  <- gm1$x
		pi <- gm1$fitted.values
		y  <- mdata$d
		w  <- pi * (1 - pi)
		
		did <- as.numeric(mdata$id)
		uid <- as.numeric(unique(mdata$id))
		
		BW1 <- BW2 <- matrix(numeric(p1*p1),p1)
		
		for(i in 1:n){
		
			wi <- which(did==uid[i])
			
			if(length(wi)==1){
			
				Ai <- pi[wi]*(1-pi[wi])*(X[wi,]%*%t(X[wi,]))
				Ui <- X[wi,]*(y[wi]-pi[wi])
				
				Aih <- mat_sqrt(Ai)
				Aiih <- mat_inv_sqrt(Ai)
				Hi <- Aih%*%Ainv%*%Aih
				Fi <- solve(diag(p1) - Hi)
				
				BW1 <- BW1 + Aiih%*%Fi%*%Ui%*%t(Ui)%*%Fi%*%Aiih
			
			}

			if(length(wi)==2){
			
				w1 <- wi[1]
				w2 <- wi[2]
				Ai <- 2*pi[w1]*(1-pi[w1])*(X[w1,]%*%t(X[w1,]))
				Ui <- X[w1,]*(y[w1]-pi[w1]) + X[w2,]*(y[w2]-pi[w2])

				Aih <- mat_sqrt(Ai)
				Aiih <- mat_inv_sqrt(Ai)
				Hi <- Aih%*%Ainv%*%Aih
				Fi <- solve(diag(p1) - Hi)
				
				BW1 <- BW1 + Aiih%*%Fi%*%Ui%*%t(Ui)%*%Fi%*%Aiih
			
			}
		
		}
		
		BW1 <- BW1/n
		
		for(i in 1:n){
		
			wi <- which(did==uid[i])
			
			if(length(wi)==1){
			
				Ai <- pi[wi]*(1-pi[wi])*(X[wi,]%*%t(X[wi,]))
				Aih <- mat_sqrt(Ai)
				BW2 <- BW2 + Aih%*%BW1%*%Aih
			
			}

			if(length(wi)==2){
			
				w1 <- wi[1]
				w2 <- wi[2]
				Ai <- 2*pi[w1]*(1-pi[w1])*(X[w1,]%*%t(X[w1,]))
				Aih <- mat_sqrt(Ai)

				BW2 <- BW2 + Aih%*%BW1%*%Aih
			
			}
		
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

mat_sqrt <- function(A, tol = 1e-12) {
  A <- (A + t(A)) / 2  # 片対称化（数値誤差の除去）
  eg <- eigen(A, symmetric = TRUE)
  vals <- eg$values
  vecs <- eg$vectors
  
  # 数値安定のため負の固有値はゼロに
  vals_pos <- pmax(vals, 0)
  sqrt_vals <- sqrt(vals_pos)
  
  return(vecs %*% diag(sqrt_vals) %*% t(vecs))
}

mat_inv_sqrt <- function(A, tol = 1e-12) {
  A <- (A + t(A)) / 2
  eg <- eigen(A, symmetric = TRUE)
  vals <- eg$values
  vecs <- eg$vectors
  
  vals_pos <- pmax(vals, 0)
  inv_sqrt_vals <- ifelse(vals_pos > tol, 1/sqrt(vals_pos), 0)

  return(vecs %*% diag(inv_sqrt_vals) %*% t(vecs))
}
