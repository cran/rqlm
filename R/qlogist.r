qlogist <- function(formula, data, eform=TRUE, cl=0.95, digits=4, var.method="MBN", id=NULL){

  call <- match.call()

  var.method <- match.arg(var.method, c("standard", "MBN", "GST", "WL"))

  mf <- stats::model.frame(formula, data=data, drop.unused.levels=TRUE)
  mt <- attr(mf, "terms")
  y <- stats::model.response(mf)
  if (!is.numeric(y) && !is.logical(y))
    stop("The outcome must be coded as numeric or logical 0/1.")
  if (!is.null(dim(y)) || anyNA(y) || !all(y %in% c(0, 1)))
    stop("The outcome must be a vector coded as 0/1 without missing values after model-frame processing.")

  n <- nrow(mf)
  id_vec <- .rqlm_id(data, substitute(id), attr(mf, "na.action"), n)
  if (!is.null(id_vec) && !(var.method %in% c("standard", "MBN")))
    stop("With id supplied, var.method must be \"standard\" or \"MBN\".")
  if (is.null(id_vec)) id_vec <- seq_len(n)
  K <- length(unique(id_vec))
  if (K < 2L) stop("At least two independent clusters are required.")

  cases <- which(y == 1)
  n1 <- length(cases)
  rows <- c(cases, seq_len(n))
  mdata <- mf[rows, , drop=FALSE]
  mdata[[1L]] <- c(rep.int(1, n1), rep.int(0, n))
  attr(mdata, "terms") <- mt
  attr(mdata, "na.action") <- NULL
  mid <- factor(id_vec[rows])
  X <- stats::model.matrix(mt, mf)
  off <- stats::model.offset(mf)
  if (!is.null(off)) off <- off[rows]
  yy <- mdata[[1L]]
  names(yy) <- rownames(mdata)

  Xm <- X[rows, , drop=FALSE]
  rownames(Xm) <- rownames(mdata)
  gm1 <- stats::glm.fit(Xm, yy,
    family=binomial("logit"), offset=off, intercept=attr(mt, "intercept") > 0L)
  gm1$model <- mdata
  gm1$x <- Xm
  gm1$terms <- mt
  gm1$formula <- formula
  gm1$call <- quote(glm(formula, data=mdata, family=binomial("logit"), x=TRUE))
  gm1$offset <- off
  gm1$control <- stats::glm.control()
  gm1$method <- "glm.fit"
  gm1$contrasts <- attr(X, "contrasts")
  gm1$xlevels <- stats::.getXlevels(mt, mf)
  class(gm1) <- c("glm", "lm")
  Vout <- NULL

	cc <- 1 - 0.5*(1 - cl)
	
	coef1 <- gm1$coefficients

	if(var.method=="standard"){

		V1 <- vcovCL(gm1, cluster = mid)
		Vout <- V1
		se1 <- sqrt(diag(V1))
	
	}
	
	if(var.method=="MBN"){
	
		V1 <- vcovCL(gm1, cluster = mid)
		Ainv <- vcov(gm1)
		A <- solve(Ainv)

		p1 <- dim(V1)[1]

		if (K <= p1 || n + n1 <= p1)
		  stop("The number of independent clusters must be larger than the number of model parameters for var.method = \"MBN\".")

		Q1 <- (n + n1 - 1)/(n + n1 - p1)
		Q2 <- K / (K - 1)
		
		Q3 <- sum(diag(V1%*%A))/p1
		
		delta <- min(0.5,p1/(K-p1))
		gamma <- max(1,Q3)
		
		V2 <- Q1*Q2*V1 + delta*gamma*Ainv
		Vout <- V2
		se1 <- sqrt(diag(V2))
	
	}
	
	if(var.method=="GST"){
	
		Ainv <- vcov(gm1)

		p1 <- dim(Ainv)[1]
		
	    X  <- gm1$x
		pi <- gm1$fitted.values
		y  <- gm1$y
		w  <- pi * (1 - pi)
		
		did <- as.numeric(mid)
		uid <- as.numeric(unique(mid))
		
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
		Vout <- V2
		se1 <- sqrt(diag(V2))
	
	}

	if(var.method=="WL"){
	
		Ainv <- vcov(gm1)

		p1 <- dim(Ainv)[1]
		
	    X  <- gm1$x
		pi <- gm1$fitted.values
		y  <- gm1$y
		w  <- pi * (1 - pi)
		
		did <- as.numeric(mid)
		uid <- as.numeric(unique(mid))
		
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
		Vout <- V2
		se1 <- sqrt(diag(V2))
	
	}
	
  dimnames(Vout) <- list(names(coef1), names(coef1))

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
    var.method  = var.method,
    vcov        = Vout,
    model       = gm1,
    n           = n,
    n.clusters  = K
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
