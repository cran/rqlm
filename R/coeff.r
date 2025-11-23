coeff <- function(gm, eform=FALSE, cl=0.95, digits=4){

	cc <- 1 - 0.5*(1 - cl)

	modegm <- mode(gm)

	if(modegm=="list"){		# lmかglmの出力オブジェクト

		call <- gm$call

		sgm <- summary(gm)
		cgm <- sgm$coefficients
		
		if(is.null(sgm$sigma)){		# glmの場合

			est <- cgm[,1]
			se <- cgm[,2]
			cl1 <- cgm[,1] - qnorm(cc)*cgm[,2]
			cl2 <- cgm[,1] + qnorm(cc)*cgm[,2]
			P <- 2*pnorm(-abs(cgm[,1]/cgm[,2]))
	
			res <- list(
				call       = call,
				est = est,
				se  = se,
				cl  = cl1,
				cu  = cl2,
				p   = P,
				eform = eform,
				digits = digits
				)
			class(res) <- "coeff"
			return(res)

		}
		
		if(is.null(sgm$sigma)==FALSE){		# lmの場合

			est <- cgm[-1,1]
			se <- cgm[-1,2]
			cl1 <- cgm[-1,1] - qt(cc,df=sgm$df)*cgm[-1,2]
			cl2 <- cgm[-1,1] + qt(cc,df=sgm$df)*cgm[-1,2]
			P <- 2*pt(-abs(cgm[-1,1]/cgm[-1,2]),df=sgm$df)
	
			res <- list(
				call = call,
				est = est,
				se  = se,
				cl  = cl1,
				cu  = cl2,
				p   = P,
				eform = eform,
				digits = digits
				)
			class(res) <- "coeff"
			return(res)

		}
	
	}

	if(modegm=="S4"){		# lmerかglmerの出力オブジェクト

		sgm <- summary(gm)

		call <- sgm$call
		cgm <- sgm$coefficients
		
		est <- cgm[,1]
		se <- cgm[,2]
		cl1 <- cgm[,1] - qnorm(cc)*cgm[,2]
		cl2 <- cgm[,1] + qnorm(cc)*cgm[,2]
		P <- 2*pnorm(-abs(cgm[,1]/cgm[,2]))
	
			res <- list(
				call = call,
				est = est,
				se  = se,
				cl  = cl1,
				cu  = cl2,
				p   = P,
				eform = eform,
				digits = digits
				)
			class(res) <- "coeff"
			return(res)

	}
	
}


print.coeff <- function(x, digits = x$digits, ...) {

  cat("Call:\n")
  print(x$call)
  cat("\n")

  cat("Regression coefficients with confidence intervals:\n", sep = "")

  ## ログスケールの推定値・SE・z・p
  est  <- x$est
  se   <- x$se
  Lower <- x$cl
  Upper <- x$cu
  pval <- x$p

 if(x$eform==FALSE){

  tab <- cbind(
    "Estimate" = est,
    "SE" = se,
	"Lower" = Lower,
	"Upper" = Upper,
    "Pr(>|z|)"  = pval
  )

 }

 if(x$eform==TRUE){

  tab <- cbind(
    "Estimate" = est,
    "SE" = se,
    "exp(coef)"   = exp(est),
    "Lower"   = exp(Lower),
    "Upper"   = exp(Upper),
    "Pr(>|z|)"  = pval
  )

 }

  tab <- round(tab, digits)
  print(tab)
  invisible(x)
  
}
