coeff <- function(gm, eform=FALSE, cl=0.95, digits=4){

	cc <- 1 - 0.5*(1 - cl)

	modegm <- mode(gm)

	if(modegm=="list"){		# lmかglmの出力オブジェクト

		sgm <- summary(gm)
		cgm <- sgm$coefficients
		
		if(is.null(sgm$sigma)){		# glmの場合

			est <- cgm[,1]
			se <- cgm[,2]
			cl1 <- cgm[,1] - qnorm(cc)*cgm[,2]
			cl2 <- cgm[,1] + qnorm(cc)*cgm[,2]
			P <- 2*pnorm(-abs(cgm[,1]/cgm[,2]))
	
			if(eform==FALSE){
	
				R <- data.frame(est,se,cl1,cl2,P)
				colnames(R) <- c("coef","SE","CL","CU","P-value")
				R <- round(R,digits)
				return(R)
	
			}

			if(eform==TRUE){
	
				est <- exp(est)
				cl1 <- exp(cl1)
				cl2 <- exp(cl2)
	
				R <- data.frame(est,se,cl1,cl2,P)
				colnames(R) <- c("exp(coef)","SE","CL","CU","P-value")
				R <- round(R,digits)
				return(R)
	
			}	

		}
		
		if(is.null(sgm$sigma)==FALSE){		# lmの場合

			est <- cgm[-1,1]
			se <- cgm[-1,2]
			cl1 <- cgm[-1,1] - qt(cc,df=sgm$df)*cgm[-1,2]
			cl2 <- cgm[-1,1] + qt(cc,df=sgm$df)*cgm[-1,2]
			P <- 2*pt(-abs(cgm[-1,1]/cgm[-1,2]),df=sgm$df)
	
			if(eform==FALSE){
	
				R <- data.frame(est,se,cl1,cl2,P)
				colnames(R) <- c("coef","SE","CL","CU","P-value")
				R <- round(R,digits)
				return(R)
	
			}

			if(eform==TRUE){
	
				est <- exp(est)
				cl1 <- exp(cl1)
				cl2 <- exp(cl2)
	
				R <- data.frame(est,se,cl1,cl2,P)
				colnames(R) <- c("exp(coef)","SE","CL","CU","P-value")
				R <- round(R,digits)
				return(R)
	
			}	

		}
	
	}

	if(modegm=="S4"){		# lmerかglmerの出力オブジェクト

		sgm <- summary(gm)
		cgm <- sgm$coefficients
		
		est <- cgm[,1]
		se <- cgm[,2]
		cl1 <- cgm[,1] - qnorm(cc)*cgm[,2]
		cl2 <- cgm[,1] + qnorm(cc)*cgm[,2]
		P <- 2*pnorm(-abs(cgm[,1]/cgm[,2]))
	
		if(eform==FALSE){
	
			R <- data.frame(est,se,cl1,cl2,P)
			colnames(R) <- c("coef","SE","CL","CU","P-value")
			R <- round(R,digits)
			return(R)
	
		}

		if(eform==TRUE){
	
			est <- exp(est)
			cl1 <- exp(cl1)
			cl2 <- exp(cl2)
	
			R <- data.frame(est,se,cl1,cl2,P)
			colnames(R) <- c("exp(coef)","SE","CL","CU","P-value")
			R <- round(R,digits)
			return(R)
	
		}	

	}
	
}
