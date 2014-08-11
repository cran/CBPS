CBPS <- function(formula, data, na.action, ATT=NULL, method="over",iterations=NULL, standardize = TRUE, twostep = TRUE, ...) {
		if (missing(data)) 
			data <- environment(formula)
		call <- match.call()
		family <- binomial()

		mf <- match.call(expand.dots = FALSE)
		m <- match(c("formula", "data", "na.action"), names(mf), 0L)
		mf <- mf[c(1L, m)]
		mf$drop.unused.levels <- TRUE
		mf[[1L]] <- as.name("model.frame")
		
			mf <- eval(mf, parent.frame())
			mt <- attr(mf, "terms")
			Y <- model.response(mf, "any")
			if (length(dim(Y)) == 1L) {
					nm <- rownames(Y)
					dim(Y) <- NULL
					if (!is.null(nm)) 
						names(Y) <- nm
			}
			  
			X <- if (!is.empty.model(mt)) model.matrix(mt, mf)#[,-2]
			else matrix(, NROW(Y), 0L)
				
			X<-cbind(1,X[,apply(X,2,sd)>0])

			
			fit <- eval(call("CBPS.fit", X = X, treat = Y, ATT=ATT, 
						           intercept = attr(mt, "intercept") > 0L, method=method, iterations=iterations, 
                       standardize = standardize, 
    				           twostep = twostep))	
				
			fit$model <- mf
			fit$na.action <- attr(mf, "na.action")
			xlevels <- .getXlevels(mt, mf)
			fit$call <- call
			fit$formula <- formula
			fit$terms<-mt
		fit
	}

	CBPS.fit<-function(treat, X, ATT, X.bal=X, method, iterations, standardize, twostep, ...){
		k=0
		if(method=="over") bal.only=FALSE
		if(method=="exact") bal.only=TRUE
		
		probs.min<-1e-6
		names.X<-colnames(X)
		names.X[apply(X,2,sd)==0]<-"(Intercept)"
		
	  #######Declare some constants and orthogonalize Xdf.
	  X.orig<-X
	  format.bal<-F
	  if(sum(dim(X.bal)[2]==dim(X)[2])) format.bal<-T
	  x.sd<-apply(as.matrix(X[,-1]),2,sd)
	  Dx.inv<-diag(c(1,x.sd))
	  diag(Dx.inv)<-1
	  x.mean<-apply(as.matrix(X[,-1]),2,mean)
	  X[,-1]<-apply(as.matrix(X[,-1]),2,FUN=function(x) (x-mean(x))/sd(x))
	  if(k==0) k<-sum(diag(t(X)%*%X%*%ginv(t(X)%*%X)))
	  k<-floor(k+.1)
	  svd1<-svd(X)
	  X<-svd1$u[,1:k]
	  XprimeX.inv<-ginv(t(X)%*%X)

	  #Determine if input presents binary, 3, 4, or continuous treatment case.	
	  no.treats<-length(levels(as.factor(treat)))

	  if (no.treats == 2)
	  {
		output<-CBPS.2Treat(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, ATT, standardize = standardize, twostep = twostep)

	  }
	  
	  if (no.treats == 3)
	  {
		output<-CBPS.3Treat(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, standardize = standardize, twostep = twostep)
	  }
	  
	  if (no.treats == 4)
	  {
		output<-CBPS.4Treat(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, standardize = standardize, twostep = twostep)
	  }

	  if (no.treats > 4)
	  {
		output<-CBPS.Continuous(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, standardize = standardize, twostep = twostep)
		
		d.inv<- svd1$d
		d.inv[d.inv> 1e-5]<-1/d.inv[d.inv> 1e-5]
		d.inv[d.inv<= 1e-5]<-0
		beta.opt<-svd1$v%*%diag(d.inv)%*%coef(output)
		beta.opt[-1,]<-beta.opt[-1,]/x.sd
	  
		beta.opt[1,]<-beta.opt[1,]-matrix(x.mean%*%beta.opt[-1,])
		output$coefficients<-as.matrix(beta.opt)
		rownames(output$coefficients)<-c(names.X)
		output$x<-X.orig
		
		var.1<-output$var
		output$var<-ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)
		rownames(output$var)<-names.X
		colnames(output$var)<-rownames(output$var)
	  }

	  if (no.treats %in% c(2,3,4))
	  {
		d.inv<- svd1$d
		d.inv[d.inv> 1e-5]<-1/d.inv[d.inv> 1e-5]
		d.inv[d.inv<= 1e-5]<-0
		beta.opt<-svd1$v%*%diag(d.inv)%*%coef(output)
		beta.opt[-1,]<-beta.opt[-1,]/x.sd
	  
		beta.opt[1,]<-beta.opt[1,]-matrix(x.mean%*%beta.opt[-1,])
		output$coefficients<-beta.opt
		rownames(output$coefficients)<-names.X
		output$x<-X.orig

		variance<-output$var
		if (no.treats == 2){
			colnames(output$coefficients)<-c("Treated")
			output$var<-ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%variance%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)
			colnames(output$var)<-names.X
			rownames(output$var)<-colnames(output$var)
		}
		
		if (no.treats == 3){
			colnames(output$coefficients)<-levels(as.factor(treat))[c(2,3)]
			var.1.1<-variance[1:k,1:k]
			var.1.2<-variance[1:k,(k+1):(2*k)]
			var.2.1<-variance[(k+1):(2*k),1:k]
			var.2.2<-variance[(k+1):(2*k),(k+1):(2*k)]
			trans.var.1.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.1.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.2.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.2.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			output$var<-rbind(cbind(trans.var.1.1,trans.var.1.2),cbind(trans.var.2.1,trans.var.2.2))
			colnames(output$var)<-c(paste0(levels(as.factor(treat))[2],": ", names.X),paste0(levels(as.factor(treat))[3], ": ", names.X))
			rownames(output$var)<-colnames(output$var)
		}
		
		if (no.treats == 4)
		{
			colnames(output$coefficients)<-levels(as.factor(treat))[c(2,3,4)]
			var.1.1<-variance[1:k,1:k]
			var.1.2<-variance[1:k,(k+1):(2*k)]
			var.1.3<-variance[1:k,(2*k+1):(3*k)]
			var.2.1<-variance[(k+1):(2*k),1:k]
			var.2.2<-variance[(k+1):(2*k),(k+1):(2*k)]
			var.2.3<-variance[(k+1):(2*k),(2*k+1):(3*k)]
			var.3.1<-variance[(2*k+1):(3*k),1:k]
			var.3.2<-variance[(2*k+1):(3*k),(k+1):(2*k)]
			var.3.3<-variance[(2*k+1):(3*k),(2*k+1):(3*k)]
			trans.var.1.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.1.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.1.3<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.3%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.2.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.2.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.2.3<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.3%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.3.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.3.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.3.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.3.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			trans.var.3.3<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.3.3%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
			output$var<-rbind(cbind(trans.var.1.1,trans.var.1.2,trans.var.1.3),cbind(trans.var.2.1,trans.var.2.2,trans.var.2.3),cbind(trans.var.3.1,trans.var.3.2,trans.var.3.3))
			colnames(output$var)<-c(paste0(levels(as.factor(treat))[2],": ", names.X),paste0(levels(as.factor(treat))[3], ": ", names.X),paste0(levels(as.factor(treat))[4], ": ", names.X))
			rownames(output$var)<-colnames(output$var)
		}
	  }
	  
	  output
	}

	print.CBPS <- function(x, digits = max(3, getOption("digits") - 3), ...) {
	  
		cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
			"\n\n", sep = "")
		if (length(coef(x))) {
			cat("Coefficients:\n")
			print.default(format(x$coefficients, digits = digits), 
				print.gap = 2, quote = FALSE)
		}
		else cat("No coefficients\n\n")
    if (max(class(x) == "CBPScontinuous"))
      cat("\nSigma-Squared: ",x$sigmasq)
		cat("\nDegrees of Freedom:", x$df, "\n")
		if (nzchar(mess <- naprint(x$na.action))) 
			cat("  (", mess, ")\n", sep = "")
		cat("Residual Deviance:\t", format(signif(x$deviance, 
			digits)), "\n")
		cat("J-Statistic:\t		", format(signif(x$J)),"\n")
		cat("Log-Likelihood:\t ",-0.5*x$deviance, "\n")
		invisible(x)
	}

	summary.CBPS<-function(object, ...){
	  ##x <- summary.glm(object, dispersion = dispersion, correlation = correlation, symbolic.cor = symbolic.cor, ...)	  
	  x<-NULL
	  names.X<-as.vector(names(object$coefficients))
	  sd.coef <- diag(object$var)^.5
	  coef.table<-(cbind(as.vector(object$coefficients),as.vector(sd.coef),as.vector(object$coefficients/sd.coef),as.vector(2-2*pnorm(abs(object$coefficients/sd.coef)))))
	  colnames(coef.table)<-c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
	  if (ncol(coef(object)) == 1)
	  {
		rownames(coef.table)<-rownames(object$coefficients)#names.X
	  }
	  if (ncol(coef(object)) > 1)
	  {
		  rnames<-array()
		  for (i in 1:ncol(coef(object)))
		  {
			rnames[((i-1)*nrow(coef(object))+1):(i*nrow(coef(object)))]<-paste0(levels(as.factor(object$y))[i],": ",rownames(coef(object)))
		  }
		  rownames(coef.table)<-rnames
	  }
	  
	  pval <- coef.table[,4]
	  symp <- symnum(pval, corr=FALSE,
					 cutpoints = c(0,  .001,.01,.05, .1, 1),
					 symbols = c("***","**","*","."," "))
	  coef.print<-cbind(signif(coef.table,3),as.vector(symp))
	  coef.print[coef.print=="0"]<-"0.000"
		
	  cat("\nCall:	\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
		  "\n", sep = "")
	  
	  cat("\nDeviance Residuals: \n")
	  
	  dev.res<-(object$residuals)/(object$fitted*(1-object$fitted))^.5
	  dev.res[is.infinite(dev.res)]<-0
	  dev.res[is.na(dev.res)]<-0
	  print(summary(dev.res))
		
	  cat("\nCoefficients:\n")

	  print(noquote(coef.print))
	  cat("---\n")
	  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
		#cat("\n	Null J:	 ",object$J)
    if(max(class(object)=="CBPScontinuous")){
      cat("\nSigma-Squared: ",object$sigmasq)
    }
	  cat("\nJ - statistic:	 ",object$J," on ", object$df, " degrees of freedom\n")
	  cat("\nLog-Likelihood: ",-0.5*object$deviance, "\n")
		
	  out<-list("call"=object$call,"coefficients"=coef.table,"J"=object$J)
	  invisible(out)
	}
  

	vcov.CBPS<-function(object,...){
		return(object$var)
	}

	plot.CBPS<-function(x, covars = NULL, silent = TRUE, ...){ 
		bal.x<-balance(x)
    
    if(is.null(covars))
    {
      covars<-1:nrow(bal.x[["balanced"]])
    }
    
		no.treats<-length(levels(as.factor(x$y)))	
		balanced.std.mean<-bal.x[["balanced"]][covars,]
		original.std.mean<-bal.x[["original"]][covars,]
		no.contrasts<-ifelse(no.treats == 2, 1, ifelse(no.treats == 3, 3, 6))
		abs.mean.ori.contrasts<-matrix(rep(0,no.contrasts*nrow(balanced.std.mean)),nrow(balanced.std.mean),no.contrasts)
		abs.mean.bal.contrasts<-matrix(rep(0,no.contrasts*nrow(balanced.std.mean)),nrow(balanced.std.mean),no.contrasts)
		contrast.names<-array()
		true.contrast.names<-array()
		ctr<-1
		for (i in 1:(no.treats-1))
		{
			for (j in (i+1):no.treats)
			{
				abs.mean.ori.contrasts[,ctr]<-abs(original.std.mean[,i+no.treats]-original.std.mean[,j+no.treats])
				abs.mean.bal.contrasts[,ctr]<-abs(balanced.std.mean[,i+no.treats]-balanced.std.mean[,j+no.treats])
				contrast.names[ctr]<-paste0(i,":",j)
				true.contrast.names[ctr]<-paste0(levels(as.factor(x$y))[i],":",levels(as.factor(x$y))[j])
				ctr<-ctr+1
			}
		}
		max.abs.contrast<-max(max(abs.mean.ori.contrasts),max(abs.mean.bal.contrasts))
		m <- matrix(c(1,1,1,2,2,2,3,3,3),nrow = 3,ncol = 3,byrow = TRUE)
		layout(mat = m,heights = c(0.4,0.4,0.3))

		par(mfrow=c(2,1))
		plot(1, type="n", xlim=c(0,max.abs.contrast), ylim=c(1,no.contrasts),xlab="",ylab="",main="",yaxt='n')
		axis(side=2, at=seq(1,no.contrasts),contrast.names)
		mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
		mtext("Contrasts",side=2,line=2)
		mtext("Before Weighting",side=3,line=0.5,font=2)
		for (i in 1:no.contrasts)
		{
			for (j in 1:nrow(balanced.std.mean))
			{
				points(abs.mean.ori.contrasts[j,i],i)
			}
		}
		plot(1, type="n", xlim=c(0,max.abs.contrast), ylim=c(1,no.contrasts),xlab="",ylab="",main="",yaxt='n')
		axis(side=2, at=seq(1,no.contrasts),contrast.names)
		mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
		mtext("Contrasts",side=2,line=2)
		mtext("After Weighting",side=3,line=0.5,font=2)
		for (i in 1:no.contrasts)
		{
			for (j in 1:nrow(balanced.std.mean))
			{
				points(abs.mean.bal.contrasts[j,i],i)
			}
		}

		if(is.null(rownames(balanced.std.mean))) rownames(balanced.std.mean)<-paste0("X",seq(1:nrow(balanced.std.mean)))
    if(!silent) return(list("original"=abs.mean.ori.contrasts,"balanced"=abs.mean.bal.contrasts))
	}
	
	plot.CBPSContinuous<-function(x, covars = NULL, silent = TRUE, ...){ 
		bal.x<-balance(x)
    if (is.null(covars))
    {
      covars<-1:nrow(bal.x[["balanced"]])
    }
		balanced.abs.cor<-abs(bal.x[["balanced"]][covars])
		original.abs.cor<-abs(bal.x[["original"]][covars])

		max.abs.cor<-max(max(original.abs.cor),max(balanced.abs.cor))

		plot(1, type="n", xlim=c(0,max.abs.cor), ylim=c(1,4),xlab="",ylab="",main="",yaxt='n')
		axis(side=2, at=seq(2,3),c("CBPS Weighted", "Unweighted"))
		mtext("Absolute Pearson Correlations",side=1,line=2.25)

		points(x=original.abs.cor, y=rep(3, length(covars)), pch=19)
		points(x=balanced.abs.cor, y=rep(2, length(covars)), pch=19)
    
    if(!silent) return(list("original"=original.abs.cor,"balanced"=balanced.abs.cor))
	}

	balance<-function(object, ...)
	{
		UseMethod("balance")
	}

	balance.CBPS<-function(object, ...){
		treats<-as.factor(object$y)
		treat.names<-levels(treats)
		X<-object$x
		bal<-matrix(rep(0,(ncol(X)-1)*2*length(treat.names)),ncol(X)-1,2*length(treat.names))
		baseline<-matrix(rep(0,(ncol(X)-1)*2*length(treat.names)),ncol(X)-1,2*length(treat.names))
		w<-object$weights
		cnames<-array()
		
		for (i in 1:length(treat.names))
		{
			for (j in 2:ncol(X))
			{
				bal[j-1,i]<-sum((treats==treat.names[i])*X[,j]*w)/sum(w*(treats==treat.names[i]))
				bal[j-1,i+length(treat.names)]<-bal[j-1,i]/sd(X[,j])
				baseline[j-1,i]<-mean(X[which(treats==treat.names[i]),j])
				baseline[j-1,i+length(treat.names)]<-baseline[j-1,i]/sd(X[,j])
			}
			cnames[i]<-paste0(treat.names[i],".mean")
			cnames[length(treat.names)+i]<-paste0(treat.names[i],".std.mean")
		}
		colnames(bal)<-cnames
		rownames(bal)<-colnames(X)[-1]
		colnames(baseline)<-cnames
		rownames(baseline)<-colnames(X)[-1]
		out<-list(balanced=bal,original=baseline)
		out
	}
	
	balance.CBPSContinuous<-function(object, ...){
		treat<-object$y
		X<-object$x
		bal<-matrix(rep(0,(ncol(X)-1)),ncol(X)-1,1)
		baseline<-matrix(rep(0,(ncol(X)-1)),ncol(X)-1,1)
		w<-object$weights
    n<-length(w)
		cnames<-array()
		
		for (j in 2:ncol(X))
		{
			bal[j-1,1]<-(mean(w*X[,j]*treat) - mean(w*X[,j])*mean(w*treat)*n/sum(w))/(sqrt(mean(w*X[,j]^2) - mean(w*X[,j])^2*n/sum(w))*sqrt(mean(w*treat^2) - mean(w*treat)^2*n/sum(w)))
			baseline[j-1,1]<-cor(treat, X[,j], method = "pearson")
		}
		
		colnames(bal)<-"Pearson Correlation"
		rownames(bal)<-colnames(X)[-1]
		colnames(baseline)<-"Pearson Correlation"
		rownames(baseline)<-colnames(X)[-1]
		out<-list(balanced=bal,original=baseline)
		out
	}	
	#########################################
