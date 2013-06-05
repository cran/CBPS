CBPS <- function(formula, data, na.action, ATT=NULL, method="over",type="propensity", iterations=NULL, ...) {

	if(type=="propensity") {
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
  
  X <- if (!is.empty.model(mt)) 
	model.matrix(mt, mf)#[,-2]
  else matrix(, NROW(Y), 0L)
	
	X<-cbind(1,X[,apply(X,2,sd)>0])

	 

  fit <- eval(call("CBPS.fit", X = X, treat = Y, ATT=ATT, 
				   intercept = attr(mt, "intercept") > 0L, method=method, iterations=iterations))	
	
  ##if (model) 
  fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  xlevels <- .getXlevels(mt, mf)
  fit$call <- call
  fit$formula <- formula
  fit$terms<-mt
				   

	}
	
	if(type=="MSM")	 {
	all.1<-sapply(formula[[1]],model.frame)
	names.X<-NULL
	for(i in 1:dim(all.1)[2]) {
		names.X<-c(names.X,paste("Intercept.",i,sep=""))
		for(j in 1:dim(all.1[2,i]$X1)[2]){
			 name.temp<-paste("X.",j,".",i,sep="")
			names.X<-c(names.X,name.temp)
		}
		}
	treat.all<-unlist(all.1[1,])
	X.all<-NULL
	for(i in 1:dim(all.1)[2])
		X.all<-rbind(X.all,as.matrix(all.1[2,i][[1]]))
	X.all<-cbind(1,X.all[,apply(X.all,2,sd)>0])
		time.all<-rep(1:dim(all.1)[2],each=length(unlist(all.1[1,1])))
	fit <- eval(call("CBMSM.fit", X = X.all, treat = treat.all, time=time.all, method=method))
	class(fit$coefficients)<-"coefficients"
	names(fit$coefficients)<-names.X
	fit$call<-formula
					}

	if(type=="MultiBin")	 {
      	all.1<-sapply(formula[[1]],model.frame)
	names.X<-NULL
	for(i in 1:dim(all.1)[2]) {
		names.X<-c(names.X,paste("Intercept.",i,sep=""))
		for(j in 1:dim(all.1[2,i]$X1)[2]){
			 name.temp<-paste("X.",j,".",i,sep="")
			names.X<-c(names.X,name.temp)
		}
		}
	treat.all<-unlist(all.1[1,])
	X.all<-NULL
	for(i in 1:dim(all.1)[2])
		X.all<-rbind(X.all,as.matrix(all.1[2,i][[1]]))
	X.all<-cbind(1,X.all[,apply(X.all,2,sd)>0])
		time.all<-rep(1:dim(all.1)[2],each=length(unlist(all.1[1,1])))
        fit <- eval(call("CBMSM.fit", X = X.all, treat = treat.all, time=time.all, method=method, MultiBin.fit = TRUE))
	class(fit$coefficients)<-"coefficients"
	names(fit$coefficients)<-names.X
	fit$call<-formula
    }    
    fit
}

CBPS.fit<-function(treat, X, ATT, X.bal=X, method, iterations, ...){
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
  x.sd<-apply(X[,-1],2,sd)
  Dx.inv<-diag(c(1,x.sd))
  diag(Dx.inv)<-1
  x.mean<-apply(X[,-1],2,mean)
  X[,-1]<-apply(X[,-1],2,FUN=function(x) (x-mean(x))/sd(x))
  if(k==0) k<-sum(diag(t(X)%*%X%*%ginv(t(X)%*%X)))
  k<-floor(k+.1)
  svd1<-svd(X)
  X<-svd1$u[,1:k]
  XprimeX.inv<-ginv(t(X)%*%X)


  #Determine if input presents binary, 3, 4, or continuous treatment case.	
  no.treats<-length(levels(as.factor(treat)))

  if (no.treats == 2)
  {
 	output<-CBPS.2Treat(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations,ATT)
  }
  
  if (no.treats == 3)
  {
	output<-CBPS.3Treat(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations)
  }
  
  if (no.treats == 4)
  {
	output<-CBPS.4Treat(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations)
  }

  if (no.treats > 4)
  {
	stop("Use only treatments which can take 2, 3, or 4 values")
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
	#output$x<-X

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
  cat("\nJ - statistic:	 ",object$J," on ", object$df, " degrees of freedom\n")
  cat("\nLog-Likelihood: ",-0.5*object$deviance, "\n")
	
  out<-list("call"=object$call,"coefficients"=coef.table,"J"=object$J)
  invisible(out)
  
}

summary.CBMSM<-summary.CBMB<-function(object,...){
 x<-NULL
  names.X<-as.vector(names(object$coefficients))
  sd.coef <- diag(object$var)^.5
  coef.table<-(cbind(as.vector(object$coefficients),as.vector(sd.coef),as.vector(object$coefficients/sd.coef),as.vector(2-2*pnorm(abs(object$coefficients/sd.coef)))))
  colnames(coef.table)<-c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if (ncol(coef(object)) == 1)
  {
	rownames(coef.table)<-names.X
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
  cat("\nJ - statistic:	 ",object$J," on ", object$df, " degrees of freedom\n")
  cat("\nLog-Likelihood: ",-0.5*object$deviance, "\n")
	
  out<-list("call"=object$call,"coefficients"=coef.table,"J"=object$J)
  invisible(out)
}

print.CBMSM<-print.CBMB<-function(x,...){
  names.X<-as.vector(names(x$coefficients))
  sd.coef <- diag(x$var)^.5
  coef.table<-(cbind(as.vector(x$coefficients),as.vector(sd.coef),as.vector(x$coefficients/sd.coef),as.vector(2-2*pnorm(abs(x$coefficients/sd.coef)))))
  colnames(coef.table)<-c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  if (ncol(coef(x)) == 1)
  {
	rownames(coef.table)<-names.X
  }
  if (ncol(coef(x)) > 1)
  {
	  rnames<-array()
	  for (i in 1:ncol(coef(x)))
	  {
		rnames[((i-1)*nrow(coef(x))+1):(i*nrow(coef(x)))]<-paste0(levels(as.factor(x$y))[i],": ",rownames(coef(x)))
	  }
	  rownames(coef.table)<-rnames
  }
  
  pval <- coef.table[,4]
  symp <- symnum(pval, corr=FALSE,
				 cutpoints = c(0,  .001,.01,.05, .1, 1),
				 symbols = c("***","**","*","."," "))
  coef.print<-cbind(signif(coef.table,3),as.vector(symp))
  coef.print[coef.print=="0"]<-"0.000"
	
  cat("\nCall:	\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
	  "\n", sep = "")
  
  cat("\nDeviance Residuals: \n")
  
  dev.res<-(x$residuals)/(x$fitted*(1-x$fitted))^.5
  dev.res[is.infinite(dev.res)]<-0
  dev.res[is.na(dev.res)]<-0
  print(summary(dev.res))
	
  cat("\nCoefficients:\n")

  print(noquote(coef.print))
  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1 \n")
	#cat("\n	Null J:	 ",x$J)
  cat("\nJ - statistic:	 ",x$J," on ", x$df, " degrees of freedom\n")
  cat("\nLog-Likelihood: ",-0.5*x$deviance, "\n")
	
  out<-list("call"=x$call,"coefficients"=coef.table,"J"=x$J)
  invisible(out)

}


vcov.CBPS<-function(object,...){
	return(object$var)
}

plot.CBPS<-function(x,...){ 
	bal.x<-balance(x)
	no.treats<-length(levels(as.factor(x$y)))	
	balanced.std.mean<-bal.x[["balanced"]]
	original.std.mean<-bal.x[["original"]]
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
			abs.mean.ori.contrasts[,ctr]<-abs(original.std.mean[,2*i]-original.std.mean[,2*j])
			abs.mean.bal.contrasts[,ctr]<-abs(balanced.std.mean[,2*i]-balanced.std.mean[,2*j])
			contrast.names[ctr]<-paste0(i,":",j)
			true.contrast.names[ctr]<-paste0(levels(as.factor(x$y))[i],":",levels(as.factor(x$y))[j])
			ctr<-ctr+1
		}
	}
	max.abs.contrast<-max(max(abs.mean.ori.contrasts),max(abs.mean.bal.contrasts))
	#par(mfrow=c(2,1))
	m <- matrix(c(1,1,1,2,2,2,3,3,3),nrow = 3,ncol = 3,byrow = TRUE)
	layout(mat = m,heights = c(0.4,0.4,0.3))

	plot(1, type="n", xlim=c(0,max.abs.contrast), ylim=c(1,no.contrasts),xlab="",ylab="",main="",yaxt='n')
	axis(side=2, at=seq(1,no.contrasts),contrast.names)
	mtext("Absolute Difference of Standardized Means",side=1,line=2.25)
	mtext("Contrasts",side=2,line=2)
	mtext("Before Weighting",side=3,line=0.5,font=2)
	for (i in 1:no.contrasts)
	{
		for (j in 1:nrow(balanced.std.mean))
		{
			par(pch=j)
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
			par(pch=j)
			points(abs.mean.bal.contrasts[j,i],i)
		}
	}

	if(is.null(rownames(balanced.std.mean))) rownames(balanced.std.mean)<-paste0("X",seq(1:nrow(balanced.std.mean)))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="")
	legend(x = "top",inset = 0,
        legend = rownames(balanced.std.mean), 
        pch=seq(from=1,to=nrow(balanced.std.mean)), horiz = TRUE)
	print(data.frame("Coded Contrasts"=contrast.names,"Named Contrasts"=true.contrast.names))
}

plot.CBMB<-function(x, stabilized = TRUE, ...){ 
	bal.x<-balance.CBMB(x, stabilized)
	if (!missing(stabilized)) bal.x<-balance.CBMB(x, stabilized)
	treats<-matrix(x$y,ncol=length(unique(x$time)),byrow=FALSE)
	decimal.treats<-rep(0,dim(treats)[1])
	for(i in 1:dim(treats)[1])
	{
		for(j in 1:dim(treats)[2])
		{
			decimal.treats[i]<-decimal.treats[i]+2^(j-1)*treats[i,j]
		}
	}
	no.treats<-2^dim(treats)[2]
	treat.names<-array()
	for(i in 0:(no.treats-1))
	{
		binary.treat.name<-treats[which(decimal.treats==i)[1],]
		treat.names[i+1]<-paste(binary.treat.name[1:dim(treats)[2]],collapse="")
	}
	balanced.std.mean<-bal.x[["balanced"]]
	original.std.mean<-bal.x[["original"]]
	no.contrasts<-no.treats-1
	abs.mean.ori.contrasts<-matrix(rep(0,no.contrasts*nrow(balanced.std.mean)),nrow(balanced.std.mean),no.contrasts)
	abs.mean.bal.contrasts<-matrix(rep(0,no.contrasts*nrow(balanced.std.mean)),nrow(balanced.std.mean),no.contrasts)
	contrast.names<-array()
	true.contrast.names<-array()
	for (i in 2:no.treats)
	{
		abs.mean.ori.contrasts[,i-1]<-abs(original.std.mean[,2]-original.std.mean[,2*i])
		abs.mean.bal.contrasts[,i-1]<-abs(balanced.std.mean[,2]-balanced.std.mean[,2*i])
		contrast.names[i-1]<-paste0(1,":",i)
		true.contrast.names[i-1]<-paste0(treat.names[1],":",treat.names[i])
	}
	max.abs.contrast<-max(max(abs.mean.ori.contrasts),max(abs.mean.bal.contrasts))
	
	#par(mfrow=c(2,1))
	m <- matrix(c(1,1,1,2,2,2,3,3,3),nrow = 3,ncol = 3,byrow = TRUE)
	layout(mat = m,heights = c(0.4,0.4,0.3))
	plot(1, type="n", xlim=c(0,max.abs.contrast), ylim=c(1,no.contrasts),xlab="",ylab="",main="",yaxt='n')
	axis(side=2, at=seq(1,no.contrasts),contrast.names)
	mtext("Absolute Difference of Standardized Means",side=1,line=2.5)
	mtext("Contrasts",side=2,line=2)
	mtext("Before Weighting",side=3,line=0.5,font=2)
	for (i in 1:no.contrasts)
	{
		for (j in 1:nrow(balanced.std.mean))
		{
			par(pch=j)
			points(abs.mean.ori.contrasts[j,i],i)
		}
	}
	plot(1, type="n", xlim=c(0,max.abs.contrast), ylim=c(1,no.contrasts),xlab="",ylab="",main="",yaxt='n')
	axis(side=2, at=seq(1,no.contrasts),contrast.names)
	mtext("Absolute Difference of Standardized Means",side=1,line=2.5)
	mtext("Contrasts",side=2,line=2)
	mtext("After Weighting",side=3,line=0.5,font=2)
	for (i in 1:no.contrasts)
	{
		for (j in 1:nrow(balanced.std.mean))
		{
			par(pch=j)
			points(abs.mean.bal.contrasts[j,i],i)
		}
	}
	if(is.null(rownames(balanced.std.mean))) rownames(balanced.std.mean)<-paste0("X",seq(1:nrow(balanced.std.mean)))
	plot(1, type = "n", axes=FALSE, xlab="", ylab="")
	legend(x = "top",inset = 0,
		legend = rownames(balanced.std.mean), 
		pch=seq(from=1,to=nrow(balanced.std.mean)), horiz = TRUE)
	print(data.frame("Coded Contrasts"=contrast.names,"Named Contrasts"=true.contrast.names))
}

plot.CBMSM<-function(x, stabilized = TRUE, ...){
	bal.x<-balance.CBMSM(x, stabilized)
	if (!missing(stabilized)) bal.x<-balance.CBMSM(x, stabilized)
	treats<-matrix(x$y,ncol=length(unique(x$time)),byrow=FALSE)
	decimal.treats<-rep(0,dim(treats)[1])
	for(i in 1:dim(treats)[1])
	{
		for(j in 1:dim(treats)[2])
		{
			decimal.treats[i]<-decimal.treats[i]+2^(j-1)*treats[i,j]
		}
	}
	no.treats<-2^dim(treats)[2]
	treat.names<-array()
	for(i in 0:(no.treats-1))
	{
		binary.treat.name<-treats[which(decimal.treats==i)[1],]
		treat.names[i+1]<-paste(binary.treat.name[1:dim(treats)[2]],collapse="")
	}
	balanced.std.mean<-array()
	original.std.mean<-array()
	balanced.abs.diff<-array()
	original.abs.diff<-array()
	no.contrasts<-no.treats-1
	original.abs.diff<-matrix(rep(0,no.contrasts*(dim(x$x)[[2]]-1)*length(bal.x)),length(bal.x)*(dim(x$x)[2]-1),no.contrasts)
	balanced.abs.diff<-matrix(rep(0,no.contrasts*(dim(x$x)[[2]]-1)*length(bal.x)),length(bal.x)*(dim(x$x)[2]-1),no.contrasts)
	
	for (i in 1:length(bal.x))
	{
		this.bal.x<-bal.x[[i]]
		this.original<-this.bal.x[["original"]]
		this.balanced<-this.bal.x[["balanced"]]
		for (j in 1:nrow(this.original))
		{
			for (k in seq(from=2,to=(dim(this.original)[2]),by=2))
			{
				original.abs.diff[(i-1)*nrow(this.original)+j,(k-2)/2]<-abs(this.original[j,k]-this.original[j,2])
				balanced.abs.diff[(i-1)*nrow(this.original)+j,(k-2)/2]<-abs(this.balanced[j,k]-this.balanced[j,2])
			}
		}
	}
	max.abs.diff<-max(max(original.abs.diff,na.rm=TRUE),max(balanced.abs.diff,na.rm=TRUE))
	plot(0,type="n",xlim=c(0,max.abs.diff),ylim=c(0,max.abs.diff),xlab="",ylab="",main="",asp=1)
	abline(0,1)
	points(original.abs.diff,balanced.abs.diff)
	mtext("Before Weighting",side=1,line=2)
	mtext("After Weighting",side=2,line=2)
	mtext("Absolute Differences of Standardized Means \n Before and After Weighting",side=3,line=0.5,font=2)
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
			bal[j-1,2*i-1]<-sum((treats==treat.names[i])*X[,j]*w)/sum(w*(treats==treat.names[i]))
			bal[j-1,2*i]<-bal[j-1,2*i-1]/sd(X[,j])
			baseline[j-1,2*i-1]<-mean(X[which(treats==treat.names[i]),j])
			baseline[j-1,2*i]<-baseline[j-1,2*i-1]/sd(X[,j])
		}
		cnames[2*i-1]<-paste0(treat.names[i],".mean")
		cnames[2*i]<-paste0(treat.names[i],".std.mean")
	}
	colnames(bal)<-cnames
	rownames(bal)<-colnames(X)[-1]
	colnames(baseline)<-cnames
	rownames(baseline)<-colnames(X)[-1]
	out<-list(balanced=bal,original=baseline)
	out
}

balance.CBMB<-function(object, stabilized = TRUE, ...)
{
	treats<-matrix(object$y,ncol=length(unique(object$time)),byrow=FALSE)
	n<-nrow(treats)
	decimal.treats<-rep(0,dim(treats)[1])
	for(i in 1:dim(treats)[1])
	{
		for(j in 1:dim(treats)[2])
		{
			decimal.treats[i]<-decimal.treats[i]+2^(j-1)*treats[i,j]
		}
	}
	no.treats<-2^dim(treats)[2]
	treat.names<-array()
	for(i in 0:(no.treats-1))
	{
		binary.treat.name<-treats[which(decimal.treats==i)[1],]
		treat.names[i+1]<-paste(binary.treat.name[1:dim(treats)[2]],collapse="")
	}
	X<-object$x[1:n,]
	bal<-matrix(rep(0,(ncol(X)-1)*2*no.treats),ncol(X)-1,2*no.treats)
	baseline<-matrix(rep(0,(ncol(X)-1)*2*no.treats),ncol(X)-1,2*no.treats)
	w<-object$weights[["stabilized"]]
	if(!stabilized) w<-object$weights[["unstabilized"]]
	cnames<-array()
	for (i in 1:no.treats)
	{
		for (j in 2:ncol(X))
		{
			bal[j-1,2*i-1]<-sum((decimal.treats==(i-1))*X[,j]*w)/sum(w*(decimal.treats==(i-1)))
			bal[j-1,2*i]<-bal[j-1,2*i-1]/sd(X[,j])
			baseline[j-1,2*i-1]<-mean(X[which(decimal.treats==(i-1)),j])
			baseline[j-1,2*i]<-baseline[j-1,2*i-1]/sd(X[,j])
		}
		cnames[2*i-1]<-paste0(treat.names[i],".mean")
		cnames[2*i]<-paste0(treat.names[i],".std.mean")
	}
	colnames(bal)<-cnames
	rownames(bal)<-colnames(X)[-1]
	colnames(baseline)<-cnames
	rownames(baseline)<-colnames(X)[-1]
	out<-list(balanced=bal,original=baseline)
	out
}

balance.CBMSM<-function(object, stabilized = TRUE, ...)
{
	if (missing(stabilized)) stabilized<-TRUE
	treats<-matrix(object$y,ncol=length(unique(object$time)),byrow=FALSE)
	n<-nrow(treats)
	decimal.treats<-rep(0,dim(treats)[1])
	for(i in 1:dim(treats)[1])
	{
		for(j in 1:dim(treats)[2])
		{
			decimal.treats[i]<-decimal.treats[i]+2^(j-1)*treats[i,j]
		}
	}
	no.treats<-2^dim(treats)[2]
	treat.names<-array()
	for(i in 0:(no.treats-1))
	{
		binary.treat.name<-treats[which(decimal.treats==i)[1],]
		treat.names[i+1]<-paste(binary.treat.name[1:dim(treats)[2]],collapse="")
	}
	X<-list()
	for(i in 1:dim(treats)[2])
	{
		X[[i]]<-matrix(object$x[(1:n)+(i-1)*n,],n,dim(object$x)[2])
	}
	w<-object$weights[["stabilized"]]
	if(!stabilized) w<-object$weights[["unstabilized"]]
	cnames<-array()
	out<-list()
	for (k in 1:length(X))
	{
		bal<-matrix(rep(0,(ncol(X[[1]])-1)*2*no.treats),ncol(X[[1]])-1,2*no.treats)
		baseline<-matrix(rep(0,(ncol(X[[1]])-1)*2*no.treats),ncol(X[[1]])-1,2*no.treats)
		this.X<-X[[k]]
		for (i in 1:no.treats)
		{
			for (j in 2:ncol(this.X))
			{
				bal[j-1,2*i-1]<-sum((decimal.treats==(i-1))*this.X[,j]*w)/sum(w*(decimal.treats==(i-1)))
				bal[j-1,2*i]<-bal[j-1,2*i-1]/sd(this.X[,j])
				baseline[j-1,2*i-1]<-mean(this.X[which(decimal.treats==(i-1)),j])
				baseline[j-1,2*i]<-baseline[j-1,2*i-1]/sd(this.X[,j])
			}
			cnames[2*i-1]<-paste0(treat.names[i],".mean")
			cnames[2*i]<-paste0(treat.names[i],".std.mean")
		}
		colnames(bal)<-cnames
		rownames(bal)<-colnames(this.X)[-1]
		colnames(baseline)<-cnames
		rownames(baseline)<-colnames(this.X)[-1]
		out[[k]]<-list(balanced=bal,original=baseline)
	}
	out
}

TE.est<-function(dv, object, M=1){
	
  match.func<-function(probs,log.odds=F,treat,M){
	out1<-sapply(which(treat==1),FUN=function(x) {
	  dist<-abs(probs[treat==0]-probs[x])
	  if(log.odds==T) dist<-abs(log(probs[treat==0]/(1-probs[treat==0])) - log(probs[x]/(1-probs[x])))
				 which(dist<= sort(unique(dist),decreasing=F)[1:M])
	})
	c(which(treat==1),which(treat==0)[unlist(out1)])
  }
  
  X <- data.frame(object$x)
  treat <- object$y
  pscore <- fitted(object)
  matches<-match.func(object$fitted,treat=object$y, M=M)
  output.match <- mean(dv[treat==1])-mean(dv[matches][treat[matches]==0])
  output.ipw <- mean(dv[treat==1])-sum((dv*object$weights)[treat==0])/sum(object$weights[treat==0])
  output.ht <- mean(dv*object$weights)
  output<-list("match"=output.match, "ipw"=output.ipw, "ht"=output.ht)
  return(output)
}

IPW<-function(outcome, treat, data=parent.frame(), pscore, k){
	IPW.inner<-function(outcome, treat, pscore, k){
	n <- length(treat)
	
	ipw.1 <- 1/n*sum(treat*outcome/pscore) - 1/n*sum((1-treat)*outcome/(1-pscore))
	
	mu.2.1 <- sum(treat/pscore)^-1*sum(treat*outcome/pscore)
	mu.2.0 <- sum((1-treat)/(1-pscore))^-1*sum((1-treat)*outcome/(1-pscore))
	ipw.2 <- mu.2.1 - mu.2.0
	
	c1 <- sum((treat-pscore)/pscore)/sum(((treat-pscore)/pscore)^2)
	c0 <- -sum((treat-pscore)/(1-pscore))/sum(((treat-pscore)/(1-pscore))^2)
	mu.3.1 <- sum(treat/pscore*(1-c1/pscore))^-1*sum(treat*outcome/pscore*(1-c1/pscore))
	mu.3.0 <- sum((1-treat)/(1-pscore)*(1-c0/(1-pscore)))^-1*sum((1-treat)*outcome/(1-pscore)*(1-c0/(1-pscore)))
	ipw.3 <- mu.3.1 - mu.3.0

	nu.1 <- -sum(treat*(outcome - mu.3.1)/pscore^2)/sum(((treat-pscore)/pscore)^2)
	nu.0 <- -sum((1-treat)*(outcome - mu.3.0)/(1-pscore)^2)/sum(((treat-pscore)/(1-pscore))^2)
			 
	W <- rep(1, k)
	H1 <- 1/n*sum(treat*outcome*(1-pscore)/pscore + (1-treat)*outcome*pscore/(1-pscore)) * W
	H2 <- 1/n*sum(treat*(outcome - mu.2.1)*(1-pscore)/pscore + (1 - treat)*(outcome-mu.2.0*pscore)/(1-pscore)) * W
	H3 <- 1/n*sum(treat*(outcome-mu.3.1+nu.1)*(1-pscore)/pscore + (1-treat)*(outcome - mu.3.0 + nu.0)*pscore/(1-pscore)) * W
	
	E.inv <- 1/n*sum(pscore*(1-pscore)) * W %*% t(W)
	
	var.ipw.1 <- 1/n^2*sum((treat*outcome/pscore - (1 - treat)*outcome/(1-pscore) - ipw.1 - (treat - pscore) * (t(H1) %*% E.inv %*% W))^2)
	var.ipw.2 <- 1/n^2*sum((treat*(outcome-mu.2.1)/pscore - (1-treat)*(outcome-mu.2.0)/(1-pscore) - (treat-pscore) * (t(H2) %*% E.inv %*% W))^2)
	var.ipw.3 <- 1/n^2*sum(((treat*(outcome-mu.3.1)+nu.1*(treat-pscore))/pscore - ((1-treat)*(outcome-mu.3.0)-nu.0*(treat-pscore))/(1-pscore) - (treat-pscore) * (t(H3) %*% E.inv %*% W))^2)
	se.ipw.1 <- sqrt(var.ipw.1)
	se.ipw.2 <- sqrt(var.ipw.2)
	se.ipw.3 <- sqrt(var.ipw.3)
	
	out <- data.frame(Point.Est=c(ipw.1, ipw.2, ipw.3), Std.Err=c(se.ipw.1, se.ipw.2, se.ipw.3), t.Statistic=c(ipw.1/se.ipw.1, ipw.2/se.ipw.2, ipw.3/se.ipw.3))
	rownames(out) <- c("IPW1", "IPW2", "IPW3")
	
	return(out)
	}
	return(eval(substitute(IPW.inner(outcome, treat, pscore, k)), data, enclos=parent.frame()))
}


DR<-function(formula, model, data, treat, pscore){
	call <- match.call()
	if (missing(data)) 
		data <- environment(formula)
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

	
	if (model == "lm"){
		m1 <- lm(formula=formula, data=data, subset=which(treat==1))
		m0 <- lm(formula=formula, data=data, subset=which(treat==0))
	}
	else if (model == "glm"){
		m1 <- glm(formula=formula, data=data, subset=which(treat==1))
		m0 <- glm(formula=formula, data=data, subset=which(treat==1))
	}
	
	DR.inner<-function(Y,m1,m0,data,treat,pscore){
		n <- length(treat)
		dr <- 1/n*sum((treat*Y - (treat - pscore)*predict(m1,data))/pscore) - 1/n*sum(((1-treat)*Y+(treat-pscore)*predict(m0,data))/(1-pscore))
		var.dr <- 1/n^2*sum(((treat*Y - predict(m1,data)*(treat-pscore))/pscore - ((1-treat)*Y + predict(m0,data)*(treat-pscore))/(1-pscore) - dr)^2)
		out <- c(dr, sqrt(var.dr), dr/sqrt(var.dr))
		names(out) <- c("Point Est", "Std. Err", "t-Statistic")
		return(out)
	}
	return(eval(substitute(DR.inner(Y,m1,m0,data,treat,pscore)), data, enclos=parent.frame()))
}


	
#########################################
