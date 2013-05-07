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
	treat.all<-unlist(all.1[1,])
	X.all<-NULL
	for(i in 1:dim(all.1)[2])
		X.all<-rbind(X.all,as.matrix(all.1[2,i][[1]]))
	X.all<-cbind(1,X.all[,apply(X.all,2,sd)>0])
		time.all<-rep(1:dim(all.1)[2],each=length(unlist(all.1[1,1])))
	fit <- eval(call("CBMSM.fit", X = X.all, treat = treat.all, time=time.all, method=method))
					}

  fit
}

CBPS.fit<-function(treat, X, ATT, X.bal=X, method, iterations, ...){
	k=0
	if(method=="over") bal.only=FALSE
	if(method=="exact") bal.only=TRUE
	
	probs.min<-1e-6
	names.X<-colnames(X)
	if(is.null(colnames(X))) names.X<-paste0("X",1:ncol(X))
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
	if (no.treats == 2) {
		colnames(output$coefficients)<-c("1")
		names(output$coefficients)<-names.X
		} else {
		colnames(output$coefficients)<-paste0("", 2:(no.treats))
		names(output$coefficients)<-rep(names.X,no.treats-1)

		}

	variance<-output$var
	if (no.treats == 2){
		output$var<-ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%variance%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)
		colnames(output$var)<-paste0("1:", names.X[1:k])
		rownames(output$var)<-colnames(output$var)
	}
	
	if (no.treats == 3){
		var.1.1<-variance[1:k,1:k]
		var.1.2<-variance[1:k,(k+1):(2*k)]
		var.2.1<-variance[(k+1):(2*k),1:k]
		var.2.2<-variance[(k+1):(2*k),(k+1):(2*k)]
		trans.var.1.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
		trans.var.1.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.1.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
		trans.var.2.1<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.1%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
		trans.var.2.2<-Dx.inv%*%ginv(t(X.orig)%*%X.orig)%*%t(X.orig)%*%X%*%svd1$v%*%ginv(diag(svd1$d))%*%var.2.2%*%ginv(diag(svd1$d))%*%t(svd1$v)%*%t(X)%*%X.orig%*%ginv(t(X.orig)%*%X.orig)%*%Dx.inv
		output$var<-rbind(cbind(trans.var.1.1,trans.var.1.2),cbind(trans.var.2.1,trans.var.2.2))
		colnames(output$var)<-c(paste0("2:", names.X[1:k]),paste0("3:", names.X[1:k]))
		rownames(output$var)<-colnames(output$var)
	}
	
	if (no.treats == 4)
	{
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
		colnames(output$var)<-c(paste0("2:", names.X[1:k]),paste0("3:", names.X[1:k]),paste0("4:", names.X[1:k]))
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
	cat("Log-Likelihood:\t ",sum(x$y*log(x$fitted.values) + (1-x$y)*log(1-x$fitted.values)), "\n")
	invisible(x)
}

summary.CBPS<-function(object, ...){
  ##x <- summary.glm(object, dispersion = dispersion, correlation = correlation, symbolic.cor = symbolic.cor, ...)
  
  x<-NULL
  names.X<-as.vector(names(object$coefficients))
  sd.coef <- diag(object$var)^.5
  coef.table<-(cbind(as.vector(object$coefficients),as.vector(sd.coef),as.vector(object$coefficients/sd.coef),as.vector(2-2*pnorm(abs(object$coefficients/sd.coef)))))
  colnames(coef.table)<-c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(coef.table)<-names(object$coefficients)#names.X
  
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
  cat("\nLog-Likelihood: ",sum(object$y*log(object$fitted.values) + (1-object$y)*log(1-object$fitted.values)), "\n")
	
  out<-list("call"=object$call,"coefficients"=coef.table,"J"=object$J)
  invisible(out)
  
}

print.CBMSM<-summary.CBMSM<-function(x){
print("A Covariate Balancing MSM Model")
}

vcov.CBPS<-function(object,...){
	return(object$var)
}

#plot.CBPSlogit<-function(x,...){ 
#	plot(1)
#}

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
