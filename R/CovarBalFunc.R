CBPS <- function(formula, data, na.action, ATT, method="over",...) {


	bayes<-FALSE
	cauchy<-FALSE
 
  call <- match.call()
  family <- binomial()
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
  
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf)#[,-2]
  else matrix(, NROW(Y), 0L)
    
	X<-cbind(1,X[,-1][,apply(X[,-1],2,sd)>0])
	
  fit <- eval(call("CBPS.fit", X = X, treat = Y, ATT=ATT, family = family, bayes=bayes, cauchy=cauchy,
                   intercept = attr(mt, "intercept") > 0L, method=method))
  ##if (model) 
  fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  xlevels <- .getXlevels(mt, mf)
  fit$call <- call
  fit$formula <- formula
  fit$terms<-mt
  fit
}

CBPS.fit<-function(treat, X, ATT, X.bal=X, method, ...){
	k=0
	if(method=="over") bal.only=FALSE
	if(method=="exact") bal.only=TRUE
		
##A tolerance, to set probabilities to zero.
	bayes<-cauchy<-FALSE
  probs.min<-1e-10
  names.X<-colnames(X)

  
  ##Generates ATT weights.  Called by loss function, etc.
	ATT.wt.func<-function(beta.curr,X.wt=X){
    X<-as.matrix(X.wt)
    n<-dim(X)[1]
    n.c<-sum(treat==0)
    n.t<-sum(treat==1)
    theta.curr<-as.vector(X%*%beta.curr)
	if (cauchy) {probs.curr<-.5+1/pi*atan(theta.curr)}
	else {probs.curr<-(1+exp(-theta.curr))^-1}
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)	
    w1<-(n/n.t*(treat-probs.curr)/(1-probs.curr))
    w1[treat==1]<-n/n.t
    w1
  }
  
  ##The gmm objective function--given a guess of beta, constructs the GMM J statistic.
  gmm.func<-function(beta.curr,X.gmm=X,ATT.gmm=ATT){
    ##Designate a few objects in the function.
    X<-as.matrix(X.gmm)
    ATT<-ATT.gmm
    
    ##Designate sample size, number of treated and control observations,
    ##theta.curr, which are used to generate probabilities.
    ##Trim probabilities, and generate weights.
    n<-dim(X)[1]
    n.c<-sum(treat==0)
    n.t<-sum(treat==1)
    theta.curr<-as.vector(X%*%beta.curr)
	if (cauchy){probs.curr<-.5+1/pi*atan(theta.curr)}
    else{probs.curr<-(1+exp(-theta.curr))^-1}
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)	
    probs.curr<-as.vector(probs.curr)
    if(ATT)
      w.curr<-ATT.wt.func(beta.curr)
    else
      w.curr<-(probs.curr-1+treat)^-1
	  
  
    ##Generate the vector of mean imbalance by weights.
    w.curr.del<-1/(n)*t(X)%*%(w.curr)
    w.curr.del<-as.vector(w.curr.del)
    w.curr<-as.vector(w.curr)

    ##Generate g-bar, as in the paper.
    gbar<-c( 1/n*t(X)%*%(treat-probs.curr)*(1-bal.only),w.curr.del)
   if (cauchy==T){
   	s<-(1+theta.curr^2)^-1*pi^-1*(treat/probs.curr-(1-treat)/(1-probs.curr))
   	gbar<-c(1/n*t(X)%*%s,w.curr.del)
   } 
    ##Generate the covariance matrix used in the GMM estimate.
    ##Was for the initial version that calculates the analytic variances.

#X.h<-X*((1-probs.curr)*probs.curr)^.5
	  hat.diag<-0#diag(X.h%*%ginv(t(X.h)%*%X.h)%*%t(X.h))
	if(ATT){
      X.1<-X*((1-probs.curr)*probs.curr)^.5
      X.2<-X*(probs.curr/(1-probs.curr))^.5
      X.1.1<-X*(probs.curr)^.5

    }else{
      X.1<-X*((1-probs.curr)*probs.curr+hat.diag)^.5
      X.2<-X*(probs.curr*(1-probs.curr)+hat.diag)^-.5		
      X.1.1<-  X
    }
    

	  V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.1),
			   1/n*cbind(t(X.1.1)%*%X.1.1,t(X.2)%*%X.2)
			   )		

#X.V<-cbind(X*(treat-probs.curr), X* w.curr)
#	  V<-cov(X.V)

    
    ##Calculate the GMM loss.
    loss1<-as.vector( t(gbar)%*%ginv(V )%*%(gbar))      
    out1<-list("loss"=max(loss1*n,loss1*n), "V"=V)
    out1
  }
  gmm.loss<-function(x,...) gmm.func(x,...)$loss
	
  ##Loss function for balance constraints, returns the squared imbalance along each dimension.
  bal.loss<-function(beta.curr){
    ##Generate theta and probabilities.
    theta.curr<-as.vector(X%*%beta.curr)
	if (cauchy){probs.curr<-.5+1/pi*atan(theta.curr)}
    else{probs.curr<-(1+exp(-theta.curr))^-1}
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)
    ##Generate weights.
    if(ATT)
      w.curr<-ATT.wt.func(beta.curr)
    else
      w.curr<-(probs.curr-1+treat)^-1
    X.2<-X
    ##Generate mean imbalance.
    loss1<-abs(t(w.curr)%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr))
    loss1
  }
	
#######Declare some constants and orthogonalize Xdf.
  format.bal<-F
  if(sum(dim(X.bal)[2]==dim(X)[2])) format.bal<-T
  x.sd<-apply(X[,-1],2,sd)
  x.mean<-apply(X[,-1],2,mean)
  X[,-1]<-apply(X[,-1],2,FUN=function(x) (x-mean(x))/sd(x))
  if(k==0) k<-sum(diag(t(X)%*%X%*%ginv(t(X)%*%X)))
  k<-floor(k+.1)
  svd1<-svd(X)
  X<-svd1$u[,1:k]
  XprimeX.inv<-ginv(t(X)%*%X)

  ##XprimeX.inv<-X%*%ginv(t(X)%*%X)%*%t(X)
  format.bal<-F
  if(format.bal){
    X.bal[,-1]<-apply(X.bal[,-1],2,FUN=function(x) (x-mean(x))/sd(x))
    XprimeX.inv.bal<-X.bal%*%ginv(t(X.bal)%*%X.bal)%*%t(X.bal)
    if(k==0) k<-sum(diag(XprimeX.inv.bal))
    k<-floor(k+.1)
    X.bal<-svd(X.bal)$u[,1:k]
    ##XprimeX.inv.bal<-X.bal%*%ginv(t(X.bal)%*%X.bal)%*%t(X.bal)
  }

  if(!format.bal) {
    X.bal<-X
    ##XprimeX.inv.bal<-XprimeX.inv
  }
  n<-length(treat)
  n.c<-sum(treat==0)
  n.t<-sum(treat==1)
  x.orig<-x<-cbind(as.matrix(X))
  
  ##GLM estimation
  if(cauchy){glm1<-glm(treat~X-1,family=binomial(link="cauchit"))}  else
  	{glm1<-glm(treat~X-1,family=binomial)}
  boundary<-glm1$boundary	
  glm1$coef[is.na(glm1$coef)]<-0
  probs.glm<-glm1$fit
  glm1$fit<-probs.glm<-pmin(1-probs.min,probs.glm)
  glm1$fit<-probs.glm<-pmax(probs.min,probs.glm)	
  beta.curr<-glm1$coef
  beta.curr[is.na(beta.curr)]<-0
  if(bayes) beta.curr<-(bayesglm(treat~X-1,family=binomial))$coef
	
  alpha.func<-function(alpha) gmm.loss(beta.curr*alpha)
  beta.curr<-beta.curr*optimize(alpha.func,interval=c(.8,1.1))$min
  
	
  ##Generate estimates for balance and CBPSE
  gmm.init<-beta.curr
 
  opt.bal<-optim(gmm.init, bal.loss, control=list("maxit"=1000), method="BFGS", hessian=TRUE)
  beta.bal<-opt.bal$par
  
	gmm.glm.init<-optim(glm1$coef,gmm.loss, control=list("maxit"=1000), method="BFGS", hessian=TRUE)
	gmm.bal.init<-optim(beta.bal,gmm.loss, control=list("maxit"=1000), method="BFGS", hessian=TRUE)
	
	if(gmm.glm.init$val<gmm.bal.init$val) opt1<-gmm.glm.init else opt1<-gmm.bal.init
  #print(opt.bal)
 
#if(gmm.loss(beta.curr)>gmm.loss(beta.curr*0)) gmm.init<-beta.curr*0
# if(gmm.loss(gmm.init)>gmm.loss(beta.bal)) gmm.init<-beta.bal

	
	if(bal.only) opt1<-opt.bal #else 
#opt1<-optim(gmm.init, gmm.loss, control=list("maxit"=1000), method="BFGS", hessian=TRUE)
#opt1<-optim(gmm.init, gmm.loss, control=list("maxit"=1000), method="BFGS", hessian=TRUE)

#print(gmm.loss(glm1$coef)) 
#	print(gmm.loss(opt.bal$par)) 
#	print(opt1) 



  ##Generate probabilities
#if(gmm.func(opt1$par)$loss-gmm.func(opt.bal$par)$loss> qchisq(.95,df=dim(X)[2])) beta.opt<-opt1$par else beta.opt<-opt.bal$par
	beta.opt<-opt1$par
	theta.opt<-as.vector(X%*%beta.opt)
  if (cauchy){probs.opt<-.5+1/pi*atan(theta.opt)}  else{probs.opt<-(1+exp(-theta.opt))^-1}
  probs.opt<-pmin(1-probs.min,probs.opt)
  probs.opt<-pmax(probs.min,probs.opt)
	
  ##Generate weights
  if(ATT){
    w.opt<-ATT.wt.func(beta.opt) 
  }else{
    w.opt<-(probs.opt-1+treat)^-1
  }
  
  J.opt<-gmm.func(beta.opt)$loss
  
  residuals<-treat-probs.opt
  deviance <- -2*c(sum(treat*log(probs.glm)+(1-treat)*log(1-probs.glm)))
  nulldeviance <- -2*c(sum(treat*log(mean(treat))+(1-treat)*log(1-mean(treat))))		
  
  d.inv<- svd1$d
  d.inv[d.inv> 1e-5]<-1/d.inv[d.inv> 1e-5]
  d.inv[d.inv<= 1e-5]<-0
  beta.opt<-svd1$v%*%diag(d.inv)%*%beta.opt
  beta.opt[-1]<-beta.opt[-1]/x.sd
  
  beta.opt[1]<-beta.opt[1]-x.mean%*%beta.opt[-1]
  
  
  Dx<-diag(c(1,1/x.sd))
  Dx.inv<-ginv(Dx)
  
  V<-gmm.func(beta.opt)$V
  
  XG.1<--X*(probs.opt)^.5*(1-probs.opt)^.5
  if(ATT==T){
  	XG.2<-X*((1-treat)*probs.opt/(1-probs.opt)*n/n.t)^.5
  } else{
  	XG.2<--X*(abs(probs.opt-treat)/(probs.opt*(1-probs.opt)))^.5
  }
    G<-cbind(t(XG.1)%*%XG.1,t(XG.2)%*%XG.2)
    var2<-ginv(G%*%ginv(V)%*%t(G))*n
    
    int.mat<-diag(dim(X)[2])

  #vcov2<-(diag(d.inv)%*%svd1$v%*%Dx.inv%*%ginv(opt1$hess)%*%Dx.inv%*%t(svd1$v)%*%diag(d.inv))
  vcov<-(Dx.inv)%*%svd1$v%*%diag(d.inv)%*%var2%*%diag(d.inv)%*%t(svd1$v)%*%(Dx.inv)


vcov[1,1]<-vcov[1,1]+x.mean^2%*%diag(vcov[-1,-1])
	
#	vcov<-diag(1/d.inv)%*%svd1$v%*%Dx.inv%*%ginv(opt1$hess)%*%Dx.inv%*%t(svd1$v)%*%diag(1/d.inv)

  #diag(vcov)<-diag(vcov)/c(1,x.sd)^2

  class(beta.opt) <- "coef"
  names(beta.opt) <- names.X
		
  output<-list("coefficients"=beta.opt,"residuals"=residuals,"fitted.values"=probs.opt,"rank"=k,"family"="CBPS",
               "deviance"=deviance,"weights"=w.opt,
               "y"=treat,"x"=X,"model"=NA,"converged"=opt1$conv,
               "data"=data, "J"=J.opt,"df"=k,"var"=vcov)
  
  class(output)<- c("CBPS","glm","lm")
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
	cat("J-Statistic:\t     ", format(signif(x$J)),"\n")
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
  rownames(coef.table)<-names.X
  
  pval <- coef.table[,4]
  symp <- symnum(pval, corr=FALSE,
                 cutpoints = c(0,  .001,.01,.05, .1, 1),
                 symbols = c("***","**","*","."," "))
  coef.print<-cbind(signif(coef.table,3),as.vector(symp))
  coef.print[coef.print=="0"]<-"0.000"
	
  cat("\nCall:  \n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
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
	#cat("\n    Null J:  ",object$J)
  cat("\nJ - statistic:  ",object$J," on ", object$df, " degrees of freedom\n")
  cat("\nLog-Likelihood: ",sum(object$y*log(object$fitted.values) + (1-object$y)*log(1-object$fitted.values)), "\n")
	
  out<-list("call"=object$call,"coefficients"=coef.table,"J"=object$J)
  invisible(out)
  
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


match.func<-function(probs,log.odds=F,treat,M=1){
    out1<-sapply(which(treat==1),FUN=function(x) {
				 dist<-abs(probs[treat==0]-probs[x])
				 if(log.odds==T) dist<-abs(log(probs[treat==0]/(1-probs[treat==0])) - log(probs[x]/(1-probs[x])))
				 which(dist<= sort(unique(dist),decreasing=F)[1:M])
				 })
    c(which(treat==1),which(treat==0)[unlist(out1)])
}


find.num.match<-function(fitted.values, y, X){
	
	object<-NULL
	object$fitted.values<-fitted.values
	object$y<-y
	object$x<-X
	
  match.func<-function(probs,log.odds=F,treat,M){
    out1<-sapply(which(treat==1),FUN=function(x) {
      dist<-abs(probs[treat==0]-probs[x])
      if(log.odds==T) dist<-abs(log(probs[treat==0]/(1-probs[treat==0])) - log(probs[x]/(1-probs[x])))
      sort(dist,decreasing=F,ind=T)$ix[1:M]  
    })
    c(which(treat==1),which(treat==0)[unlist(out1)])
  }
  
  match.check<-function(ind,probs,treat,X,ATT=T){
    X<-X[ind,]
    treat<-treat[ind]
    probs.curr<-probs[ind]
    
    ##Generates ATT weights 
    ATT.wt.func<-function(probs.curr,X.wt=X){
      X<-as.matrix(X.wt)
      n<-dim(X)[1]
      n.c<-sum(treat==0)
      n.t<-sum(treat==1)
      w1<-(n/n.t*(treat-probs.curr)/(1-probs.curr))
      w1[treat==1]<-n/n.t
      w1
    }
    
    n<-dim(X)[1]
    n.c<-sum(treat==0)
    n.t<-sum(treat==1)
    
    if(ATT) w.curr<-ATT.wt.func(probs.curr) else
    w.curr<-(probs.curr-1+treat)^-1
    w.curr.del<-1/(n)*t(X)%*%(w.curr)
    w.curr.del<-as.vector(w.curr.del)
    w.curr<-as.vector(w.curr)
    
    gbar<-c(1/n*t(X)%*%(treat-probs.curr),w.curr.del)
    
	  
	#Old version w analytic variances
#if(ATT){
#     X.1<-X*((1-probs.curr)*probs.curr)^.5 
#     X.2<-X*(probs.curr/(1-probs.curr))^.5
#     X.1.1<-X*(probs.curr)^.5}
#   else{
#      X.1<-X*((1-probs.curr)*probs.curr)^.5
#      X.2<-X*0#*(w.opt)  
#      X.1.1<-X*(1+probs.curr*(1-probs.curr))^-.5
#    }
#    V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.1), 1/n*cbind(t(X.1.1)%*%X.1.1,t(X.2)%*%X.2)) 

	  X.1<-X*(treat-probs.curr)
	  X.2<-X*w.curr
	  
	  X.1<-apply(X.1,2,FUN=function(x) x-mean(x))
	  X.2<-apply(X.2,2,FUN=function(x) x-mean(x))
	  
	  V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1)%*%X.2),
			   1/n*cbind(t(X.2)%*%X.1,t(X.2)%*%X.2))	
	  
	  
	  
	  
    loss1<-as.vector( t(gbar)%*%ginv(V)%*%(gbar))
    
    max(0,loss1*n.t)
  }
  
  m1.run <- array()
  m1.2 <- array()
  for(i in 1:10){
    ind.2<-match.func(object$fitted.values,log.odds=T,treat=object$y,M=i)
    m1<-match.check(ind=ind.2,object$fitted.values,treat=object$y,X=object$x,ATT=T)
    m1.run[i]<-m1
  }
  out<-NULL
  out$N<-m1.2[1]<-which(m1.run==min(m1.run))
  out$matches<-match.func(object$fitted.values,log.odds=T,treat=object$y,M=out$N)
  
  return(out)
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
