
CBPSlogit <- function(formula, data, na.action, ATT, k = 0, bayes = FALSE, ...) {

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
    
  fit <- eval(call("CBPSlogit.fit", X = X, treat = Y, ATT=ATT, family = family, k=k, bayes=bayes,
                   intercept = attr(mt, "intercept") > 0L))
  ##if (model) 
  fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  xlevels <- .getXlevels(mt, mf)
  fit$call <- call
  fit$formula <- formula
  fit$terms<-mt
  fit
}

CBPSlogit.fit<-function(treat, X, ATT, X.bal=X, k, bayes, ...){
  ##A tolerance, to set probabilities to zero.
  probs.min<-1e-10
  names.X<-colnames(X)
  
  ##Generates ATT weights.  Called by loss function, etc.
  ATT.wt.func<-function(beta.curr,X.wt=X){
    X<-as.matrix(X.wt)
    n<-dim(X)[1]
    n.c<-sum(treat==0)
    n.t<-sum(treat==1)
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)	
    w1<-(n/n.t*(treat-probs.curr)/(1-probs.curr))
    w1[treat==1]<-n/n.t
    w1
  }
  
  ##The gmm objective function--given a guess of beta, constructs the GMM J statistic.
  gmm.loss<-function(beta.curr,X.gmm=X,ATT.gmm=ATT){
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
    probs.curr<-(1+exp(-theta.curr))^-1
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
    gbar<-c(1/n*t(X)%*%(treat-probs.curr),w.curr.del)
    
    ##Generate the covariance matrix used in the GMM estimate.
    if(ATT){
      X.1<-X*((1-probs.curr)*probs.curr)^.5 
      X.2<-X*(probs.curr/(1-probs.curr))^.5
      X.1.1<-X*(probs.curr)^.5      
    }else{
      X.1<-X*((1-probs.curr)*probs.curr)^.5
      X.2<-X*0#*(w.opt)		
      X.1.1<-X*(1+probs.curr*(1-probs.curr))^-.5
    }
		
    V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.1),
             1/n*cbind(t(X.1.1)%*%X.1.1,t(X.2)%*%X.2))	
    
    ##Calculate the GMM loss.
    loss1<-as.vector( t(gbar)%*%ginv(V)%*%(gbar))      
    max(0,loss1*n)
    
  }
	
  ##Loss function for balance constraints, returns the squared imbalance along each dimension.
  bal.loss<-function(beta.curr){
    ##Generate theta and probabilities.
    theta.curr<-as.vector(X%*%beta.curr)
    probs.curr<-(1+exp(-theta.curr))^-1
    probs.curr<-pmin(1-probs.min,probs.curr)
    probs.curr<-pmax(probs.min,probs.curr)
    ##Generate weights.
    if(ATT)
      w.curr<-ATT.wt.func(beta.curr)
    else
      w.curr<-(probs.curr-1+treat)^-1
    X.2<-X
    ##Generate mean imbalance.
    ##loss1<-abs(t(w.curr)%*%XprimeX.inv%*%(w.curr))
    loss1<-0
    loss1
  }
	
#######Declare some constants and orthogonalize Xdf.
  X0<-X
  
  format.bal<-F
  if(sum(dim(X.bal)[2]==dim(X)[2])) format.bal<-T
  x.sd<-apply(X[,-1],2,sd)
  x.mean<-apply(X[,-1],2,mean)
  X[,-1]<-apply(X[,-1],2,FUN=function(x) (x-mean(x))/sd(x))
  ##XprimeX.inv<-X%*%ginv(t(X)%*%X)%*%t(X)
  if(k==0) k<-sum(diag(t(X)%*%X%*%ginv(t(X)%*%X)))
  k<-floor(k+.1)
  svd1<-svd(X)
  X<-svd1$u[,1:k]
  
  ##XprimeX.inv<-X%*%ginv(t(X)%*%X)%*%t(X)

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
  glm1<-glm(treat~X-1,family=binomial)
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
  if(gmm.loss(beta.curr)>gmm.loss(beta.curr*0)) gmm.init<-beta.curr*0
  opt1<-optim(gmm.init, gmm.loss, control=list("maxit"=1000), method="BFGS", hessian=TRUE)

  ##Generate probabilities
  beta.opt<-opt1$par
  theta.opt<-as.vector(X%*%beta.opt)
  probs.opt<-(1+exp(-theta.opt))^-1
  probs.opt<-pmin(1-probs.min,probs.opt)
  probs.opt<-pmax(probs.min,probs.opt)
  
  ##beta.bal<-opt.bal$par
  ##theta.bal<-as.vector(X%*%beta.bal)
  ##probs.bal<-(1+exp(-theta.bal))^-1
  ##probs.bal<-pmin(1-probs.min,probs.bal)
  ##probs.bal<-pmax(probs.min,probs.bal)
	
  ##Generate weights
  if(ATT){
    w.opt<-ATT.wt.func(beta.opt) 
  }else{
    w.opt<-(probs.opt-1+treat)^-1
  }
  
  J.opt<-opt1$val*n
  
  residuals<-treat-probs.opt
  deviance <- -2*c(sum(treat*log(probs.glm)+(1-treat)*log(1-probs.glm)))
  nulldeviance <- -2*c(sum(treat*log(mean(treat))+(1-treat)*log(1-mean(treat))))		
  
  d.inv<-svd1$d
  d.inv[d.inv> 1e-5]<-1/d.inv[d.inv> 1e-5]
  d.inv[d.inv<= 1e-5]<-0
  beta.opt<-svd1$v%*%diag(d.inv)%*%beta.opt
  beta.opt[-1]<-beta.opt[-1]/x.sd
  
  Dx<-diag(c(1,1/x.sd))
  Dx.inv<-ginv(Dx)
  
  vcov<-ginv(diag(d.inv)%*%svd1$v%*%Dx.inv%*%opt1$hess%*%Dx.inv%*%t(svd1$v)%*%diag(d.inv))
  ##diag(vcov[-1,-1])<-diag(vcov[-1,-1])/x.sd^2

  class(beta.opt) <- "coef"
  names(beta.opt) <- names.X
		
  output<-list("coefficients"=beta.opt,"residuals"=residuals,"fitted.values"=probs.opt,"rank"=k,"family"="CBPS",
               "deviance"=deviance,"null.deviance"=NA,"weights"=w.opt,
               "df.residual"=1,"df.null"=NA,"y"=treat,"x"=X,"model"=NA,"converged"=opt1$conv,
               "data"=data, "J"=J.opt,"df"=k,"var"=vcov,"sd.coef"=diag(vcov)^.5)
  
  class(output)<- c("CBPSlogit","glm","lm")
  output	
}

print.CBPSlogit <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), 
        "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients:\n")
        print.default(format(x$coefficients, digits = digits), 
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", 
        x$df.residual, "Residual\n")
    if (nzchar(mess <- naprint(x$na.action))) 
        cat("  (", mess, ")\n", sep = "")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, 
        digits)), "\nResidual Deviance:", format(signif(x$deviance, 
        digits)), "\n")
	cat("J-Statistic:\t     ", format(signif(x$J)))
	cat("\n")
    invisible(x)
}

summary.CBPSlogit<-function(object, ...){
  ##x <- summary.glm(object, dispersion = dispersion, correlation = correlation, symbolic.cor = symbolic.cor, ...)
  
  x<-NULL
  names.X<-as.vector(names(object$coefficients))
  coef.table<-(cbind(as.vector(object$coefficients),as.vector(object$sd.coef),as.vector(object$coefficients/object$sd.coef),as.vector(2-2*pnorm(abs(object$coefficients/object$sd.coef)))))
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
	
  out<-list("call"=object$call,"coefficients"=coef.table,"J"=object$J)
  invisible(out)
  
}

#plot.CBPSlogit<-function(x,...){ 
#	plot(1)
#}

TE.est<-function(dv, object, M=1){
	
  match.func<-function(probs,log.odds=F,treat,M){
    out1<-sapply(which(treat==1),FUN=function(x) {
      dist<-abs(probs[treat==0]-probs[x])
      if(log.odds==T) dist<-abs(log(probs[treat==0]/(1-probs[treat==0])) - log(probs[x]/(1-probs[x])))
      sort(dist,decreasing=F,ind=T)$ix[1:M]  
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



find.num.match<-function(object){
  
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
    
    if(ATT){
      X.1<-X*((1-probs.curr)*probs.curr)^.5 
      X.2<-X*(probs.curr/(1-probs.curr))^.5
      X.1.1<-X*(probs.curr)^.5}
    else{
      X.1<-X*((1-probs.curr)*probs.curr)^.5
      X.2<-X*0#*(w.opt)  
      X.1.1<-X*(1+probs.curr*(1-probs.curr))^-.5
    }
    
    V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.1), 1/n*cbind(t(X.1.1)%*%X.1.1,t(X.2)%*%X.2)) 
    
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
