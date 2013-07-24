CBPS.2Treat<-function(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, ATT, standardize){
	probs.min<- 1e-6
	if (is.null(iterations)) iterations<-1000
	
	treat<-sapply(treat,function(x) ifelse(x==TRUE,1,ifelse(x==FALSE,0,-1)))
	
	if (min(treat) < 0)
	{
		stop("For binary, treatments must be TRUE/FALSE or 1/0")
	}
	
  ##Generates ATT weights.	Called by loss function, etc.
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
	gbar<-c( 1/n*t(X)%*%(treat-probs.curr),w.curr.del)

	##Generate the covariance matrix used in the GMM estimate.
	##Was for the initial version that calculates the analytic variances.

	if(ATT){
	  X.1<-X*((1-probs.curr)*probs.curr)^.5
	  X.2<-X*(probs.curr/(1-probs.curr))^.5
	  X.1.1<-X*(probs.curr)^.5

	}else{
	  X.1<-X*((1-probs.curr)*probs.curr)^.5
	  X.2<-X*(probs.curr*(1-probs.curr))^-.5		
	  X.1.1<-  X
	}
	

	  V<-rbind(1/n*cbind(t(X.1)%*%X.1,t(X.1.1)%*%X.1.1)*n/sum(treat),
			   1/n*cbind(t(X.1.1)%*%X.1.1*n/sum(treat),t(X.2)%*%X.2*n^2/sum(treat)^2)
			   )		
	
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
	loss1<-abs(t(w.curr)%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr))
	loss1
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
	
  alpha.func<-function(alpha) gmm.loss(beta.curr*alpha)
  beta.curr<-beta.curr*optimize(alpha.func,interval=c(.8,1.1))$min
  
	
  ##Generate estimates for balance and CBPSE
  gmm.init<-beta.curr
 
  opt.bal<-optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
  beta.bal<-opt.bal$par
  
  if(bal.only) opt1<-opt.bal
  
  if(!bal.only)
  {
	gmm.glm.init<-optim(glm1$coef,gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
	gmm.bal.init<-optim(beta.bal,gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
	
	if(gmm.glm.init$val<gmm.bal.init$val) opt1<-gmm.glm.init else opt1<-gmm.bal.init
	}
	
	

  ##Generate probabilities
	beta.opt<-opt1$par
	theta.opt<-as.vector(X%*%beta.opt)
  probs.opt<-(1+exp(-theta.opt))^-1
  probs.opt<-pmin(1-probs.min,probs.opt)
  probs.opt<-pmax(probs.min,probs.opt)
	
  ##Generate weights
  if(ATT){
	w.opt<-abs(ATT.wt.func(beta.opt)) 
  }else{
	w.opt<-abs((probs.opt-1+treat)^-1)
  }
  
  norm1<-norm2<-1
  if (standardize)
  {
  	if (ATT)
  	{
  		norm1<-n/sum(treat==1)
  		norm2<-n/sum(treat==1)*sum((1-treat)*probs.opt/(1-probs.opt))
  	}
  	else
  	{
  		norm1<-sum(treat/probs.opt)
  		norm2<-sum((1-treat)/(1-probs.opt))
  	}
  }
  w.opt[which(treat==1)]<-w.opt[which(treat==1)]/norm1
  w.opt[which(treat==0)]<-w.opt[which(treat==0)]/norm2
  
  
  J.opt<-gmm.func(beta.opt)$loss
  
  residuals<-treat-probs.opt
  deviance <- -2*c(sum(treat*log(probs.opt)+(1-treat)*log(1-probs.opt)))
  nulldeviance <- -2*c(sum(treat*log(mean(treat))+(1-treat)*log(1-mean(treat))))

  XG.1<- -X*(probs.opt)^.5*(1-probs.opt)^.5
  XW.1<- X*(treat-probs.opt)
  if(ATT==T){
  	XW.2<-X*(treat-probs.opt)/(1-probs.opt)*n/n.t
	XG.2<-X*((1-treat)*probs.opt/(1-probs.opt)*n/n.t)^.5
  } else{
  	XW.2<- X*(probs.opt-1+treat)^-1
	XG.2<- -X*probs.opt^.5*(1-probs.opt)^.5*abs((probs.opt-1+treat)^-1)#*(abs(probs.opt-treat)/(probs.opt*(1-probs.opt)))^.5
  }
   W<-ginv(gmm.func(beta.opt)$V)
  	W1<-rbind(t(XW.1),t(XW.2))
  	Omega<-(W1%*%t(W1)/n)
	G<-cbind(t(XG.1)%*%XG.1,t(XG.2)%*%XG.2)/n
	vcov<-ginv(G%*%W%*%t(G))%*%G%*%W%*%Omega%*%W%*%t(G)%*%ginv(G%*%W%*%t(G))

	

	beta.opt<-opt1$par
  #class(beta.opt) <- "coef"
  #names(beta.opt) <- names.X
		
  output<-list("coefficients"=beta.opt,"residuals"=residuals,"fitted.values"=probs.opt,"rank"=k,"family"="CBPS",
			   "deviance"=deviance,"weights"=w.opt,
			   "y"=treat,"x"=X,"model"=NA,"converged"=opt1$conv,
			   "data"=data, "J"=J.opt,"df"=k,"var"=vcov)
  
  class(output)<- c("CBPS","glm","lm")
  output
  }
