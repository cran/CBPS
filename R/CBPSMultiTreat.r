CBPS.3Treat<-function(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, ...)
{
	probs.min<-1e-6
	names.X<-colnames(X)
	if (is.null(iterations)) iterations<-10000
	
	no.treats<-length(levels(as.factor(treat)))
	treat.names<-levels(as.factor(treat))
	T1<-as.numeric(treat==treat.names[1])
	T2<-as.numeric(treat==treat.names[2])
	T3<-as.numeric(treat==treat.names[3])
	  
	##The gmm objective function--given a guess of beta, constructs the GMM J statistic.
	  gmm.func<-function(beta.curr,X.gmm=X){
		##Designate a few objects in the function.
		beta.curr<-matrix(beta.curr,k,no.treats-1)
		X<-as.matrix(X.gmm)
		
		##Designate sample size, number of treated and control observations,
		##theta.curr, which are used to generate probabilities.
		##Trim probabilities, and generate weights.
		n<-c(sum(T1), sum(T2), sum(T3))
		theta.curr<-X%*%beta.curr
		baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
		probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob)
		probs.curr[,1]<-pmin(1-probs.min,probs.curr[,1])
		probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
		probs.curr[,2]<-pmin(1-probs.min,probs.curr[,2])
		probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
		probs.curr[,3]<-pmin(1-probs.min,probs.curr[,3])
		probs.curr[,3]<-pmax(probs.min,probs.curr[,3])

		w.curr<-cbind(T1/probs.curr[,1] - T2/probs.curr[,2],
					  T1/probs.curr[,1] - T3/probs.curr[,3],
					  T2/probs.curr[,2] - T3/probs.curr[,3])
		  
	  
		##Generate the vector of mean imbalance by weights.
		w.curr.del<-1/sum(n)*t(X)%*%w.curr
		w.curr.del<-as.matrix(w.curr.del)
		w.curr<-as.matrix(w.curr)

		##Generate g-bar, as in the paper.
		gbar<-c(1/sum(n)*t(X)%*%(T2-probs.curr[,2])*(1-bal.only), 
				1/sum(n)*t(X)%*%(T3-probs.curr[,3])*(1-bal.only),
				w.curr.del)

		##Generate the covariance matrix used in the GMM estimate.
		##Was for the initial version that calculates the analytic variances.		
		X.1.1<-X*(probs.curr[,2]*(1-probs.curr[,2]))
		X.1.2<-X*(-probs.curr[,2]*probs.curr[,3])
		X.1.3<-X*-1
		X.1.4<-X*0
		X.1.5<-X*1
		X.2.2<-X*(probs.curr[,3]*(1-probs.curr[,3]))
		X.2.3<-X*0
		X.2.4<-X*-1
		X.2.5<-X*-1
		X.3.3<-X*(probs.curr[,1]^-1+probs.curr[,2]^-1)
		X.3.4<-X*probs.curr[,1]^-1
		X.3.5<-X*-probs.curr[,2]^-1
		X.4.4<-X*(probs.curr[,1]^-1+probs.curr[,3]^-1)
		X.4.5<-X*probs.curr[,3]^-1
		X.5.5<-X*(probs.curr[,2]^-1+probs.curr[,3]^-1)
		
		V<-1/sum(n)*rbind(cbind(t(X.1.1)%*%X,t(X.1.2)%*%X,t(X.1.3)%*%X,t(X.1.4)%*%X,t(X.1.5)%*%X),
						  cbind(t(X.1.2)%*%X,t(X.2.2)%*%X,t(X.2.3)%*%X,t(X.2.4)%*%X,t(X.2.5)%*%X),
						  cbind(t(X.1.3)%*%X,t(X.2.3)%*%X,t(X.3.3)%*%X,t(X.3.4)%*%X,t(X.3.5)%*%X),
						  cbind(t(X.1.4)%*%X,t(X.2.4)%*%X,t(X.3.4)%*%X,t(X.4.4)%*%X,t(X.4.5)%*%X),
						  cbind(t(X.1.5)%*%X,t(X.2.5)%*%X,t(X.3.5)%*%X,t(X.4.4)%*%X,t(X.5.5)%*%X))
		
		##Calculate the GMM loss.
		loss1<-as.vector(t(gbar)%*%ginv(V)%*%(gbar))      
		out1<-list("loss"=max(loss1*sum(n),loss1*sum(n)), "V"=V)
		out1
	  }
	  gmm.loss<-function(x,...) gmm.func(x,...)$loss
		
	  ##Loss function for balance constraints, returns the squared imbalance along each dimension.
	  bal.loss<-function(beta.curr){
		beta.curr<-matrix(beta.curr,k,no.treats-1)
		theta.curr<-X%*%beta.curr
		baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
		probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob)
		probs.curr[,1]<-pmin(1-probs.min,probs.curr[,1])
		probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
		probs.curr[,2]<-pmin(1-probs.min,probs.curr[,2])
		probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
		probs.curr[,3]<-pmin(1-probs.min,probs.curr[,3])
		probs.curr[,3]<-pmax(probs.min,probs.curr[,3])

		w.curr<-cbind(T1/probs.curr[,1] - T2/probs.curr[,2],
					  T1/probs.curr[,1] - T3/probs.curr[,3],
					  T2/probs.curr[,2] - T3/probs.curr[,3])
		X.2<-X
		##Generate mean imbalance.
		loss1<-sum(diag(abs(t(w.curr)%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr))))
		loss1
	  }
	  n<-length(treat)
	  n.t<-c(sum(T1),sum(T2),sum(T3))
	  x.orig<-x<-cbind(as.matrix(X))
  
	##Run multionmial logit
	  dat.dummy<-data.frame(treat=treat,X)
	  #Need to generalize for different dimensioned X's
	  xnam<- paste0("X", 1:dim(X)[2])
	  fmla <- as.formula(paste("as.factor(treat) ~ -1 + ", paste(xnam, collapse= "+")))
	  mnl1<-multinom(fmla, data=dat.dummy, trace=FALSE)
	  mcoef<-matrix(coef(mnl1),k,no.treats-1,byrow=TRUE)
	  mcoef[is.na(mcoef[,1]),1]<-0
	  mcoef[is.na(mcoef[,2]),2]<-0
	  probs.mnl<-cbind(1/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])),
				exp(X%*%mcoef[,1])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])),
				exp(X%*%mcoef[,2])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])))
	  colnames(probs.mnl)<-c("p1","p2","p3")
	  probs.mnl[,1]<-pmin(1-probs.min,probs.mnl[,1])
	  probs.mnl[,1]<-pmax(probs.min,probs.mnl[,1])
	  probs.mnl[,2]<-pmin(1-probs.min,probs.mnl[,2])
	  probs.mnl[,2]<-pmax(probs.min,probs.mnl[,2])
	  probs.mnl[,3]<-pmin(1-probs.min,probs.mnl[,3])
	  probs.mnl[,3]<-pmax(probs.min,probs.mnl[,3])
	  mnl1$fit<-matrix(probs.mnl,nrow=n,ncol=no.treats)
	  beta.curr<-mcoef
	  beta.curr[is.na(beta.curr)]<-0
  		
	  alpha.func<-function(alpha) gmm.loss(beta.curr*alpha)
	  beta.curr<-beta.curr*optimize(alpha.func,interval=c(.8,1.1))$min
		
	  ##Generate estimates for balance and CBPSE
	  gmm.init<-beta.curr
	  
	  opt.bal<-tryCatch({
							optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
						},
						error = function(err)
						{
							return(optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
						})
	  beta.bal<-opt.bal$par
	  
	  gmm.glm.init<-tryCatch({
								optim(mcoef,gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
							},
							error = function(err)
							{
								return(optim(mcoef,gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
							})
	  gmm.bal.init<-tryCatch({
								optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
							},
							error = function(err)
							{
								return(optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
							})
		
		if(gmm.glm.init$val<gmm.bal.init$val) opt1<-gmm.glm.init else opt1<-gmm.bal.init
		if(bal.only) opt1<-opt.bal 
				
	  ##Generate probabilities
		beta.opt<-matrix(opt1$par,nrow=k,ncol=no.treats-1)
		theta.opt<-X%*%beta.opt
		baseline.prob<-apply(theta.opt,1,function(x) (1+sum(exp(x)))^-1)
		probs.opt<-cbind(baseline.prob, exp(theta.opt[,1])*baseline.prob, exp(theta.opt[,2])*baseline.prob)
		probs.opt[,1]<-pmin(1-probs.min,probs.opt[,1])
		probs.opt[,1]<-pmax(probs.min,probs.opt[,1])
		probs.opt[,2]<-pmin(1-probs.min,probs.opt[,2])
		probs.opt[,2]<-pmax(probs.min,probs.opt[,2])
		probs.opt[,3]<-pmin(1-probs.min,probs.opt[,3])
		probs.opt[,3]<-pmax(probs.min,probs.opt[,3])
			  
	  J.opt<-gmm.func(beta.opt)$loss
	  
	  if ((J.opt > gmm.loss(mcoef)) & (bal.loss(beta.opt) > bal.loss(mcoef)))
		{	  
			beta.opt<-mcoef
			probs.opt<-probs.mnl
			J.opt <- gmm.loss(mcoef)
			warning("Optimization failed.  Results returned are for MLE.")
		}
		
	  ##Generate weights
	    w.opt<-rbind(T1/probs.opt[,1] - T2/probs.opt[,2],
					 T1/probs.opt[,1] - T3/probs.opt[,3],
					 T2/probs.opt[,2] - T3/probs.opt[,3])
	  
	  residuals<-cbind(T1-probs.opt[,1],T2-probs.opt[,2],T3-probs.opt[,3])
	  deviance <- -2*c(sum(T1*log(probs.opt[,1])+T2*log(probs.opt[,2])+T3*log(probs.opt[,3])))
	  nulldeviance <- -2*c(sum(T1*log(mean(T1))+T2*log(mean(T2))+T3*log(mean(T3))))
	  
	  V<-gmm.func(beta.opt)$V
	  X.G.1<-cbind(-X*probs.opt[,2]*(1-probs.opt[,2]),X*probs.opt[,2]*probs.opt[,3])
	  X.G.2<-cbind(X*probs.opt[,2]*probs.opt[,3],-X*probs.opt[,3]*(1-probs.opt[,3]))
	  X.G.3<-cbind(X*(T1*probs.opt[,2]/probs.opt[,1] + T2*(1-probs.opt[,2])/probs.opt[2]),X*probs.opt[,3]*(T1/probs.opt[,1] - T2/probs.opt[,2]))
	  X.G.4<-cbind(X*probs.opt[,2]*(T1/probs.opt[,1] - T3/probs.opt[,3]),X*(T1*probs.opt[,3]/probs.opt[,1] + T3*(1-probs.opt[,3]/probs.opt[,3])))
	  X.G.5<-cbind(-X*(T2*(1-probs.opt[,2])/probs.opt[,2] + T3/probs.opt[,3]),X*(T2*probs.opt[,3]/probs.opt[2] + T3*(1-probs.opt[,3])/probs.opt[,3]))
	  
	  G<-cbind(t(X.G.1)%*%X, t(X.G.2)%*%X, t(X.G.3)%*%X, t(X.G.4)%*%X, t(X.G.5)%*%X)
	  vcov<-ginv(G%*%ginv(V)%*%t(G))*n
	  
	  colnames(vcov)<-c(paste0("2:X", 1:k),paste0("3:X", 1:k))
	  rownames(vcov)<-colnames(vcov)
	  colnames(probs.opt)<-c("p1","p2","p3")

	  
	    class(beta.opt) <- "coef"
		names(beta.opt) <- names.X
		
		output<-list("coefficients"=beta.opt,"residuals"=residuals,"fitted.values"=probs.opt,"rank"=k,"family"="CBPS",
					 "deviance"=deviance,"weights"=w.opt,
					 "y"=treat,"x"=X,"model"=NA,"converged"=opt1$conv,
					 "data"=data, "J"=J.opt,"df"=k,"var"=vcov)
  
		class(output)<- c("CBPS")
		output
}

CBPS.4Treat<-function(treat, X, X.bal, method, k, XprimeX.inv, bal.only, iterations, ...)
{
	probs.min<-1e-6
	names.X<-colnames(X)
	if (is.null(iterations)) iterations<-10000

	no.treats<-length(levels(as.factor(treat)))
	treat.names<-levels(as.factor(treat))
	T1<-as.numeric(treat==treat.names[1])
	T2<-as.numeric(treat==treat.names[2])
	T3<-as.numeric(treat==treat.names[3])
	T4<-as.numeric(treat==treat.names[4])
	  
	##The gmm objective function--given a guess of beta, constructs the GMM J statistic.
	  gmm.func<-function(beta.curr,X.gmm=X){
		##Designate a few objects in the function.
		beta.curr<-matrix(beta.curr,k,no.treats-1)
		X<-as.matrix(X.gmm)
		
		##Designate sample size, number of treated and control observations,
		##theta.curr, which are used to generate probabilities.
		##Trim probabilities, and generate weights.
		n<-c(sum(T1), sum(T2), sum(T3), sum(T4))
		theta.curr<-X%*%beta.curr
		baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
		probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob, exp(theta.curr[,3])*baseline.prob)
		probs.curr[,1]<-pmin(1-probs.min,probs.curr[,1])
		probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
		probs.curr[,2]<-pmin(1-probs.min,probs.curr[,2])
		probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
		probs.curr[,3]<-pmin(1-probs.min,probs.curr[,3])
		probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
		probs.curr[,4]<-pmin(1-probs.min,probs.curr[,4])
		probs.curr[,4]<-pmax(probs.min,probs.curr[,4])

		w.curr<-cbind(T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] - T4/probs.curr[,4],
					  T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4],
					  -T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4])
		  
	  
		##Generate the vector of mean imbalance by weights.
		w.curr.del<-1/sum(n)*t(X)%*%w.curr
		w.curr.del<-as.matrix(w.curr.del)
		w.curr<-as.matrix(w.curr)

		##Generate g-bar, as in the paper.
		gbar<-c(1/sum(n)*t(X)%*%(T2-probs.curr[,2]), 
				1/sum(n)*t(X)%*%(T3-probs.curr[,3]),
				1/sum(n)*t(X)%*%(T4-probs.curr[,4]),
				w.curr.del)

		##Generate the covariance matrix used in the GMM estimate.
		##Was for the initial version that calculates the analytic variances.		
		X.1.1<-X*(probs.curr[,2]*(1-probs.curr[,2]))
		X.1.2<-X*(-probs.curr[,2]*probs.curr[,3])
		X.1.3<-X*(-probs.curr[,2]*probs.curr[,4])
		X.1.4<-X
		X.1.5<-X*(-1)
		X.1.6<-X
		X.2.2<-X*(probs.curr[,3]*(1-probs.curr[,3]))
		X.2.3<-X*(-probs.curr[,3]*probs.curr[,4])
		X.2.4<-X*(-1)
		X.2.5<-X*(-1)
		X.2.6<-X*(-1)
		X.3.3<-X*(probs.curr[,4]*(1-probs.curr[,4]))
		X.3.4<-X*(-1)
		X.3.5<-X
		X.3.6<-X
		X.4.4<-X*(probs.curr[,1]^-1+probs.curr[,2]^-1+probs.curr[,3]^-1+probs.curr[,4]^-1)
		X.4.5<-X*(probs.curr[,1]^-1-probs.curr[,2]^-1+probs.curr[,3]^-1-probs.curr[,4]^-1)
		X.4.6<-X*(-probs.curr[,1]^-1+probs.curr[,2]^-1+probs.curr[,3]^-1-probs.curr[,4]^-1)
		X.5.5<-X.4.4
		X.5.6<-X*(-probs.curr[,1]^-1-probs.curr[,2]^-1+probs.curr[,3]^-1+probs.curr[,4]^-1)
		X.6.6<-X.4.4
		
		V<-1/sum(n)*rbind(cbind(t(X.1.1)%*%X,t(X.1.2)%*%X,t(X.1.3)%*%X,t(X.1.4)%*%X,t(X.1.5)%*%X,t(X.1.6)%*%X),
						  cbind(t(X.1.2)%*%X,t(X.2.2)%*%X,t(X.2.3)%*%X,t(X.2.4)%*%X,t(X.2.5)%*%X,t(X.2.6)%*%X),
						  cbind(t(X.1.3)%*%X,t(X.2.3)%*%X,t(X.3.3)%*%X,t(X.3.4)%*%X,t(X.3.5)%*%X,t(X.3.6)%*%X),
						  cbind(t(X.1.4)%*%X,t(X.2.4)%*%X,t(X.3.4)%*%X,t(X.4.4)%*%X,t(X.4.5)%*%X,t(X.4.6)%*%X),
						  cbind(t(X.1.5)%*%X,t(X.2.5)%*%X,t(X.3.5)%*%X,t(X.4.5)%*%X,t(X.5.5)%*%X,t(X.5.6)%*%X),
						  cbind(t(X.1.6)%*%X,t(X.2.6)%*%X,t(X.3.6)%*%X,t(X.4.6)%*%X,t(X.5.6)%*%X,t(X.6.6)%*%X))
		
		##Calculate the GMM loss.
		loss1<-as.vector(t(gbar)%*%ginv(V)%*%(gbar))      
		out1<-list("loss"=max(loss1*sum(n),loss1*sum(n)), "V"=V)
		out1
	  }
	  gmm.loss<-function(x,...) gmm.func(x,...)$loss
		
	  ##Loss function for balance constraints, returns the squared imbalance along each dimension.
	  bal.loss<-function(beta.curr){
		beta.curr<-matrix(beta.curr,k,no.treats-1)
		theta.curr<-X%*%beta.curr
		baseline.prob<-apply(theta.curr,1,function(x) (1+sum(exp(x)))^-1)
		probs.curr<-cbind(baseline.prob, exp(theta.curr[,1])*baseline.prob, exp(theta.curr[,2])*baseline.prob, exp(theta.curr[,3]*baseline.prob))
		probs.curr[,1]<-pmin(1-probs.min,probs.curr[,1])
		probs.curr[,1]<-pmax(probs.min,probs.curr[,1])
		probs.curr[,2]<-pmin(1-probs.min,probs.curr[,2])
		probs.curr[,2]<-pmax(probs.min,probs.curr[,2])
		probs.curr[,3]<-pmin(1-probs.min,probs.curr[,3])
		probs.curr[,3]<-pmax(probs.min,probs.curr[,3])
		probs.curr[,4]<-pmin(1-probs.min,probs.curr[,4])
		probs.curr[,4]<-pmax(probs.min,probs.curr[,4])

		w.curr<-cbind(T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] - T4/probs.curr[,4],
					  T1/probs.curr[,1] - T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4],
					  -T1/probs.curr[,1] + T2/probs.curr[,2] - T3/probs.curr[,3] + T4/probs.curr[,4])
					  
		X.2<-X
		##Generate mean imbalance.
		loss1<-sum(diag(abs(t(w.curr)%*%X%*%XprimeX.inv%*%t(X)%*%(w.curr))))
		loss1
	  }
	  n<-length(treat)
	  n.t<-c(sum(T1),sum(T2),sum(T3),sum(T4))
	  x.orig<-x<-cbind(as.matrix(X))
  
	##Run multionmial logit
	  dat.dummy<-data.frame(treat=treat,X)
	  #Need to generalize for different dimensioned X's
	  xnam<- paste0("X", 1:dim(X)[2])
	  fmla <- as.formula(paste("as.factor(treat) ~ -1 + ", paste(xnam, collapse= "+")))
	  mnl1<-multinom(fmla, data=dat.dummy,trace=FALSE)
	  mcoef<-matrix(coef(mnl1),k,no.treats-1,byrow=TRUE)
	  mcoef[is.na(mcoef[,1]),1]<-0
	  mcoef[is.na(mcoef[,2]),2]<-0
	  mcoef[is.na(mcoef[,3]),3]<-0
	  probs.mnl<-cbind(1/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])+exp(X%*%mcoef[,3])),
				exp(X%*%mcoef[,1])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])+exp(X%*%mcoef[,3])),
				exp(X%*%mcoef[,2])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])+exp(X%*%mcoef[,3])),
				exp(X%*%mcoef[,3])/(1+exp(X%*%mcoef[,1])+exp(X%*%mcoef[,2])+exp(X%*%mcoef[,3])))
	  colnames(probs.mnl)<-c("p1","p2","p3","p4")
	  probs.mnl[,1]<-pmin(1-probs.min,probs.mnl[,1])
	  probs.mnl[,1]<-pmax(probs.min,probs.mnl[,1])
	  probs.mnl[,2]<-pmin(1-probs.min,probs.mnl[,2])
	  probs.mnl[,2]<-pmax(probs.min,probs.mnl[,2])
	  probs.mnl[,3]<-pmin(1-probs.min,probs.mnl[,3])
	  probs.mnl[,3]<-pmax(probs.min,probs.mnl[,3])
	  probs.mnl[,4]<-pmin(1-probs.min,probs.mnl[,4])
	  probs.mnl[,4]<-pmax(probs.min,probs.mnl[,4])
	  mnl1$fit<-matrix(probs.mnl,nrow=n,ncol=no.treats)
	  beta.curr<-mcoef
	  beta.curr[is.na(beta.curr)]<-0
	  
	  alpha.func<-function(alpha) gmm.loss(beta.curr*alpha)
	  beta.curr<-beta.curr*optimize(alpha.func,interval=c(.8,1.1))$min
		
	  ##Generate estimates for balance and CBPSE
	  gmm.init<-beta.curr
	 
	  opt.bal<-tryCatch({
							optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
						},
						error = function(err)
						{
							return(optim(gmm.init, bal.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
						})
	  beta.bal<-opt.bal$par
	  
	  gmm.glm.init<-tryCatch({
								optim(mcoef,gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
							},
							error = function(err)
							{
								return(optim(mcoef,gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
							})
	  gmm.bal.init<-tryCatch({
								optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="BFGS", hessian=TRUE)
							},
							error = function(err)
							{
								return(optim(beta.bal, gmm.loss, control=list("maxit"=iterations), method="Nelder-Mead", hessian=TRUE))
							})
							
		if(gmm.glm.init$val<gmm.bal.init$val) opt1<-gmm.glm.init else opt1<-gmm.bal.init
		if(bal.only) opt1<-opt.bal 
				
	  ##Generate probabilities
		beta.opt<-matrix(opt1$par,nrow=k,ncol=no.treats-1)
		theta.opt<-X%*%beta.opt
		baseline.prob<-apply(theta.opt,1,function(x) (1+sum(exp(x)))^-1)
		probs.opt<-cbind(baseline.prob, exp(theta.opt[,1])*baseline.prob, exp(theta.opt[,2])*baseline.prob, exp(theta.opt[,3])*baseline.prob)
		probs.opt[,1]<-pmin(1-probs.min,probs.opt[,1])
		probs.opt[,1]<-pmax(probs.min,probs.opt[,1])
		probs.opt[,2]<-pmin(1-probs.min,probs.opt[,2])
		probs.opt[,2]<-pmax(probs.min,probs.opt[,2])
		probs.opt[,3]<-pmin(1-probs.min,probs.opt[,3])
		probs.opt[,3]<-pmax(probs.min,probs.opt[,3])
		probs.opt[,4]<-pmin(1-probs.min,probs.opt[,4])
		probs.opt[,4]<-pmax(probs.min,probs.opt[,4])
	  
	  J.opt<-gmm.func(beta.opt)$loss
	  
	  if ((J.opt > gmm.loss(mcoef)) & (bal.loss(beta.opt) > bal.loss(mcoef)))
		{	  
			beta.opt<-mcoef
			probs.opt<-probs.mnl
			J.opt <- gmm.loss(mcoef)
			warning("Optimization failed.  Results returned are for MLE.")
		}
		
	  ##Generate weights
		w.opt<-cbind(T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] - T4/probs.opt[,4],
					  T1/probs.opt[,1] - T2/probs.opt[,2] - T3/probs.opt[,3] + T4/probs.opt[,4],
					  -T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] + T4/probs.opt[,4])
		
	  #How are residuals now defined?
	  residuals<-cbind(T1-probs.opt[,1],T2-probs.opt[,2],T3-probs.opt[,3],T4-probs.opt[,4])
	  deviance <- -2*c(sum(T1*log(probs.opt[,1])+T2*log(probs.opt[,2])+T3*log(probs.opt[,3])+T4*log(probs.opt[,4])))
	  nulldeviance <- -2*c(sum(T1*log(mean(T1))+T2*log(mean(T2))+T3*log(mean(T3))+T4*log(mean(T4))))
	  
	  V<-gmm.func(beta.opt)$V
	  X.G.1<-cbind(-X*probs.opt[,2]*(1-probs.opt[,2]),X*probs.opt[,2]*probs.opt[,3],X*probs.opt[,2]*probs.opt[,4])
	  X.G.2<-cbind(X*probs.opt[,2]*probs.opt[,3],-X*probs.opt[,3]*(1-probs.opt[,3]),X*probs.opt[,3]*probs.opt[,4])
	  X.G.3<-cbind(X*probs.opt[,2]*probs.opt[,4],X*probs.opt[3,]*probs.opt[,4],-X*probs.opt[,4]*(1-probs.opt[,4]))
	  X.G.4<-cbind(X*probs.opt[,2]*(T1/probs.opt[,1] - T2*(1-probs.opt[,2])/probs.opt[,2] - T3/probs.opt[,3] - T4/probs.opt[,4]),
				   X*probs.opt[,3]*(T1/probs.opt[,1] + T2/probs.opt[,2] + T3*(1-probs.opt[,3])/probs.opt[,3] - T4/probs.opt[,4]),
				   X*probs.opt[,4]*(T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] + T4*(1-probs.opt[,4])/probs.opt[,4]))
	  X.G.5<-cbind(X*probs.opt[,2]*(T1/probs.opt[,1] + T2*(1-probs.opt[,2])/probs.opt[,2] - T3/probs.opt[,3] + T4/probs.opt[,4]),
				   X*probs.opt[,3]*(T1/probs.opt[,1] - T2/probs.opt[,2] + T3*(1-probs.opt[,3])/probs.opt[,3] + T4/probs.opt[,4]),
				   X*probs.opt[,4]*(T1/probs.opt[,1] - T2/probs.opt[,2] - T3/probs.opt[,3] - T4*(1-probs.opt[,4])/probs.opt[,4]))
	  X.G.6<-cbind(X*probs.opt[,2]*(-T1/probs.opt[,1] - T2*(1-probs.opt[,2])/probs.opt[,2] - T3/probs.opt[,3] + T4/probs.opt[,4]),
				   X*probs.opt[,3]*(-T1/probs.opt[,1] + T2/probs.opt[,2] + T3*(1-probs.opt[,3])/probs.opt[,3] + T4/probs.opt[,4]),
				   X*probs.opt[,4]*(-T1/probs.opt[,1] + T2/probs.opt[,2] - T3/probs.opt[,3] - T4*(1-probs.opt[,4])/probs.opt[,4]))				   
				   
	  G<-cbind(t(X.G.1)%*%X, t(X.G.2)%*%X, t(X.G.3)%*%X, t(X.G.4)%*%X, t(X.G.5)%*%X, t(X.G.6)%*%X)
	  vcov<-ginv(G%*%ginv(V)%*%t(G))*n
	  
	  colnames(vcov)<-c(paste0("2:X",1:k),paste0("3:X",1:k),paste0("4:X",1:k))
	  rownames(vcov)<-colnames(vcov)
	  colnames(probs.opt)<-c("p1","p2","p3","p4")
	  
	  	class(beta.opt) <- "coef"
		names(beta.opt) <- names.X
		
		output<-list("coefficients"=beta.opt,"residuals"=residuals,"fitted.values"=probs.opt,"rank"=k,"family"="CBPS",
					 "deviance"=deviance,"weights"=w.opt,
					 "y"=treat,"x"=X,"model"=NA,"converged"=opt1$conv,
					 "data"=data, "J"=J.opt,"df"=k,"var"=vcov)
  
		class(output)<- c("CBPS")
		output
}

