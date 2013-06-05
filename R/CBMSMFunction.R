
library(MASS)


########################
###Calls loss function
########################
CBMSM.fit<-function(treat,X, time,method, MultiBin.fit=FALSE, ...){

	if(method=="over") bal.only=FALSE
	if(method=="exact") bal.only=TRUE

	X0<-X

	time0=time
	n.t<-sum(time==time[1])
	X<-X[,apply(X,2,sd)>0]
	X<-cbind(1,X)
	svd1<-x.sd<-NULL
	for(i in unique(time)){
		x.sd<-c(x.sd,1,apply(X[time==i,-1],2,FUN=function(x) sd(x)  ))
		X[time==i,-1]<-apply(X[time==i,-1],2,FUN=function(x) (x-mean(x))/sd(x)  )
		X[is.na(X)]<-0
		svd.X<-svd(X[time==i,])
		svd1$d<-c(svd1$d,svd.X$d)
		svd1$v<-cbind(svd1$v,svd.X$v)
		X[time==i,]<-svd.X$u%*%diag(abs(svd.X$d)>0.0001) 
		}
	
		time<-as.numeric(as.factor(time))	
		n<-sum(time==time[1])
		
		betas.glm<-as.vector(
							 sapply(sort(unique(time)),FUN=function(x) glm(treat~X-1,subset=which(time0==x),family="binomial")$coef)
							 )
		betas.glm[is.na(betas.glm)]<-0
		
		


		bin.obs.mat<-matrix(NA, nrow=dim(X)[1],ncol=length(unique(time)))

		for(i in unique(time)){
			bin.obs.mat[,i]<-treat[time==i]

		}
		twos<-2^(seq(1:length(unique(time)))-1)
		bin.obs<-bin.obs.mat%*%twos
		
	
		probs.uncond<-NULL
		for(i in unique(bin.obs)){
			probs.uncond[bin.obs==i]<-mean(bin.obs==i)
		
		}
		
		
		probs.uncond[is.na(probs.uncond)]<-0
		
		
		
		#treat<-c(treat.1,treat.2,treat.3)
		
		msm.loss.bal<-function(x) msm.loss.func(x,X=X,treat=treat,bal.only=TRUE,n=n,n.t=length(unique(time0)),time.loss=time0,MultiBin.loss=MultiBin.fit)
#msm.loss.bal(betas.glm)
		
		
		opt.bal<-optim(betas.glm,msm.loss.bal,method="BFGS",hessian=TRUE)
		
		
		if(bal.only=="FALSE"){
		if(msm.loss.func(as.vector(betas.glm),X=X,treat=treat,n.t=n.t,n=n,time.loss=time0,MultiBin.loss=MultiBin.fit)>
		   msm.loss.func(as.vector(opt.bal$par),X=X,treat=treat,n.t=n.t,n=n,time.loss=time0,MultiBin.loss=MultiBin.fit)
		   )
		opt1<-optim(opt.bal$par,msm.loss.func,method="BFGS",X=X,treat=treat,n=n,n.t=length(unique(time)),time.loss=time0,MultiBin.loss=MultiBin.fit,hessian=TRUE) else
			opt1<-optim(betas.glm,msm.loss.func,method="BFGS",X=X,treat=treat,n=n,n.t=length(unique(time)),time.loss=time0,MultiBin.loss=MultiBin.fit,hessian=TRUE)
		} else {
			opt1<-opt.bal
		}
		
		#wts.glm<-make.wts(betas.glm,X,treat)
		w.opt<-make.wts(opt1$par,X,treat,probs.uncond,time0,MultiBin.fit)
		probs.opt<-make.probs(opt1$par,X,treat,time=time0)
		#w.each<-matrix(treat/probs.opt+(1-treat)/(1-probs.opt),nrow=n.t,byrow=FALSE)
		#w.opt[["stabilized"]]<-apply(w.each,1,prod)*probs.uncond[1:n.t]
		#w.opt[["unstabilized"]]<-apply(w.each,1,prod)
		#w.opt[["unconditional"]]<-probs.uncond[1:n.t]
		#wts.bal<-make.wts(opt.bal$par,X,treat)
		beta.opt<-opt1$par
		residuals<-treat-probs.opt
		k<-round(sum(diag(t(X)%*%X%*%ginv(t(X)%*%X))))
		deviance <- -2*c(sum(treat*log(probs.opt)+(1-treat)*log(1-probs.opt)))
		nulldeviance <- -2*c(sum(treat*log(mean(treat))+(1-treat)*log(1-mean(treat))))	
		data<-X
		J.opt<-opt1$val
		var<-opt1$hessian
		jacob<-jacobian(msm.loss.func.g, opt1$par,X=X,treat=treat,n=n,n.t=length(unique(time)),time.loss=time0,MultiBin.loss=MultiBin.fit)
		V<-msm.loss.func.V(opt1$par,X=X,treat=treat,n=n,n.t=length(unique(time)),time.loss=time0,MultiBin.loss=MultiBin.fit)
		Omega0<-msm.loss.func.Omega(opt1$par,X=X,treat=treat,n=n,n.t=length(unique(time)),time.loss=time0,MultiBin.loss=MultiBin.fit)


	G<-jacob
	W0<-V
	W<-Omega<-matrix(0,nrow=dim(G)[1],ncol=dim(G)[1])
	k1<-dim(G)[1]
	k2<-dim(Omega0)[1]

	for(i.mat in 1:(k1/k2)) {
		Omega[(1:k2)+k2*(i.mat-1),(1:k2)+k2*(i.mat-1)]<-Omega0
		W[(1:k2)+k2*(i.mat-1),(1:k2)+k2*(i.mat-1)]<-W0
		}
		

	X2<-X0.2<-NULL
	for(i.t in unique(time)) {
		X2<-cbind(X2,X[time==i.t,])
		X0.2<-cbind(X0.2,cbind(X0[time==i.t,]))
	}


		variance<-ginv(t(G)%*%W%*%G)%*%t(G)%*%W%*%Omega%*%W%*%G%*%ginv(t(G)%*%W%*%G)/n^2
	
		svd.v<-NULL
		for(i.mat in 1:(dim(variance)[1]/dim(svd1$v)[1]))		svd.v<-rbind(svd.v,svd1$v)
		d2<-rep(svd1$d,dim(variance)[1]/length(svd1$d))
		
		var<-ginv(t(X0.2)%*%X0.2)%*%t(X0.2)%*%X2%*%t(svd.v)%*%ginv(diag(d2))%*%variance%*%ginv(diag(d2))%*%
		(svd.v)%*%t(X2)%*%X0.2%*%ginv(t(X0.2)%*%X0.2)
		
		svd.v.X<-matrix(0,nrow=length(opt1$par),ncol=length(opt1$par))
		k2<-length(opt1$par)/length(unique(time))
		for(i in 1:length(unique(time))) svd.v.X[(1:k2)+k2*(i-1),(1:k2)+k2*(i-1)]<-svd1$v[,(1:k2)+k2*(i-1)]
		
			beta.opt<-svd.v.X%*%ginv(diag(svd1$d))%*%opt1$par
	beta.opt<-beta.opt/x.sd
	beta.opt[is.na(beta.opt)]<-0
		 ##Output from CBPS object
		 #output<-list("coefficients"=beta.opt,"residuals"=residuals,"fitted.values"=probs.opt,"rank"=k,"family"="CBPS",
         #    "deviance"=deviance,"weights"=w.opt,
         #   "y"=treat,"x"=X,"model"=NA,"converged"=opt1$conv,
         #  "data"=data, "J"=J.opt,"df"=k,"var"=vcov)
         J.time<-length(unique(time))

		df.out<-k*(J.time*2^J.time- (2^(J.time-1)-1))
	output<-list("coefficients"=beta.opt,"residuals"=residuals,"fitted.values"=probs.opt,"rank"=k,"family"="CBMSM",
         "deviance"=deviance,"weights"=w.opt,"variance"=var,
         "y"=treat,"x"=X0,"model"=NA,"converged"=opt1$conv,
               "data"=data, "J"=J.opt,"df"=df.out,"var"=var,"time"=time)
  
  class(output)<- c("CBMSM","glm","lm")
  if(MultiBin.fit) class(output)<-c("CBMB","glm","lm")
  output
		
		
}

########################
###Makes weights and probabilities
########################

make.wts<-function(betas,X,treat,probs.uncond,time,MultiBin.fit){
treat.use<-betas.use<-thetas.use<-NULL
for(i in sort(unique(time))){
thetas.use<-cbind(thetas.use,X[time==i,]%*%betas[1:dim(X)[2]+(i-1)*dim(X)[2] ])
treat.use<-cbind(treat.use,treat[time==i])
	}
	probs.use<-(1+exp(-thetas.use))^-1
	probs.trim<-.0001
	probs.use<-apply(probs.use,2,FUN=function(x) pmax(x,probs.trim))
	probs.use<-apply(probs.use,2,FUN=function(x) pmin(x,1-probs.trim))


probs.obs<-treat.use/probs.use+(1-treat.use)/(1-probs.use)#+(treat.use-probs.use)^2/(probs.use*(1-probs.use))
	wts<-list()
	wts[["stabilized"]]<-apply(probs.obs,1,prod)*probs.uncond[time==1]
	wts[["unstabilized"]]<-apply(probs.obs,1,prod)
	wts[["unconditional"]]<-probs.uncond[time==1]
	wts
}



make.probs<-function(betas,X,treat,time){
treat.use<-betas.use<-thetas.use<-NULL
for(i in sort(unique(time))){
thetas.use<-c(thetas.use,X[time==i,]%*%betas[1:dim(X)[2]+(i-1)*dim(X)[2] ])
	}
	probs.use<-(1+exp(-thetas.use))^-1
	probs.trim<-.0001
	probs.use<-pmax(probs.use,probs.trim)
	probs.use<-pmin(probs.use,1-probs.trim)
probs.use

}


########################
###Loss function for MSM
########################
		
msm.loss.func.0<-function(betas,X=X,treat=treat,bal.only=F,time.sub=0,n.t=n.t,n=n,time.loss, MultiBin.loss){

	time<-as.vector(time.loss)
	unique.time<-sort(unique(time.loss))
	treat.use<-betas.use<-NULL
	n.t<-length(unique(time.loss))
	for(i in 1:n.t){
	betas.use<-cbind(betas.use,betas[1:dim(X)[2]+(i-1)*dim(X)[2]  ])
	treat.use<-cbind(treat.use,treat[time==unique.time[i]])
	}
	betas<-betas.use
	betas[is.na(betas)]<-0
	treat<-treat.use
	thetas<-NULL
	for(i in 1:n.t)
	thetas<-cbind(thetas,X[time==i,]%*%betas[,i] )

	probs.trim<-.0001
	probs<-(1+exp(-thetas))^(-1)
	probs<-pmax(probs,probs.trim)
	probs<-pmin(probs,1-probs.trim)
	probs.obs<-treat*probs+(1-treat)*(1-probs)



	w.each<-treat/probs+(1-treat)/(1-probs)#+(treat-probs)^2/(probs*(1-probs))
	w.all<-apply(w.each,1,prod)#*probs.uncond

	bin.mat<-matrix(0,nrow=(2^n.t-1),ncol=n.t)
	for(i in 1:(2^n.t-1)) bin.mat[i,(n.t-length(integer.base.b(i))+1):n.t]<-
	integer.base.b(i)

    #if(MultiBin.loss) bin.mat<-bin.mat*0+1
    
	num.valid.outer<-constr.mat.outer<-NULL
	is.valid.outer<-0
	for(i.time in 1:n.t){
		num.valid<-rep(0,dim(treat)[1])
		constr.mat.prop<-constr.mat<-matrix(0,nrow=dim(treat)[1],ncol=dim(bin.mat)[1])
		for(i in 1:dim(bin.mat)[1]){
			is.valid<-sum(bin.mat[i,(i.time):dim(bin.mat)[2]])>0
            if(MultiBin.loss) is.valid<-T
			if(is.valid){
				#for(i.wt in i.time:n.t) w.all.now<-w.all.now*1/(1+3*probs[,i.wt]*(1-probs[,i.wt]))
				constr.mat[,i]<-(w.all*(-1)^(treat%*%bin.mat[i,]))
				num.valid<-num.valid+1
				}else{
				constr.mat[,i]<-0
				}

			}	
			num.valid.outer<-c(num.valid.outer,num.valid)
			constr.mat.outer<-rbind(constr.mat.outer,constr.mat)
		}
		
		if(MultiBin.loss) outer.count<-rep(unique(num.valid.outer),n.t) else
			outer.count<-unique(num.valid.outer)

	X.wt<-X.prop<-g.wt<-g.prop<-omega.prop<-omega.wt<-NULL
	for(i in 1:n.t){
		g.prop<-c(g.prop, 1/n*t(X[time==i,])%*%(treat[,i]-probs[,i]))
		#g.prop<-rbind(g.prop,1/n*t(X[time==i,])%*%cbind((treat[,i]-probs[,i])/probs.obs[,i]*probs.uncond*abs(constr.mat.outer[time==i,]))*(i>time.sub))
		g.wt<-rbind(g.wt,1/n*t(X[time==i,])%*%cbind(constr.mat.outer[time==i,])*(i>time.sub))
		X.prop.curr<-matrix(0,ncol=n,nrow=dim(X)[2])
		X.wt.curr<-matrix(0,ncol=n,nrow=dim(X)[2])
		X.prop<-rbind(X.prop,1/n^.5*t((X[time==i,]*(probs.obs[,i]*(1-probs.obs[,i]))^.5)))
		X.wt<-rbind(X.wt,1/n^.5*t(X[time==i,]*w.all^.5*outer.count[i]^.5))

	}
	mat.prop<-matrix(0,nrow=n, ncol=dim(X.wt)[2])
	mat.prop[,1]<-1
	g.prop.all<-0*g.wt
	g.prop.all[,1]<-g.prop
	if(bal.only==T) g.prop.all<-0*g.prop.all

	g.all<-rbind(g.prop.all,g.wt)
	X.all<-rbind(X.prop*(1-bal.only),X.wt)
	var.X<-(X.all)%*%t(X.all)
	length.zero<-dim(g.prop.all)[2]#length(g.prop)

	
	Omega<-0
	for(j in 1:dim(constr.mat.outer)[2]){
	XG.1<-X*(as.vector(treat-probs))*1/dim(constr.mat.outer)[2]^.5
	XW.1<-X*as.vector(constr.mat.outer[,j])
	
	XG.2<-XW.2<-NULL
	for(i in 1:n.t){
		XG.2<-cbind(XG.2,XG.1[time==i,])
		XW.2<-cbind(XW.2,XW.1[time==i,])
	}	
		Omega<-Omega+rbind(t(XG.2),t(XW.2))%*%t(rbind(t(XG.2),t(XW.2)))
	}
	
	
	loss<-t(g.all)%*%ginv(var.X)%*%g.all
	out<-list("loss"=(sum(diag(loss)))*n, "g"=g.all, "V"=var.X, "Omega"=Omega,"df"=sum(num.valid.outer))
	out
	
}#closes msm.loss.func

msm.loss.func<-function(betas,...) msm.loss.func.0(betas,...)$loss
msm.loss.func.g<-function(betas,...) msm.loss.func.0(betas,...)$g
msm.loss.func.V<-function(betas,...) msm.loss.func.0(betas,...)$V
msm.loss.func.wt<-function(betas,...) msm.loss.func.0(betas,...)$wt.constr
msm.loss.func.Omega<-function(betas,...) msm.loss.func.0(betas,...)$Omega
msm.loss.func.df<-function(betas,...) msm.loss.func.0(betas,...)$df

########################
###Makes binary representation
########################
integer.base.b <-
function(x, b=2){
        xi <- as.integer(x)
        if(any(is.na(xi) | ((x-xi)!=0)))
                print(list(ERROR="x not integer", x=x))
        N <- length(x)
        xMax <- max(x)
        ndigits <- (floor(logb(xMax, base=2))+1)
        Base.b <- array(NA, dim=c(N, ndigits))
        for(i in 1:ndigits){#i <- 1
                Base.b[, ndigits-i+1] <- (x %% b)
                x <- (x %/% b)
        }
        if(N ==1) Base.b[1, ] else Base.b
} 


MSMdata<-function(n){
	nX<-4
	
beta.true<-rep(c(1,-.5,.25,.1,-1),length=4*nX+4)
#beta.true[c(1:3*(nX+1))]<-2
beta.true<-matrix(beta.true,nrow=nX+1,ncol=4)
beta.Y<-rep(c(27.4,13.7,13.7,13.7,0),times=3)
#beta.Y[c((1:4)*(nX+1))]<-0

#beta.true[10,]<-beta.true[10,]*2
time<-rep(1:3,each=n)
unique.time<-unique(time)
n.t<-length(unique(time))


X1<-matrix(0,nrow=n,ncol=nX+1)
X1[,1:nX]<-matrix(rnorm(nX*n),nrow=n,ncol=nX)
X1[,3:4]<-abs(X1[,3:4])
probs.true.1<-((1+exp(-X1%*%beta.true[,1]-.5))^-1)
treat.1<-rbinom(probs.true.1,size=1,n=n)


X2<-matrix(0,nrow=n,ncol=nX+1)
X2[,1:nX]<-matrix(rnorm(nX*n),nrow=n,ncol=nX)*(2+(2*treat.1-1))/(3)
X2[,nX+1]<-treat.1
X2[,3:4]<-abs(X2[,3:4])
probs.true.2<-((1+exp(-(X2)%*%beta.true[,2]+.5))^-1)
treat.2<-rbinom(probs.true.2,size=1,n=n)



X3<-matrix(0,nrow=n,ncol=nX+1)
X3[,1:nX]<-matrix(rnorm(nX*n),nrow=n,ncol=nX) *(2+(2*treat.2-1))/(3)
X3[,nX+1]<-treat.2
X3[,3:4]<-abs(X3[,3:4])
probs.true.3<-((1+exp(-(X3%*%beta.true[,3]-.5)))^-1)
treat.3<-rbinom(probs.true.3,size=1,n=n)



probs.joint<- (probs.true.1*treat.1+(1-probs.true.1)*(1-treat.1))*
(probs.true.2*treat.2+(1-probs.true.2)*(1-treat.2))*
(probs.true.3*treat.3+(1-probs.true.3)*(1-treat.3))

probs.uncond<- (mean(treat.1)*treat.1+(1-mean(treat.1))*(1-treat.1))*
(mean(treat.2)*treat.2+(1-mean(treat.2))*(1-treat.2))*
(mean(treat.3)*treat.3+(1-mean(treat.3))*(1-treat.3))

bin.obs<-4*treat.1+2*treat.2+1*treat.3

probs.uncond.0<-(table(bin.obs)/length(bin.obs))

probs.uncond<-NULL
for(i in 0:7){
	probs.uncond[bin.obs==i]<-probs.uncond.0[i+1]
	
}


#y<-2*treat.1+3*treat.2-3*treat.3+
#((cbind(X1[,-(nX+1)],X2[,-(nX+1)],X3[,-(nX+1)])))%*%as.vector(-beta.true)[-c(1:3*(nX+1))]*2
y<-210+30*(treat.1+treat.2+treat.3==1)+20*(treat.1+treat.2+treat.3==2)+10*(treat.1+treat.2+treat.3==3)+ cbind(X1,X2,X3)%*%as.vector(beta.Y)
y<-y+rnorm(length(y),sd=5)


output<-list("treat.1"=treat.1,"treat.2"=treat.2,"treat.3"=treat.3,"X1"=X1,"X2"=X2,"X3"=X3,
"probs.joint"=probs.joint,"probs.uncond"=probs.uncond,"y"=y,"n.t"=length(unique(time)))


}