runall<-FALSE
if(runall==TRUE){
#load("../Data/BlackwellData")
form1<-"d.gone.neg ~ d.gone.neg.l1 + d.gone.neg.l2 + d.neg.frac.l3 + camp.length + camp.length + deminc + base.poll + as.factor(year) +   base.und + office"

#write.table(file="../Data/Blackwell.tab",blackwell)
blackwell<-read.table("../Data/Blackwell.tab")

numer<-form1
X0<-model.frame(form1,data=blackwell)
X.mat<-model.matrix(X0,data=X0)
X.mat<-X.mat[,apply(X.mat,2,sd)>0]
id<-as.numeric(as.factor(blackwell$demName))

treat<-blackwell$d.gone.neg
time<-blackwell$time



fit1<-CBMSM(formula = numer, time=time,id=id,data=blackwell, type="MSM",  iterations = NULL, twostep = TRUE, msm.variance = "approx", time.vary = FALSE)
}
