
# https://rpubs.com/aaronsc32/newton-raphson-method

library(rootSolve)
#library(nleqslv)

myfun<-function(f){
  
  # defining parameters
  phi<-5
  aM<-0.1
  aN<-0.2
  d<-0.5
  bmax<-0.8
  s<-0.1
  e<-0.5
  u<-0.4
  ps<-0.3
  K<-10
  KA<-5
  
  # get Req
  Req1<-((1-f)*K*d)/(e*(aN-aM*(1-f)))
  Req2<-(1/((1-f)*(bmax-d))) - (1/((bmax*(1-s))-d))
  Req<-Req1*Req2
  
  
  yeq<- (K*d)/((bmax*(1-s))-d)-(e*aM*Req)
  tempo<-(1-ps)/(u*yeq)
  
  Meq<-0.5*(tempo+sqrt((tempo^2)+(4*tempo*KA)))
  
  # solve the equantion when Neq = 0 for f=fmax
  Neq<-((phi/Req) -(aM*Meq))/aN
  
  return(list(Req=Req,
              Meq=Meq,
              Neq=Neq))
  #return(Neq) # to get uniroot solution
}

d<-0.5
bmax<-0.8
s<-0.1
fmin<-(s*bmax)/(bmax-d)
fmin
curve(myfun, xlim=c(fmin,1), col='blue', lwd=2, lty=2, ylab='Neq',xlab='f',ylim=c(-0.1,1))
abline(h=0)
uniroot(myfun, c(0.4,0.8)) #0.4947031


ans<-myfun(f=0.3)
(Meq<-ans$Meq)
(Neq<-ans$Neq)
(Req<-ans$Req)

Meqs<-c()
Neqs<-c()
frange<-seq(from=fmin,to=1,by=0.01)
for(f in frange){
  #cat("f=",f,"\n")
  ans<-myfun(f=f)
  Meqs<-c(Meqs,ans$Meq)
  Neqs<-c(Neqs,ans$Neq)
}

ylm<-round(max(Meqs,Neqs[-1]))
plot(frange,Meqs,type="l",ylab="",xlab="f",ylim=c(-1,ylm),xlim=c(fmin,1))
lines(frange,Neqs,lty="dashed")
abline(h=0,col="red")



