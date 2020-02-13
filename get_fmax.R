
# https://rpubs.com/aaronsc32/newton-raphson-method

library(rootSolve)
#library(nleqslv)

myfun<-function(f,s=0.1,aM=0.1,aN=0.2){
  
  # defining parameters
  phi<-5
  #aM<-0.1 
  #aN<-0.2
  d<-0.5
  bmax<-0.8
  #s<-0.1
  e<-0.5
  u<-0.4
  ps<-0.3
  K<-10
  KA<-5
  
  # get Req
  Req1<-((1-f)*K*d)/(e*(aN-aM*(1-f)))
  Req2<-(1/((1-f)*(bmax-d))) - (1/((bmax*(1-s))-d))
  Req<-Req1*Req2
  
  #get Meq
  yeq<- (K*d)/((bmax*(1-s))-d)-(e*aM*Req)
  tempo<-(1-ps)/(u*yeq)
  
  Meq<-0.5*(tempo+sqrt((tempo^2)+(4*tempo*KA)))
  
  # get Neq
  Neq<-((phi/Req) -(aM*Meq))/aN
  
  #get Aeq
  PMeq<-Meq/(Meq+Neq)
  alphaeq<-(Meq+Neq)*(1-f+(f*PMeq))
  Aeq<-yeq*alphaeq
  
  #return(list(Aeq=Aeq,
  #            Req=Req,
  #            Meq=Meq,
  #            Neq=Neq))
  
  return(Neq) # to get uniroot solution
}

#===============================
# Now, call the function
d<-0.5
bmax<-0.8

s<-0.1
fmin<-(s*bmax)/(bmax-d)
#fmin
#curve(myfun, xlim=c(fmin,1), lwd=2, lty=2, ylab='Neq',xlab='f',ylim=c(-0.1,1))
#abline(h=0)
myfmax<-uniroot(myfun, interval=c(fmin,0.9)) 
fmax<-myfmax$root
cat("range of fidelity for coexistence: (fmin,fmax)=(",fmin,",",fmax,")\n")
(fmax-fmin)

#when aN=aM=0.1, s=0.1
s<-0.1
fmin<-(s*bmax)/(bmax-d)
myfmax<-uniroot(myfun, interval=c(fmin,0.9), s=s, aM=0.1, aN=0.1) #fmax=0.3146655
fmax<-myfmax$root
cat("range of fidelity for coexistence: (fmin,fmax)=(",fmin,",",fmax,")\n")
(fmax-fmin)

# when s=0, aN>aM 
s<-0
fmin<-(s*bmax)/(bmax-d)
myfmax<-uniroot(myfun, interval=c(fmin,0.9), s=s, aM=0.1, aN=0.2) #fmax=0.2103166
fmax<-myfmax$root
cat("range of fidelity for coexistence: (fmin,fmax)=(",fmin,",",fmax,")\n")
(fmax-fmin)

#======================================================================================================
# To see plot comment the uniroot return line in the myfun and uncomment the other return statement
ans<-myfun(f=fmin)
(Aeq<-ans$Aeq)
(Req<-ans$Req)
(Meq<-ans$Meq)
(Neq<-ans$Neq)

# Analytical plot of eqm soln vs fidelity
Aeqs<-c()
Reqs<-c()
Meqs<-c()
Neqs<-c()
frange<-seq(from=fmin,to=1,by=0.01)
for(f in frange){
  #cat("f=",f,"\n")
  ans<-myfun(f=f)
  Aeqs<-c(Aeqs,ans$Aeq)
  Reqs<-c(Reqs,ans$Req)
  Meqs<-c(Meqs,ans$Meq)
  Neqs<-c(Neqs,ans$Neq)
}

pdf("./ARMN_Results/analytical_MNeqm_vs_f_ps_0.3_KM_10_KN_10_phi_5.pdf",width=8,height=8)
op<-par(mar=c(6,6,2,2),pty="s")
ylm<-round(max(Meqs,Neqs[-1]))
plot(frange,Meqs,type="l",ylab="",xlab="f",ylim=c(-1,ylm),xlim=c(fmin,1),lwd=2,cex.lab=2.5,cex.axis=2)
lines(frange,Neqs,lty="dashed",lwd=2)
abline(h=0,col="red")
points(x=0.4947031,y=0,cex=1.5,pch=19) # this is the fmax

legend("topright", c(expression(hat(M)),expression(hat(N))), col = c("black", "black"),
       cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1.2,x.intersp = 0.2,
       bty = "n")

par(op)
dev.off()

pdf("./ARMN_Results/analytical_AReqm_vs_f_ps_0.3_KM_10_KN_10_phi_5.pdf",width=8,height=8)
op<-par(mar=c(6,6,2,2),pty="s")
ylm<-round(max(Reqs,Aeqs[-1]))
plot(frange,Aeqs,type="l",ylab="",xlab="f",ylim=c(-1,ylm),xlim=c(fmin,1),lwd=2,cex.lab=2.5,cex.axis=2)
lines(frange,Reqs,lty="dashed",lwd=2)# Beyond fmax~0.49, Reqm should be constant as D/aM*Meqm 
                                    # as Neqm goes to 0 then the expression for Reqm is not valid there.
abline(h=0,col="red")

legend("topright", c(expression(hat(A)),expression(hat(R))), col = c("black", "black"),
       cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1.2,x.intersp = 0.2,
       bty = "n")

par(op)
dev.off()
#==========================================================================
# Now, get max N values when Meq=0

# defining parameters
f<-0.3
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

fmin<-(s*bmax)/(bmax-d)
fmin

# get Req
Req1<-((1-f)*K*d)/(e*(aN-aM*(1-f)))
Req2<-(1/((1-f)*(bmax-d))) - (1/((bmax*(1-s))-d))
Req<-Req1*Req2

(Neq_max<-(1/aN)*(phi/Req))

#Neq_max<-((e*phi)/(aN*K*d))*((1-fmin)*(bmax-d)*(aN-(aM*(1-f))))/(f-fmin)
#Neq_max





















