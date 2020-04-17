
# https://rpubs.com/aaronsc32/newton-raphson-method

library(rootSolve)
#library(nleqslv)

source("get_MNAR_eqm_analytical.R")
#===============================
# Now, call the function
d<-0.5
bmax<-0.8

s<-0.1
fmin<-(s*bmax)/(bmax-d)
#fmin
#curve(get_MNAR_eqm_analytical, xlim=c(fmin,1), lwd=2, lty=2, ylab='Neq',xlab='f',ylim=c(-0.1,1))
#abline(h=0)
myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,0.9),ps=0.3,s=0.1,aM=0.1,aN=0.2,phi=5,getalleqmval=F) 
fmax<-myfmax$root #
cat("range of fidelity for coexistence: (fmin,fmax)=(",fmin,",",fmax,")\n")
(fmax-fmin)

#------------------ when aN=aM=0.1, s=0.1
s<-0.1
fmin<-(s*bmax)/(bmax-d)
myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,0.9), ps=0.3,s=0.1,aM=0.1,aN=0.1,phi=5,getalleqmval=F) #fmax=0.3146655
fmax<-myfmax$root
cat("range of fidelity for coexistence: (fmin,fmax)=(",fmin,",",fmax,")\n")
(fmax-fmin)

#----------------- when s=0, aN>aM 
s<-0
fmin<-(s*bmax)/(bmax-d)
myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,0.9), ps=0.3, s=s, aM=0.1, aN=0.2, phi=5, getalleqmval=F) #fmax=0.2103166
fmax<-myfmax$root
cat("range of fidelity for coexistence: (fmin,fmax)=(",fmin,",",fmax,")\n")
(fmax-fmin)


#========================== generate plot with variation of s, aM<aN, phi=5 and phi=1 ==================

bmax<-0.8
d<-0.5
#smax<-1-(d/bmax)
smax<-0.31
s<-0
fmins<-c() 
fmaxs_phi5_aM_0.1_aN_0.2<-c()
fmaxs_phi1_aM_0.1_aN_0.2<-c()

s_vec<-c()
while(s<=smax){
  cat("s=",s,"\n")
  s_vec<-c(s_vec,s)
  fmin<-(s*bmax)/(bmax-d) # indep. of phi, aM or aN
  fmins<-c(fmins,fmin)
  
  #---------- aM=0.1, aN=0.2 ---------------
  
  myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin+0.0001,0.999), ps=0.3, s=s, aM=0.1, aN=0.2, phi=5, getalleqmval=F) 
  fmax_phi5<-myfmax$root
  fmaxs_phi5_aM_0.1_aN_0.2<-c(fmaxs_phi5_aM_0.1_aN_0.2,fmax_phi5)
  
  myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin+0.0001,0.999), ps=0.3, s=s, aM=0.1, aN=0.2, phi=1, getalleqmval=F) 
  fmax_phi1<-myfmax$root
  fmaxs_phi1_aM_0.1_aN_0.2<-c(fmaxs_phi1_aM_0.1_aN_0.2,fmax_phi1)
  
  s<-s+0.01
}

resloc<-"./ARMN_Results/"
pdf(paste(resloc,"range_of_fidelity_vs_s_aM_0.1_aN_0.2.pdf",sep=""),width=8,height=8)

op<-par(mar=c(6,6,2,2),pty="s",family="serif")

plot(s_vec,fmins,ylim=c(0,1),type="l",lwd=2,ylab="",xlab="s",cex.lab=2.5,cex.axis=2)
lines(s_vec,fmaxs_phi1_aM_0.1_aN_0.2,col="red",lwd=2)
lines(s_vec,fmaxs_phi5_aM_0.1_aN_0.2,col="blue",lwd=2)
legend("topleft", 
       c(expression(paste("f"["max"], " , D=5")),
         expression(paste("f"["max"], " , D=1")),
         expression("f"["min"])), 
       col = c("blue", "red", "black"),
       cex = 2.5, lty = c(1, 1, 1), lwd=c(2,2,2), xpd = TRUE, 
       horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
       bty = "n") 
legend("bottomright", c(expression(paste("a"["M"],"=0.1, ")),expression(paste("a"["N"],"=0.2"))),
       cex = 2.5,
       horiz = T, x.intersp = -0.6,
       bty = "n")

par(op)
dev.off()

#---------------------------------------------------------------

#---------- aM=0.1, aN=0.1 ---------------

bmax<-0.8
d<-0.5
fmins<-c() 
fmaxs_phi5_aM_0.1_aN_0.1<-c()
fmaxs_phi1_aM_0.1_aN_0.1<-c()
smax<-0.31
s<-0
s_vec<-c()
while(s<=smax){
#  if(s>0){
    cat("s=",s,"\n")
    s_vec<-c(s_vec,s)
    
    fmin<-(s*bmax)/(bmax-d) # indep. of phi, aM or aN
    fmins<-c(fmins,fmin)
    
    myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin+0.001,0.9), ps=0.3, s=s, aM=0.1, aN=0.1, phi=5, getalleqmval=F) 
    fmax_phi5<-myfmax$root
    fmaxs_phi5_aM_0.1_aN_0.1<-c(fmaxs_phi5_aM_0.1_aN_0.1,fmax_phi5)
    
    myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin+0.001,0.9), ps=0.3, s=s, aM=0.1, aN=0.1, phi=1, getalleqmval=F) 
    fmax_phi1<-myfmax$root
    fmaxs_phi1_aM_0.1_aN_0.1<-c(fmaxs_phi1_aM_0.1_aN_0.1,fmax_phi1)
#  }
  s<-s+0.01
}

resloc<-"./ARMN_Results/"
pdf(paste(resloc,"range_of_fidelity_vs_s_aM_0.1_aN_0.1.pdf",sep=""),width=8,height=8)

op<-par(mar=c(6,6,2,2),pty="s",family="serif")

plot(s_vec,fmins,ylim=c(0,1),type="l",lwd=2,ylab="",xlab="s",cex.lab=2.5,cex.axis=2)
lines(s_vec,fmaxs_phi1_aM_0.1_aN_0.1,col="red",lwd=2)
lines(s_vec,fmaxs_phi5_aM_0.1_aN_0.1,col="blue",lwd=2)
legend("topleft", 
       c(expression(paste("f"["max"], " , D=5")),
         expression(paste("f"["max"], " , D=1")),
         expression("f"["min"])), 
       col = c("blue", "red", "black"),
       cex = 2.5, lty = c(1, 1, 1), lwd=c(2,2,2), xpd = TRUE, 
       horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
       bty = "n") 
legend("bottomright", c(expression(paste("a"["M"],"=0.1,")),expression(paste("a"["N"],"=0.1"))),
       cex = 2.5,
       horiz = T, x.intersp = -0.6,
       bty = "n")


par(op)
dev.off()



#======================================================================================================

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





















