source("get_MNAR_eqm_analytical.R")
library(rootSolve)
#---------------------------------------------------------------------------------------------------

# This plotter function genarates the plot for A vs. R for both symbionts (the ZNGIs for Eqn 8, 9) 

# Input
#     f = fidelity of plant C allocation to mutualist (M)
#     KM,KN = half saturation constant for M, N; KM=KN
#     Meq,Neq = eqm. values for M, N
#     eM, eN = energy allocation rate by the plant to M, N; eM=eN
#     aM, aN = colonization rate of M, N; aM<=aN
#     bmax, d = per capita birth and death rate; bmax>d
#     s = cost of mutualism; s>0
#     x1 = data files generated which has A,R values with time
#     xlm,ylm = x and y axes-limit in the plot
#     n=500 default values to add arrows in the trajectories
#     resloc = folder name to save the plot
#     nametag = additional info to file name
#     plot_MZNGI,plot_NZNGI = logical 

Plotter_AR<-function(f,KM=10,KN=10,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,x1,xlm,ylm,n,resloc,nametag,plot_MZNGI,plot_NZNGI){
  
  alpha<-(Meq+Neq)*(1-f+((f*Meq)/(Meq+Neq)))
  
  iM0<-(KM*d)/((bmax*(1-s))-d) #intercept for mutualist when plot A/alpha vs. R
  iN0<-(KN*d)/((bmax-d)*(1-f)) #intercept for non-mutualist when plot A/alpha vs. R
  
  sM0<- -(eM*aM) #slope for mutualist when plot A/alpha vs. R
  sN0<- -(eN*aN)/(1-f) #slope for non-mutualist when plot A/alpha vs. R
  
  
  iM<- iM0*alpha #intercept for mutualist when plot A vs. R
  iN<- iN0*alpha #intercept for non-mutualist when plot A vs. R
  
  sM<- sM0*alpha #slope for mutualist when plot A vs. R
  sN<- sN0*alpha #slope for non-mutualist when plot A vs. R
  
  #------------------------- A vs. R plot -------------------------------------------------
  #pdf(paste(resloc,nametag,"f_",f,"_KM_",KM,"_KN_",KN,"_A_vs_R.pdf",sep=""),width=8,height=8)
  
  #op<-par(mar=c(6,6,2,2),pty="s")
  
  #plot(NA,xlim=xlm,ylim=ylm,
  #     xlab="x",ylab="y",
       #xlab="R",ylab="A",
  #     cex.lab=2.5,cex.axis=2)
  
#  if(plot_MZNGI==T & plot_NZNGI==T){
#    abline(a=iM, b=sM,col="red",lwd=2)
#    abline(a =iN, b=sN,col="blue",lwd=2)
#    legend("topright", c("Mutualist ZNGI (Eqn. 8)","Non-mutualist ZNGI (Eqn. 9)"), col = c("red", "blue"),
#           cex = 2.5, lty = c(1, 1), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
#           bty = "n") 
#  }

 # if(plot_MZNGI==T){
#    abline(a=iM, b=sM,col="red",lwd=2)
#    legend("topright", c("Mutualist ZNGI (Eqn. 8)"), col = c("red"),
#           cex = 2.5, lty = c(1), lwd=c(2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
#           bty = "n") 
#  }
  
 #if(plot_NZNGI==T){
#    abline(a =iN, b=sN,col="blue",lwd=2)
#    legend("topright", c("Non-mutualist ZNGI (Eqn. 9)"), col = c("blue"),
#           cex = 2.5, lty = c(1), lwd=c(2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
#           bty = "n") 
#  }
  
  
  #abline(h=0,col="dimgrey",lty="dotted")
 # abline(v=0,col="dimgrey",lty="dotted")
  #grid()
  
#  colnames(x1)<-c("time","A","R")
 # lines(x1$R,x1$A,col="black",lty="dashed",lwd=2)
 
  
  #arrows(x1$R[which(1:nrow(x1) %% n == 0)-0.1], x1$A[which(1:nrow(x1) %% n == 0)-0.5], 
  #       x1$R[1:nrow(x1) %% n == 0], x1$A[1:nrow(x1) %% n == 0] - 0.01, angle=40, 
  #       length=0.1, col="black",lwd=2)

  #par(op)
  #dev.off()
  
  #------------------------------ A/alpha vs. R plot -----------------------------------------------
  
  pdf(paste(resloc,nametag,"f_",f,"_KM_",KM,"_KN_",KN,"_A_by_alpha_vs_R.pdf",sep=""),width=8,height=8)
  
  op<-par(mar=c(6,7,2,2),pty="s",mgp=c(3.5,1,0),family="serif")
  
  plot(NA,xlim=xlm,ylim=ylm/alpha,
       xlab="x",ylab="y",
       #xlab=expression(hat(R)),ylab=bquote(I == hat(A)/hat(alpha)),
       cex.lab=2.5,cex.axis=2)
 
  
  if(plot_MZNGI==T && plot_NZNGI==T){
    abline(a=iM0, b=sM0,col="red",lwd=2)
    abline(a =iN0, b=sN0,col="blue",lwd=2)
    legend("top", c("Eq. 8, from the \n mutualist model equation.",
                         "Eq. 9, from the \n non-mutualist model equation."), 
           col = c("red", "blue"),seg.len = c(0.8,0.8),
           cex = 2.5, lty = c(1, 1), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),
           y.intersp = 2,x.intersp = 0.1,
           bty = "n") 
  }
  
  if(plot_MZNGI==T && plot_NZNGI==F){
    abline(a=iM0, b=sM0,col="red",lwd=2)
    legend("topright", c("Mutualist ZNGI (Eqn. 8)"), col = c("red"),
           cex = 2, lty = c(1), lwd=c(2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
           bty = "n") 
  }
  
  if(plot_MZNGI==F && plot_NZNGI==T){
    abline(a =iN0, b=sN0,col="blue",lwd=2)
    legend("topright", c("Non-mutualist ZNGI (Eqn. 9)"), col = c("blue"),
           cex = 2, lty = c(1), lwd=c(2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
           bty = "n") 
  }
  
  abline(h=0,col="dimgrey",lty="dotted")
  abline(v=0,col="dimgrey",lty="dotted")
  #grid()
  
  colnames(x1)<-c("time","A","R")
  x1$A_by_alpha<- x1$A/alpha
  lines(x1$R,x1$A_by_alpha,col="black",lty="dashed",lwd=2)
  
  
  #arrows(x1$R[which(1:nrow(x1) %% n == 0)-0.1], x1$A[which(1:nrow(x1) %% n == 0)-0.5], 
  #       x1$R[1:nrow(x1) %% n == 0], x1$A[1:nrow(x1) %% n == 0] - 0.01, angle=40, 
  #       length=0.1, col="black",lwd=2)
  
  par(op)
  dev.off()
  
  #-------------------------------------------------------
}
#----------------------------------------
resloc<-"./ARMN_Results/"

# ------------------ call the function for --------------- KM=KN ---------------

# ----------------- when f < fmin ----------
f<-0.2
xlm<-c(0,600)
ylm<-c(0,40000)
#ylm=c(0,16000)
x1<-read.delim("./ARMN_Results/ARMN_dat/tAR_f_0.2_ps_0.3_km_10_kn_10.dat",sep="",header = F)
x<-read.delim("./ARMN_Results/ARMN_dat/tMN_f_0.2_ps_0.3_km_10_kn_10.dat",sep="",header = F)
colnames(x)<-c("t","M","N")
Meq<-x$M
Meq<-tail(Meq,1)
Neq<-x$N
Neq<-tail(Neq,1)
Plotter_AR(f=f,KM=10,KN=10,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,
           x1,xlm=xlm,ylm=ylm,n=5000,resloc=resloc,nametag="phi_5_",plot_MZNGI = T, plot_NZNGI = T)

# ----------------- when fmin < f < fmax ----------
f<-0.3
x1<-read.delim("./ARMN_Results/ARMN_dat/tAR_f_0.3_ps_0.3_km_10_kn_10.dat",sep="",header = F)
xlm<-c(0,600)
ylm<-c(0,120)
x<-read.delim("./ARMN_Results/ARMN_dat/tMN_f_0.3_ps_0.3_km_10_kn_10.dat",sep="",header = F)
colnames(x)<-c("t","M","N")
Meq<-x$M
Meq<-tail(Meq,1)
Neq<-x$N
Neq<-tail(Neq,1)
Plotter_AR(f=f,KM=10,KN=10,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,
           x1,xlm,ylm,n=200,resloc,nametag="phi_5_",plot_MZNGI = T, plot_NZNGI = T)

# ----------------- when f > fmax ----------
f<-0.6
x1<-read.delim("./ARMN_Results/ARMN_dat/tAR_f_0.6_ps_0.3_km_10_kn_10.dat",sep="",header = F)
xlm<-c(0,600)
ylm<-c(0,48)
x<-read.delim("./ARMN_Results/ARMN_dat/tMN_f_0.6_ps_0.3_km_10_kn_10.dat",sep="",header = F)
colnames(x)<-c("t","M","N")
Meq<-x$M
Meq<-tail(Meq,1)
Neq<-x$N
Neq<-tail(Neq,1)
Plotter_AR(f=f,KM=10,KN=10,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,
           x1,xlm,ylm,n=200,resloc,nametag="phi_5_",plot_MZNGI = T, plot_NZNGI = T)


#======================================================================================================

# ==== plotter function to plot variables (A,R or M,N) against time from numerical solution of 4 ODEs =========
Plotter_ARMN_vs_t<-function(x1,xlm,ylm,nametag,taglegend,resloc){
  
  pdf(paste(resloc,nametag,"_vs_t.pdf",sep=""),width=8,height=8)
  
  op<-par(mar=c(6,6,2,2),pty="s",family="serif")
  plot(x1[,1],x1[,2],xlab="time",ylab="",cex.lab=2.5,cex.axis=2,
       col="black",type="l",xlim=xlm,
       ylim=ylm,lwd=2)
  abline(h=0,col="grey")
  #title(main=bquote(C[c]^0 == .(phi0)),cex.main=2.5,line=-8,adj=0.8)
  lines(x1[,1],x1[,3],col="black",lty="dashed",lwd=2)
  legend("topright", taglegend, col = c("black", "black"),
         cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.2,
         bty = "n")
  
  par(op)
  dev.off()
}

#----------------------------

resloc<-"./ARMN_Results/"

# ------------------ call the function for --------------- KM = KN ---------------
# ----- for f < fmin ------
xC<-read.delim("./ARMN_Results/ARMN_dat/tAR_f_0.2_ps_0.3_km_10_kn_10.dat",sep="",header = F)
xS<-read.delim("./ARMN_Results/ARMN_dat/tMN_f_0.2_ps_0.3_km_10_kn_10.dat",sep="",header = F)

Plotter_ARMN_vs_t(x1=xC,xlm=c(0,1000),ylm=c(0,100),nametag="phi_5_f_0.2_KM_10_KN_10_AR",
                  taglegend=c("A","R"), resloc)
Plotter_ARMN_vs_t(x1=xS,xlm=c(0,1000),ylm=c(0,3),nametag="phi_5_f_0.2_KM_10_KN_10_MN",
                  taglegend=c("M","N"),resloc)

#----------for fmin < f < fmax ---------

xC<-read.delim("./ARMN_Results/ARMN_dat/tAR_f_0.3_ps_0.3_km_10_kn_10.dat",sep="",header = F)
xS<-read.delim("./ARMN_Results/ARMN_dat/tMN_f_0.3_ps_0.3_km_10_kn_10.dat",sep="",header = F)

Plotter_ARMN_vs_t(x1=xC,xlm=c(0,1000),ylm=c(0,100),nametag="phi_5_f_0.3_KM_10_KN_10_AR",
                  taglegend=c("A","R"), resloc)
Plotter_ARMN_vs_t(x1=xS,xlm=c(0,1000),ylm=c(0,3),nametag="phi_5_f_0.3_KM_10_KN_10_MN",
                  taglegend=c("M","N"),resloc)

#----------for f > fmax ---------

xC<-read.delim("./ARMN_Results/ARMN_dat/tAR_f_0.6_ps_0.3_km_10_kn_10.dat",sep="",header = F)
xS<-read.delim("./ARMN_Results/ARMN_dat/tMN_f_0.6_ps_0.3_km_10_kn_10.dat",sep="",header = F)

Plotter_ARMN_vs_t(x1=xC,xlm=c(0,1000),ylm=c(0,100),nametag="phi_5_f_0.6_KM_10_KN_10_AR",
                  taglegend=c("A","R"), resloc)
Plotter_ARMN_vs_t(x1=xS,xlm=c(0,1000),ylm=c(0,3),nametag="phi_5_f_0.6_KM_10_KN_10_MN",
                  taglegend=c("M","N"),resloc)


#==============================================================================================

################################################################################################

########## PLOT FROM ANALYTICAL RESULTS #########################

################################################################################################

# Now get analytical plot for eqm values vs. Ps (soil phosphorous)

Aeqs<-c()
Reqs<-c()
Meqs<-c()
Neqs<-c()

f<-0.3
aM<-0.1
aN<-0.2
phi<-5

ps_range<-seq(from=0,to=1,by=0.01)

for(ps in ps_range){
  #cat("ps=",ps,"\n")
  ans<-get_MNAR_eqm_analytical(f=f,ps=ps,s=0.1,aM=0.1,aN=0.2,phi=5,getalleqmval=T)
  Aeqs<-c(Aeqs,ans$Aeq)
  Reqs<-c(Reqs,ans$Req)
  Meqs<-c(Meqs,ans$Meq)
  Neqs<-c(Neqs,ans$Neq)
}

pdf("./ARMN_Results/analytical_MNeqm_vs_ps_f_0.3_KM_10_KN_10_phi_5.pdf",width=8,height=8)
op<-par(mar=c(6,6,2,2),pty="s",family="serif")
plot(ps_range,Meqs,type="l",ylab="",xlab=expression(P[s]),ylim=c(0,3),xlim=range(ps_range),lwd=2,cex.lab=2.5,cex.axis=2)
lines(ps_range,Neqs,lty="dashed",lwd=2)
abline(h=0,col="gray")

legend("topright", c(expression(hat(M)),expression(hat(N))), col = c("black", "black"),
       cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1.2,x.intersp = 0.2,
       bty = "n")

par(op)
dev.off()

pdf("./ARMN_Results/analytical_AReqm_vs_ps_f_0.3_KM_10_KN_10_phi_5.pdf",width=8,height=8)
op<-par(mar=c(6,6,2,2),pty="s",family="serif")
plot(ps_range,Aeqs,type="l",ylab="",xlab=expression(P[s]),ylim=c(0,60),xlim=range(ps_range),lwd=2,cex.lab=2.5,cex.axis=2)
lines(ps_range,Reqs,lty="dashed",lwd=2)
abline(h=0,col="gray")

legend("topright", c(expression(hat(A)),expression(hat(R))), col = c("black", "black"),
       cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1.2,x.intersp = 0.2,
       bty = "n")

par(op)
dev.off()

#===============================================================================================

# Analytical plot of eqm soln vs fidelity, f, ranging in between (fmin to fmax)

Aeqs<-c()
Reqs<-c()
Meqs<-c()
Neqs<-c()

s<-0.1
bmax<-0.8
d<-0.5
fmin<-(s*bmax)/(bmax-d)
ps<-0.3
aM<-0.1
aN<-0.2
phi<-5
myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,0.9), ps=ps,s=s,aM=aM,aN=aN,phi=phi,getalleqmval=F) #fmax=0.3146655
fmax<-myfmax$root


frange<-seq(from=fmin,to=fmax,by=0.01)

for(f in frange){
  #cat("f=",f,"\n")
  ans<-get_MNAR_eqm_analytical(f=f,ps=ps,s=s,aM=aM,aN=aN,phi=phi,getalleqmval=T)
  Aeqs<-c(Aeqs,ans$Aeq)
  Reqs<-c(Reqs,ans$Req)
  Meqs<-c(Meqs,ans$Meq)
  Neqs<-c(Neqs,ans$Neq)
}

pdf("./ARMN_Results/analytical_MNeqm_vs_f_ps_0.3_KM_10_KN_10_phi_5.pdf",width=8,height=8)
op<-par(mar=c(6,6,2,2),pty="s",family="serif")
ylm<-round(max(Meqs,Neqs[-1]))
plot(frange,Meqs,type="l",ylab="",xlab="f",ylim=c(-1,ylm),xlim=c(fmin,fmax),lwd=2,cex.lab=2.5,cex.axis=2)
lines(frange,Neqs,lty="dashed",lwd=2)
abline(h=0,col="gray")
#points(x=0.4947031,y=0,cex=1.5,pch=19) # this is the fmax

legend("topright", c(expression(hat(M)),expression(hat(N))), col = c("black", "black"),
       cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1.2,x.intersp = 0.2,
       bty = "n")

par(op)
dev.off()

pdf("./ARMN_Results/analytical_AReqm_vs_f_ps_0.3_KM_10_KN_10_phi_5.pdf",width=8,height=8)
op<-par(mar=c(6,6,2,2),pty="s",family="serif")
ylm<-round(max(Reqs,Aeqs[-1]))
plot(frange,Aeqs,type="l",ylab="",xlab="f",ylim=c(-1,ylm),xlim=c(fmin,fmax),lwd=2,cex.lab=2.5,cex.axis=2)
lines(frange,Reqs,lty="dashed",lwd=2)# Beyond fmax~0.49, Reqm should be constant as D/aM*Meqm 
# as Neqm goes to 0 then the expression for Reqm is not valid there.
abline(h=0,col="gray")

legend("topright", c(expression(hat(A)),expression(hat(R))), col = c("black", "black"),
       cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1.2,x.intersp = 0.2,
       bty = "n")

par(op)
dev.off()

#==================================================================================================

# Plot for PM=Meq/(Meq+Neq) vs. Ps (soil phosphorous) vs. fidelity,f 

s<-0.1
bmax<-0.8
d<-0.5
fmin<-(s*bmax)/(bmax-d)

ps_range<-seq(from=0,to=1,by=1/50)

ps_and_fmax<-data.frame(ps=ps_range,fmax=NA)

for(i in 1:length(ps_range)){
  ps<-ps_range[i]
  cat("ps=",ps,"\n")
  
  if(ps<1){
    myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,200), # ideally the interval we should look (fmin,1) 
                    ps=ps,s=s,aM=0.1,aN=0.2,phi=5,getalleqmval=F) 
    fmax<-myfmax$root
    
    if(fmax>1){
      fmax<-1 # force set fmax = 1 
      #stop("Stop code: fmax is not less that 1 for a given ps",call.=T)
    }
  }else{ # because at ps=1, Neq never be zero, so no root finding 
    fmax<-1
  }
  
  ps_and_fmax$fmax[i]<-fmax
}


# Now make a table for all possible combo of ps, f for co-existence 
# We varied f in between >= fmin to < fmax for a given ps

len<-3000 # give it a big number

ps_f_PM<-data.frame(ps=NA*numeric(len),
                    f=NA*numeric(len),
                    PM=NA*numeric(len),
                    Aeq=NA*numeric(len),
                    Req=NA*numeric(len),
                    Meq=NA*numeric(len),
                    Neq=NA*numeric(len))

k<-1
for(i in c(1:nrow(ps_and_fmax))){
  ps<-ps_and_fmax$ps[i]
  fmax<-ps_and_fmax$fmax[i] 
  
  #cat("ps=",ps,"fmax=",fmax,"\n")
  
  f_range<-seq(from=fmin,to=fmax,by=(fmax-fmin)/50)
  
  for(j in c(1:length(f_range))){
    
    f<-f_range[j]
    
    ans<-get_MNAR_eqm_analytical(f=f,ps=ps,s=s,aM=0.1,aN=0.2,phi=5,getalleqmval=T)
    PM<-ans$Meq/(ans$Meq+ans$Neq)
    
    ps_f_PM$f[k]<-f
    ps_f_PM$ps[k]<-ps
    ps_f_PM$PM[k]<-PM
    ps_f_PM$Aeq[k]<-ans$Aeq
    ps_f_PM$Req[k]<-ans$Req
    ps_f_PM$Meq[k]<-ans$Meq
    ps_f_PM$Neq[k]<-ans$Neq
    cat("k=",k,"f=",f,"ps=",ps,"PM=",PM,"\n")
    k<-k+1
  }
}

(z<-ps_f_PM[which(ps_f_PM$PM>1),]) # this table shows at f=fmax, PM should be 1 but instead it's slightly >1
# I think, this is because uniroot function just finds the root (that could be 0,>0,<0)
# but for meaningful biological variable Neq can't be negative, so we can consider 
# Neq goes to zero and PM = 1 for given ps, f combination

ind<-which(ps_f_PM$PM>1)
ps_f_PM[ind,]$PM<-1

# Now, at f=1, analytically Req  = NaN, so we need to omit those rows: these are the ps, f combo for which
# no fmax found within [fmin,1] in ps_and_fmax table.
ind<-which(is.nan(ps_f_PM$PM))
(ps_f_PM[ind,])

ps_f_PM<-na.omit(ps_f_PM) # delete unnecessary rows
dim(ps_f_PM)
range(ps_f_PM$PM)

ps_f_PM<-ps_f_PM[,c("ps","f","PM","Meq","Neq")]

write.csv(ps_f_PM,"./ARMN_Results/ps_f_PM.csv", row.names = F) # I made a contour 
                                                        # plot tthe first 3 columns of this csv files 
                                                        # using origin pro software.

# ==================================================================================================

  multi_plotter<-function(resloc,figname){
    
    if(figname=="M_by_N_eqm_vs_phi"){
      
      #------------------ analytical expression results -----------------------------------------------
      
      phi_seq<-seq(from=1,to=10,by=0.1)
      
      x1<-data.frame(phi=phi_seq,M_by_N_eqm=NA)
      
      for(i in c(1:length(phi_seq))){
        phi<-phi_seq[i]
        #cat(phi,"\n")
        ans<-get_MNAR_eqm_analytical(f=0.3,ps=0.3,s=0.1,aM=0.1,aN=0.2,phi=phi,getalleqmval=T)
        Meq<-ans$Meq
        Neq<-ans$Neq
        x1$M_by_N_eqm[i]<-Meq/Neq
      }
      
      pdf(paste(resloc,"M_by_N_eqm_vs_phi.pdf",sep=""),width=8,height=8)
      op<-par(mar=c(6,6,2,2),mgp=c(3,1,0),pty="s",family="serif")
      plot(x1[,1],x1[,2],xlab="D",ylab=c(expression(hat(M)/hat(N))),cex.lab=2.5,cex.axis=2,col="black",type="l",xlim=c(1,10),
           ylim=c(0,7),lwd=2)
      abline(h=1,col="black",lty="dotted",lwd=2)
      par(op)
      dev.off()
    
    }else if(figname=="Puptake_vs_M_N"){
      
      # A schematic plot using expression 
      # Plotting P uptake function by AMF
      # 3D figure
      
      pdf(paste(resloc,"Puptake_vs_M_N.pdf",sep=""),width=8,height=8)
     
      op<-par(mar=c(2,2,2,2),mgp=c(3,0.5,0),pty="s",family="serif")
      f<-0.3
      u<-0.4
      KA<-5
      M    <- seq(from=0,to=100,by=2)
      N    <- seq(from=0,to=100,by=2)
      PUfun <- function(M,N){(M/(M+KA))*u*((M/(M+N))/(1-f+(f*(M/(M+N)))))}
      PU    <- outer(M,N, FUN="PUfun")
      
      persp(M,N,PU,theta = -50, phi = 25,col = "grey",xlab="Mutualist (M)",
            ylab="Non-mutualist (N)",
            zlab="P-uptake via AMF (F)",ticktype = "detailed",
            cex.lab=1.5,cex.axis=1.4)
      par(op)
      dev.off()
    }else if(figname=="schematic_diagram"){
      
      # schematic diagram for nullclines(?)
      
      pdf(paste(resloc,"schematic_diagram.pdf",sep=""),width=8,height=8)
      linepos<- -34
      linepos2<- -11.5
     
      op<-par(mar=c(3,3,2,2),mgp=c(1,1,0),pty="s",family="serif")
      plot(-1,-2,xlim=c(0,0.8),ylim=c(0,1),xlab="x", ylab="y",
           xaxt="n",yaxt="n",cex.lab=3)
      abline(a=0.5,b=-0.7,lwd=2)
      abline(a=0.7,b=-2.5,lty="dashed",lwd=2)
      #text1<-bquote("-ea"[M])
      legend("top", c("Equation (8), from the \n mutualist model equation.",
                      expression("Slope = -ea"[M]*"."),
                      "Equation (9), from the \n non-mutualist model equation.",
                      expression("Slope = -ea"[N]*"/(1-f).")), 
             cex = 2, lty = c(1, NA, 2, NA), lwd=c(2,NA,2,NA), xpd = TRUE, horiz = F, inset = c(0,0,0,0),
             y.intersp = c(2,1.4,2,1.8),x.intersp = 0.2,
             bty = "n") 
     # mtext(adj=0.65,line=-15,(bquote("A"[M]^"*"~"= C"[cM]^"*")),cex=1.5)
    #  mtext(adj=0.7,line=-18,(bquote("A"[N]^"*"~"= C"[cN]^"*"~"/(1-f)")),cex=1.5)
    #  mtext(adj=0.96,line=linepos,(bquote("R"[M]^"*")),cex=1.5)
    #  mtext(adj=0.4,line=linepos,(bquote("R"[N]^"*")),cex=1.5)
    #  mtext(adj=-0.08,line=linepos2,(bquote("A"[M]^"*")),cex=1.5)
    #  mtext(adj=-0.08,line=-2.6,(bquote("A"[N]^"*")),cex=1.5)
      par(op)
      dev.off()
    }
  }
  
  #------------------------------------------
  # Now call the plotter function
  
resloc <- "./ARMN_Results/"
multi_plotter(resloc, figname = "M_by_N_eqm_vs_phi")
multi_plotter(resloc, figname = "Puptake_vs_M_N")
multi_plotter(resloc, figname = "schematic_diagram")













