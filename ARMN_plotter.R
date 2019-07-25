# This plotter function genarates the plot for A vs. R for both symbionts 

# Input
#     f = fidelity of plant C allocation to mutualist (M)
#     KM,KN = half saturation constant for M, N
#     Meq,Neq = eqm. values for M, N
#     x1 = data files generated which has A,R values with time
#     xlm,ylm = x and y axes-limit in the plot
#     n=500 default values to add arrows in the trajectories
#     resloc = folder name to save the plot
#     nametag = additional info to file name

Plotter_AR<-function(f,KM,KN,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,x1,xlm,ylm,n,resloc,nametag){
  
  alpha<-(Meq+Neq)*(1-f+((f*Meq)/(Meq+Neq)))
  
  iM<-(KM*d*alpha)/((bmax*(1-s))-d) #intercept for mutualist
  iN<-(KN*d*alpha)/((bmax-d)*(1-f)) #intercept for non-mutualist
  
  sM<- -(eM*aM*alpha) #slope for mutualist
  sN<- -(eN*aN*alpha)/(1-f) #slope for non-mutualist
  
  pdf(paste(resloc,nametag,"f_",f,"_KM_",KM,"_KN_",KN,"_A_vs_R.pdf",sep=""),width=8,height=8)
  
  op<-par(mar=c(6,6,2,2),pty="s")
  
  plot(NA,xlim=xlm,ylim=ylm,xlab="R",ylab="A",cex.lab=2.5,cex.axis=2)
  
  abline(a=iM, b=sM,col="red",lwd=2)
  abline(a =iN, b=sN,col="blue",lwd=2)
  legend("topright", c("Mutualist","Non-mutualist"), col = c("red", "blue"),
         cex = 2.5, lty = c(1, 1), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
         bty = "n") 
  abline(h=0,col="dimgrey",lty="dotted")
  abline(v=0,col="dimgrey",lty="dotted")
  #grid()
  
  colnames(x1)<-c("time","A","R")
  lines(x1$R,x1$A,col="black",lty="dashed",lwd=2)
 
  
  #arrows(x1$R[which(1:nrow(x1) %% n == 0)-0.1], x1$A[which(1:nrow(x1) %% n == 0)-0.5], 
  #       x1$R[1:nrow(x1) %% n == 0], x1$A[1:nrow(x1) %% n == 0] - 0.01, angle=40, 
  #       length=0.1, col="black",lwd=2)

  par(op)
  dev.off()
}
#----------------------------------------
resloc<-"./ARMN_Results/"

# ------------------ call the function for --------------- KM=KN ---------------

# ----------------- when f < fmin ----------
f<-0.2
xlm<-c(0,600)
ylm<-c(0,1000)
x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.2_km_10_kn_10_phi_5.dat",sep="")
x<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.2_km_10_kn_10_phi_5.dat",sep="")
colnames(x)<-c("t","M","N")
Meq<-x$M
Meq<-tail(Meq,1)
Neq<-x$N
Neq<-tail(Neq,1)
Plotter_AR(f=f,KM=10,KN=10,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,x1,xlm,ylm,n=5000,resloc=resloc,nametag="phi_5_")

# ----------------- when fmin < f < fmax ----------
f<-0.3
x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.3_km_10_kn_10_phi_0.5.dat",sep="")
xlm<-c(0,450)
ylm<-c(0,50)
x<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.3_km_10_kn_10_phi_0.5.dat",sep="")
colnames(x)<-c("t","M","N")
Meq<-x$M
Meq<-tail(Meq,1)
Neq<-x$N
Neq<-tail(Neq,1)
Plotter_AR(f=f,KM=10,KN=10,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,x1,xlm,ylm,n=200,resloc,nametag="phi_0.5_")

x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.3_km_10_kn_10_phi_5.dat",sep="")
xlm<-c(0,450)
ylm<-c(0,50)
x<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.3_km_10_kn_10_phi_5.dat",sep="")
colnames(x)<-c("t","M","N")
Meq<-x$M
Meq<-tail(Meq,1)
Neq<-x$N
Neq<-tail(Neq,1)
Plotter_AR(f=f,KM=10,KN=10,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,x1,xlm,ylm,n=200,resloc,nametag="phi_5_")

# ----------------- when f > fmax ----------
f<-0.6
x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.6_km_10_kn_10_phi_5.dat",sep="")
xlm<-c(0,450)
ylm<-c(0,50)
x<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.6_km_10_kn_10_phi_5.dat",sep="")
colnames(x)<-c("t","M","N")
Meq<-x$M
Meq<-tail(Meq,1)
Neq<-x$N
Neq<-tail(Neq,1)
Plotter_AR(f=f,KM=10,KN=10,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,x1,xlm,ylm,n=200,resloc,nametag="phi_5_")

# ------------------ call the function for --------------- KM is not equal to KN ---------------

# ----------------- when fmin < f < fmax ----------
f<-0.6
xlm<-c(0,450)
ylm<-c(0,50)
x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.6_km_10_kn_6_phi_5.dat",sep="")
x<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.6_km_10_kn_6_phi_5.dat",sep="")
colnames(x)<-c("t","M","N")
Meq<-x$M
Meq<-tail(Meq,1)
Neq<-x$N
Neq<-tail(Neq,1)
Plotter_AR(f=f,KM=10,KN=6,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,x1,xlm,ylm,n=10000,resloc=resloc,nametag="phi_5_")


# ----------------- when f > fmax ----------

f<-0.9
xlm<-c(0,450)
ylm<-c(0,100)
x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.9_km_10_kn_6_phi_5.dat",sep="")
x<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.9_km_10_kn_6_phi_5.dat",sep="")
colnames(x)<-c("t","M","N")
Meq<-x$M
Meq<-tail(Meq,1)
Neq<-x$N
Neq<-tail(Neq,1)
Plotter_AR(f=f,KM=10,KN=6,Meq,Neq,eM=0.5,eN=0.5,aM=0.1,aN=0.2,bmax=0.8,d=0.5,s=0.1,x1,xlm,ylm,n=10000,resloc=resloc,nametag="phi_5_")

#======================================================================================================

# ================== plotter function to plot variables (A,R or M,N) against time ==================
Plotter_ARMN_vs_t<-function(x1,axlim,nametag,taglegend,resloc){
  
  pdf(paste(resloc,nametag,"_vs_t.pdf",sep=""),width=8,height=8)
  
  op<-par(mar=c(6,6,2,2),pty="s")
  plot(x1[,1],x1[,2],xlab="time",ylab="",cex.lab=2.5,cex.axis=2,
       col="darkgrey",type="l",
       ylim=axlim,lwd=2)
  #title(main=bquote(C[c]^0 == .(phi0)),cex.main=2.5,line=-8,adj=0.8)
  lines(x1[,1],x1[,3],col="black",lty="dashed",lwd=2)
  legend("topright", taglegend, col = c("darkgrey", "black"),
         cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.2,
         bty = "n")
  
  par(op)
  dev.off()
}

#----------------------------

resloc<-"./ARMN_Results/"

# ------------------ call the function for --------------- KM = KN ---------------
# ----- for f < fmin ------
xC<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.2_km_10_kn_10_phi_5.dat",sep="")
xS<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.2_km_10_kn_10_phi_5.dat",sep="")

Plotter_ARMN_vs_t(x1=xC,axlim=c(0,800),nametag="phi_5_f_0.2_KM_10_KN_10_AR",
                  taglegend=c("A","R"), resloc)
Plotter_ARMN_vs_t(x1=xS,axlim=c(0,50),nametag="phi_5_f_0.2_KM_10_KN_10_MN",
                  taglegend=c("M","N"),resloc)

#----------for fmin < f < fmax ---------

xC<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.3_km_10_kn_10_phi_5.dat",sep="")
xS<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.3_km_10_kn_10_phi_5.dat",sep="")

Plotter_ARMN_vs_t(x1=xC,axlim=c(0,60),nametag="phi_5_f_0.3_KM_10_KN_10_AR",
                  taglegend=c("A","R"), resloc)
Plotter_ARMN_vs_t(x1=xS,axlim=c(0,3),nametag="phi_5_f_0.3_KM_10_KN_10_MN",
                  taglegend=c("M","N"),resloc)

#----------for f > fmax ---------

xC<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.6_km_10_kn_10_phi_5.dat",sep="")
xS<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.6_km_10_kn_10_phi_5.dat",sep="")

Plotter_ARMN_vs_t(x1=xC,axlim=c(0,100),nametag="phi_5_f_0.6_KM_10_KN_10_AR",
                  taglegend=c("A","R"), resloc)
Plotter_ARMN_vs_t(x1=xS,axlim=c(0,1),nametag="phi_5_f_0.6_KM_10_KN_10_MN",
                  taglegend=c("M","N"),resloc)

# ------------------ call the function for --------------- KM not equal to KN ---------------

# ----- for fmin < f < fmax ------
xC<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.6_km_10_kn_6_phi_5.dat",sep="")
xS<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.6_km_10_kn_6_phi_5.dat",sep="")

Plotter_ARMN_vs_t(x1=xC,axlim=c(0,60),nametag="phi_5_f_0.6_KM_10_KN_6_AR",
                  taglegend=c("A","R"), resloc)
Plotter_ARMN_vs_t(x1=xS,axlim=c(0,3),nametag="phi_5_f_0.6_KM_10_KN_6_MN",
                  taglegend=c("M","N"),resloc)

#----------for f > fmax ---------

xC<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tAR_f_0.9_km_10_kn_6_phi_5.dat",sep="")
xS<-read.delim("./ARMN_Results/ARMN_dat/ARMN_tMN_f_0.9_km_10_kn_6_phi_5.dat",sep="")

Plotter_ARMN_vs_t(x1=xC,axlim=c(0,100),nametag="phi_5_f_0.9_KM_10_KN_10_AR",
                  taglegend=c("A","R"), resloc)
Plotter_ARMN_vs_t(x1=xS,axlim=c(0,1),nametag="phi_5_f_0.9_KM_10_KN_6_MN",
                  taglegend=c("M","N"),resloc)

#==============================================================================================
# plotter fn to get plots of (A,R) or (M,N) vs. ps and f

plot_ARMN_vs_ps_f<-function(x1,resloc,nametag,xlb,taglegend,axlim){
  
  pdf(paste(resloc,nametag,".pdf",sep=""),width=8,height=8)
  
  op<-par(mar=c(6,6,2,2),pty="s")
  plot(x1[,1],x1[,2],xlab=xlb,ylab="",cex.lab=2.5,cex.axis=2,
       col="darkgrey",type="l",
       ylim=axlim,lwd=2)
  
  lines(x1[,1],x1[,3],col="black",lty="dashed",lwd=2)
  legend("topright", taglegend, col = c("darkgrey", "black"),
         cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1.2,x.intersp = 0.2,
         bty = "n")
  
  par(op)
  dev.off()
}

# ------------------ call the function for --------------- KM = KN ---------------
resloc<-"./ARMN_Results/"

#-------- variation against fidelity ---------------  
x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_fAR_ps_0.3_km_10_kn_10_phi_5.dat",sep="")
nametag<-"AR_vs_f_ps_0.3_KM_10_KN_10_phi_5"
xlb<-"f"
taglegend<-c(expression(hat(A)),expression(hat(R)))
axlim<-range(c(x1[,2],x1[,3]))
plot_ARMN_vs_ps_f(x1,resloc,nametag,xlb,taglegend,axlim)

x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_fMN_ps_0.3_km_10_kn_10_phi_5.dat",sep="")
nametag<-"MN_vs_f_ps_0.3_KM_10_KN_10_phi_5"
xlb<-"f"
taglegend<-c(expression(hat(M)),expression(hat(N)))
axlim<-range(c(x1[,2],x1[,3]))
plot_ARMN_vs_ps_f(x1,resloc,nametag,xlb,taglegend,axlim)

# --------------- variation against Ps -----------------------------
x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_psAR_f_0.3_km_10_kn_10_phi_5.dat",sep="")
nametag<-"AR_vs_ps_f_0.3_KM_10_KN_10_phi_5"
xlb<-expression(P[s])
taglegend<-c(expression(hat(A)),expression(hat(R)))
axlim<-c(0,60)#range(c(0,x1[,2],x1[,3]))
plot_ARMN_vs_ps_f(x1,resloc,nametag,xlb,taglegend,axlim)

x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_psMN_f_0.3_km_10_kn_10_phi_5.dat",sep="")
nametag<-"MN_vs_ps_f_0.3_KM_10_KN_10_phi_5"
xlb<-expression(P[s])
taglegend<-c(expression(hat(M)),expression(hat(N)))
axlim<-c(0,3)
plot_ARMN_vs_ps_f(x1,resloc,nametag,xlb,taglegend,axlim)

# ------------------ call the function for --------------- KM is not equal to KN ---------------

#-------- variation against fidelity ---------------  
x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_fAR_ps_0.3_km_10_kn_6_phi_5.dat",sep="")
nametag<-"AR_vs_f_ps_0.3_KM_10_KN_6_phi_5"
xlb<-"f"
taglegend<-c(expression(hat(A)),expression(hat(R)))
axlim<-range(c(x1[,2],x1[,3]))
plot_ARMN_vs_ps_f(x1,resloc,nametag,xlb,taglegend,axlim)

x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_fMN_ps_0.3_km_10_kn_6_phi_5.dat",sep="")
nametag<-"MN_vs_f_ps_0.3_KM_10_KN_6_phi_5"
xlb<-"f"
taglegend<-c(expression(hat(M)),expression(hat(N)))
axlim<-range(c(x1[,2],x1[,3]))
plot_ARMN_vs_ps_f(x1,resloc,nametag,xlb,taglegend,axlim)

# --------------- variation against Ps -----------------------------
x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_psAR_f_0.6_km_10_kn_6_phi_5.dat",sep="")
nametag<-"AR_vs_ps_f_0.6_KM_10_KN_6_phi_5"
xlb<-expression(P[s])
taglegend<-c(expression(hat(A)),expression(hat(R)))
axlim<-c(0,60)#range(c(x1[,2],x1[,3]))
plot_ARMN_vs_ps_f(x1,resloc,nametag,xlb,taglegend,axlim)

x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_psMN_f_0.6_km_10_kn_6_phi_5.dat",sep="")
nametag<-"MN_vs_ps_f_0.6_KM_10_KN_6_phi_5"
xlb<-expression(P[s])
taglegend<-c(expression(hat(M)),expression(hat(N)))
axlim<-c(0,3)
plot_ARMN_vs_ps_f(x1,resloc,nametag,xlb,taglegend,axlim)


# ==================================================================================================
multi_plotter<-function(resloc,figname){
  if(figname=="CrSr_vs_phi"){
   
    pdf(paste(resloc,"CrSr_vs_phi.pdf",sep=""),width=8,height=8)
   
    x1<-read.delim("./ARMN_Results/ARMN_dat/ARMN_phiCrSr_f_0.3_km_10_kn_10.dat",sep="")
    op<-par(mar=c(6,6,2,2),mgp=c(4,1,0),pty="s")
    plot(x1[,1],x1[,2],xlab="D",ylab="",cex.lab=2.5,cex.axis=2,col="darkgreen",type="l",xlim=c(0,50),
         ylim=c(0,2+max(x1[,2],x1[,3])),lwd=2)
    lines(x1[,1],x1[,3],col="darkmagenta",lty="dashed",lwd=2)
    abline(h=1,col="black",lty="dotted",lwd=2)
    legend("topleft", c("R / A", "N / M"), col = c("darkgreen", "darkmagenta"),
           cex = 2.3, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1,x.intersp = 0.2,
           bty = "n") 
    par(op)
    dev.off()
  }else if(figname=="Puptake_vs_M_N"){
    # Plotting P uptake function by AMF
    # 3D figure
    
    pdf(paste(resloc,"Puptake_vs_M_N.pdf",sep=""),width=8,height=8)
   
    op<-par(mar=c(6,6,2,2),pty="s")
    f<-0.3
    u<-0.4
    KA<-5
    M    <- seq(from=0,to=100,by=2)
    N    <- seq(from=0,to=100,by=2)
    PUfun <- function(M,N){(M/(M+KA))*u*((M/(M+N))/(1-f+(f*(M/(M+N)))))}
    PU    <- outer(M,N, FUN="PUfun")
    
    persp(M,N,PU,theta = -45, phi = 25,col = "grey",xlab="Mutualist (M)",
          ylab="Non-mutualist (N)",
          zlab="P-uptake via AMF (F)",ticktype = "detailed",
          cex.lab=1.5,cex.axis=1.5)
    par(op)
    dev.off()
  }else if(figname=="schematic_diagram"){
    # schematic diagram
    
    pdf(paste(resloc,"schematic_diagram.pdf",sep=""),width=8,height=8)
    linepos<- -34
    linepos2<- -11.5
   
    op<-par(mar=c(6,6,2,2),pty="s")
    plot(-1,-2,xlim=c(0,0.8),ylim=c(0,0.8),xlab="R", ylab="A",
         xaxt="n",yaxt="n",cex.lab=1.8)
    abline(a=0.5,b=-1,lwd=2)
    abline(a=0.7,b=-2.5,lty="dashed",lwd=2)
    legend("topright", c("Mutualist","Non-mutualist"), 
           cex = 1.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1,x.intersp = 0.2,
           bty = "n") 
   # mtext(adj=0.65,line=-15,(bquote("A"[M]^"*"~"= C"[cM]^"*")),cex=1.5)
  #  mtext(adj=0.7,line=-18,(bquote("A"[N]^"*"~"= C"[cN]^"*"~"/(1-f)")),cex=1.5)
    mtext(adj=0.7,line=linepos,(bquote("R"[M]^"*")),cex=1.5)
    mtext(adj=0.4,line=linepos,(bquote("R"[N]^"*")),cex=1.5)
    mtext(adj=-0.08,line=linepos2,(bquote("A"[M]^"*")),cex=1.5)
    mtext(adj=-0.08,line=-2.6,(bquote("A"[N]^"*")),cex=1.5)
    par(op)
    dev.off()
  }
}

#------------------------------------------
# Now call the plotter function

resloc <- "./ARMN_Results/"
multi_plotter(resloc, figname = "CrSr_vs_phi")
multi_plotter(resloc, figname = "Puptake_vs_M_N")
multi_plotter(resloc, figname = "schematic_diagram")













