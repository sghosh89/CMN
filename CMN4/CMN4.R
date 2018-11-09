# This plotter function genarates the plot for Ca vs. Cc for both symbionts

# Input
#     f = fidelity of plant C allocation to mutualist (M)
#     x1,x2 = data files generated from CMN4_tCaCc.dat for different phi
#     axlim = x and y axes-limit in the plot
Plotter_CMN4<-function(f,x1,x2,axlim,n=500,resloc){
  # specify the other parameters
  k_m<- 10# half saturation constant for mutualist 
  k_n<- 30 # half saturation constant for non-mutualist
  d<- 0.5 # death rate of mutualist and non-mutualist
  bmax<- 0.8 #maximum growth rate of symbionts
  s<- 0.3 # cost of mutualism
  
  a1 <-((k_m*d)/((bmax*(1-s))-d)) #intercept for mutualist
  a2<-((k_n*d)/((bmax-d)*((1-f)^2))) #intercept for non-mutualist
  
  pdf(paste(resloc,"f_",f,"_Ca_vs_Cc.pdf",sep=""),width=8,height=8)
  op<-par(mar=c(5,5,2,2))
  
  plot(NA,xlim=axlim,ylim=axlim,xlab="Construction C",ylab="Allocated C",cex.lab=2,cex.axis=2)
  #lines(x = c(0,max(a1,a2)+20),y=c(0,max(a1,a2)+20),col="green")
  abline(a=a1,b=-1,col="red")
  abline(a =a2 , b=-1/(1-f),col="blue")
  legend("topright", c("Mutualist","Non-mutualist"), col = c("red", "blue"),
         cex = 1.5, lty = c(1, 1), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,
         bty = "n") 
  lines(x=axlim,y=c(0,0),col="dimgrey",lty="dotted")
  #grid()
  
  colnames(x1)<-c("time","C_alloc","C_cons")
  colnames(x2)<-colnames(x1)
  lines(x1$C_cons,x1$C_alloc,col="green4",lty="dashed")
  lines(x2$C_cons,x2$C_alloc,col="black",lty="dashed")
  
  arrows(x1$C_cons[which(1:nrow(x1) %% n == 0)-0.1], x1$C_alloc[which(1:nrow(x1) %% n == 0)-0.5], 
         x1$C_cons[1:nrow(x1) %% n == 0], x1$C_alloc[1:nrow(x1) %% n == 0] - 0.01, angle=40,
         length=0.1, col="green4")
  
  arrows(x2$C_cons[which(1:nrow(x2) %% n == 0)-0.1], x2$C_alloc[which(1:nrow(x2) %% n == 0)-0.5], 
         x2$C_cons[1:nrow(x2) %% n == 0], x2$C_alloc[1:nrow(x2) %% n == 0] - 0.01, angle=40,
         length=0.1, col="black")
  par(op)
  dev.off()
}
#----------------------------------------
# call the function
f<-0.2
axlim<-c(0,165)
x1<-read.delim("CMN4_tCaCc_f_0.2_phi_5.dat",sep="")
x2<-read.delim("CMN4_tCaCc_f_0.2_phi_5.dat",sep="")
resloc<-"./Results/"
Plotter_CMN4(f=f,x1=x1,x2=x2,axlim=axlim,resloc=resloc)


# call the function
f<-0.3
axlim<-c(0,165)
x1<-read.delim("CMN4_tCaCc_f_0.3_phi_5.dat",sep="")
x2<-read.delim("CMN4_tCaCc_f_0.3_phi_25.dat",sep="")
resloc<-"./Results/"
Plotter_CMN4(f=f,x1=x1,x2=x2,axlim=axlim,resloc=resloc)

# call the function
f<-0.4
axlim<-c(0,165)
x1<-read.delim("CMN4_tCaCc_f_0.4_phi_5.dat",sep="")
x2<-read.delim("CMN4_tCaCc_f_0.4_phi_25.dat",sep="")
resloc<-"./Results/"
Plotter_CMN4(f=f,x1=x1,x2=x2,axlim=axlim,resloc=resloc)

#--------------------------------------------------------------------
Plotter_CMN_vs_t<-function(x1,axlim,nametag,taglegend,resloc){
  pdf(paste(resloc,nametag,"_vs_t.pdf",sep=""),width=8,height=8)
  op<-par(mar=c(5,5,2,2))
  plot(x1[,1],x1[,2],xlab="time",ylab="",cex.lab=2,cex.axis=2,col="darkgreen",type="l",ylim=axlim)
  lines(x1[,1],x1[,3],col="darkmagenta",lty="dashed")
  legend("topright", taglegend, col = c("darkgreen", "darkmagenta"),
         cex = 1.5, lty = c(1, 2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,
         bty = "n") 
  par(op)
  dev.off()
}

#------------------------------
# for f=0.2
x1<-read.delim("CMN4_tCaCc_f_0.2_phi_5.dat",sep="")
nametag<-"f_0.2_phi_5_CaCc"
taglegend<-c("Allocated C","Construction C")
resloc<-"./Results/"
axlim<-c(-5,500)
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

x1<-read.delim("CMN4_tMN_f_0.2_phi_5.dat",sep="")
nametag<-"f_0.2_phi_5_MN"
taglegend<-c("Mutualist","Non-mutualist")
resloc<-"./Results/"
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

#-----------------------------
# for f=0.3, phi=5
x1<-read.delim("CMN4_tCaCc_f_0.3_phi_5.dat",sep="")
nametag<-"f_0.3_phi_5_CaCc"
taglegend<-c("Allocated C","Construction C")
resloc<-"./Results/"
axlim=c(0,max(x1[,2],x1[,3]))
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

x1<-read.delim("CMN4_tMN_f_0.3_phi_5.dat",sep="")
nametag<-"f_0.3_phi_5_MN"
taglegend<-c("Mutualist","Non-mutualist")
resloc<-"./Results/"
axlim=c(0,max(x1[,2],x1[,3]))
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

# for f=0.3, phi=25
x1<-read.delim("CMN4_tCaCc_f_0.3_phi_25.dat",sep="")
nametag<-"f_0.3_phi_25_CaCc"
taglegend<-c("Allocated C","Construction C")
resloc<-"./Results/"
axlim=c(0,max(x1[,2],x1[,3]))
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

x1<-read.delim("CMN4_tMN_f_0.3_phi_25.dat",sep="")
nametag<-"f_0.3_phi_25_MN"
taglegend<-c("Mutualist","Non-mutualist")
resloc<-"./Results/"
axlim=c(0,max(x1[,2],x1[,3]))
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

#-----------------------------
# for f=0.4,phi=5
x1<-read.delim("CMN4_tCaCc_f_0.4_phi_5.dat",sep="")
nametag<-"f_0.4_phi_5_CaCc"
taglegend<-c("Allocated C","Construction C")
resloc<-"./Results/"
axlim=c(0,max(x1[,2],x1[,3]))
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

x1<-read.delim("CMN4_tMN_f_0.4_phi_5.dat",sep="")
nametag<-"f_0.4_phi_5_MN"
taglegend<-c("Mutualist","Non-mutualist")
resloc<-"./Results/"
axlim=c(0,max(x1[,2],x1[,3]))
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

# for f=0.4, phi=25
x1<-read.delim("CMN4_tCaCc_f_0.4_phi_25.dat",sep="")
nametag<-"f_0.4_phi_25_CaCc"
taglegend<-c("Allocated C","Construction C")
resloc<-"./Results/"
axlim=c(0,max(x1[,2],x1[,3]))
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

x1<-read.delim("CMN4_tMN_f_0.4_phi_25.dat",sep="")
nametag<-"f_0.4_phi_25_MN"
taglegend<-c("Mutualist","Non-mutualist")
resloc<-"./Results/"
axlim=c(0,max(x1[,2],x1[,3]))
Plotter_CMN_vs_t(x1=x1,axlim=axlim,nametag=nametag,taglegend=taglegend,resloc=resloc)

#--------------------------------------

x1<-read.delim("CMN4_phiCrSr.dat",sep="")
pdf("./Results/CrSr_vs_phi.pdf",width=8,height=8)
op<-par(mar=c(5,5,2,2))
plot(x1[,1],x1[,2],xlab=expression(phi),ylab="",cex.lab=2.5,cex.axis=2,col="darkgreen",type="l",ylim=c(0,max(x1[,2],x1[,3])))
lines(x1[,1],x1[,3],col="darkmagenta",lty="dashed")
abline(h=1,col="black",lty="dotted")
legend("topleft", c("Construction C / Allocated C", "Non-mutualist / Mutualist"), col = c("darkgreen", "darkmagenta"),
       cex = 1.5, lty = c(1, 2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1,
       bty = "n") 
par(op)
dev.off()

#------------------------------------------
# Plotting P uptake function by AMF

# 3D figure
pdf("./Results/Puptake_vs_M_N.pdf",width=8,height=8)
op<-par(mar=c(5,5,2,2))
f<-0.3
u<-0.4
kc<-5.0
M    <- seq(from=1,to=3,by=0.1)
N    <- seq(from=0,to=3,by=0.1)
PUfun <- function(M,N){(M/(M+kc))*u*((M/(M+N))/(1-f+(f*(M/(M+N)))))}
PU    <- outer(M,N, FUN="PUfun")

persp(M,N,PU,theta = -45, phi = 25,col = "grey",xlab="Mutualist",ylab="Non-mutualist",zlab="P-uptake by AMF",cex.lab=1.5,cex.axis=1.5)
par(op)
dev.off()







