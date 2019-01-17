# This plotter function genarates the plot for Ca vs. Cc for both symbionts ()

# Input
#     f = fidelity of plant C allocation to mutualist (M)
#     x1,x2 = data files generated from CMN4_tCaCc.dat for different phi
#     xlm,ylm = x and y axes-limit in the plot
#     n=500 default values to add arrows in the trajectories
#     resloc = folder name to save the plot
#     figformat = default value "pdf" (other options : "eps", "jpeg")
Plotter_CMN4<-function(f,km,kn,x1,x2,xlm,ylm,n,resloc,figformat,tagon){
  # specify the other parameters
  #k_m<- 10# half saturation constant for mutualist 
  #k_n<- 10 # half saturation constant for non-mutualist
  d<- 0.5 # death rate of mutualist and non-mutualist
  bmax<- 0.8 #maximum growth rate of symbionts
  s<- 0.1 # cost of mutualism
  
  a1 <-((km*d)/((bmax*(1-s))-d)) #intercept for mutualist
  a2<-(kn*d)/((bmax-d)*(1-f)) #intercept for non-mutualist
  
  if(figformat=="eps"){
    setEPS()
    postscript(paste(resloc,"f_",f,"_km_",km,"_kn_",kn,"_Ca_vs_Cc.eps",sep=""),width=8,height=8)
  }else if(figformat=="jpeg"){
    jpeg(paste(resloc,"f_",f,"_km_",km,"_kn_",kn,"_Ca_vs_Cc.jpeg",sep=""),width=640,height=640)
  }else if(figformat=="pdf"){
    pdf(paste(resloc,"f_",f,"_km_",km,"_kn_",kn,"_Ca_vs_Cc.pdf",sep=""),width=8,height=8)
  }else{
    print("----------Error : figformat is not specified---------")
  }
  
  op<-par(mar=c(6,6,2,2))
  
  plot(NA,xlim=xlm,ylim=ylm,xlab="Construction C",ylab="Allocated C",cex.lab=2.5,cex.axis=2)
  #lines(x = c(0,max(a1,a2)+20),y=c(0,max(a1,a2)+20),col="green")
  abline(a=a1,b=-1,col="red",lwd=2)
  abline(a =a2 , b=-1/(1-f),col="blue",lwd=2)
  legend("topright", c("Mutualist","Non-mutualist"), col = c("red", "blue"),
         cex = 2.5, lty = c(1, 1), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
         bty = "n") 
  lines(x=axlim,y=c(0,0),col="dimgrey",lty="dotted")
  #grid()
  
  colnames(x1)<-c("time","C_alloc","C_cons")
  colnames(x2)<-colnames(x1)
  lines(x1$C_cons,x1$C_alloc,col="green4",lty="dashed",lwd=2)
  lines(x2$C_cons,x2$C_alloc,col="black",lty="dashed",lwd=2)
  
  arrows(x1$C_cons[which(1:nrow(x1) %% n == 0)-0.1], x1$C_alloc[which(1:nrow(x1) %% n == 0)-0.5], 
         x1$C_cons[1:nrow(x1) %% n == 0], x1$C_alloc[1:nrow(x1) %% n == 0] - 0.01, angle=45, 
         length=0.1, col="green4",lwd=2)
  
  arrows(x2$C_cons[which(1:nrow(x2) %% n == 0)-0.1], x2$C_alloc[which(1:nrow(x2) %% n == 0)-0.5], 
         x2$C_cons[1:nrow(x2) %% n == 0], x2$C_alloc[1:nrow(x2) %% n == 0] - 0.01, angle=45,
         length=0.1, col="black",lwd=2)
  
  if(tagon==T){
    legend(x=-2.5,y=ylm[2]+4, c(expression(paste("C"[c]^0, " = ", 0.5)),expression(paste("C"[c]^0, " = ", 5))),
           col = c("green4", "black"),
           cex = 2.5, lty = c(2, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),
           y.intersp = 1.3,x.intersp = 0.1,
           bty = "n") 
  }
  par(op)
  dev.off()
}
#----------------------------------------
# call the function
f<-0.2
xlm<-c(0,30)
ylm<-c(0,55)
x1<-read.delim("CMN4_tCaCc_f_0.2_phi_5.0_km_10.0_kn_10.0.dat",sep="")
x2<-read.delim("CMN4_tCaCc_f_0.2_phi_5.0_km_10.0_kn_10.0.dat",sep="")
Plotter_CMN4(f=f,km=10,kn=10,x1=x1,x2=x2,xlm=xlm,ylm=ylm,n=300,resloc="./Results/pdf_fig/",figformat = "pdf",tagon = F)
Plotter_CMN4(f=f,km=10,kn=10,x1=x1,x2=x2,xlm=xlm,ylm=ylm,n=300,resloc="./Results/jpeg_fig/",figformat = "jpeg",tagon = F)
Plotter_CMN4(f=f,km=10,kn=10,x1=x1,x2=x2,xlm=xlm,ylm=ylm,n=300,resloc="./Results/eps_fig/",figformat = "eps",tagon = F)

# call the function
f<-0.3
xlm<-c(0,30)
ylm<-c(0,55)
x1<-read.delim("CMN4_tCaCc_f_0.3_phi_0.5_km_10.0_kn_10.0.dat",sep="")
x2<-read.delim("CMN4_tCaCc_f_0.3_phi_5.0_km_10.0_kn_10.0.dat",sep="")
Plotter_CMN4(f=f,km=10,kn=10,x1=x1,x2=x2,xlm=xlm,ylm=ylm,n=10000,resloc="./Results/pdf_fig/",figformat = "pdf",tagon = T)
Plotter_CMN4(f=f,km=10,kn=10,x1=x1,x2=x2,xlm=xlm,ylm=ylm,n=10000,resloc="./Results/jpeg_fig/",figformat = "jpeg",tagon = T)
Plotter_CMN4(f=f,km=10,kn=10,x1=x1,x2=x2,xlm=xlm,ylm=ylm,n=10000,resloc="./Results/eps_fig/",figformat = "eps",tagon = T)

# call the function
f<-0.7
xlm<-c(0,30)
ylm<-c(0,55)
x1<-read.delim("CMN4_tCaCc_f_0.7_phi_0.5_km_10.0_kn_6.0.dat",sep="")
x2<-read.delim("CMN4_tCaCc_f_0.7_phi_5.0_km_10.0_kn_6.0.dat",sep="")
Plotter_CMN4(f=f,km=10,kn=6,x1=x1,x2=x2,xlm=xlm,ylm=ylm,n=10000,resloc="./Results/pdf_fig/",figformat = "pdf",tagon = T)
Plotter_CMN4(f=f,km=10,kn=6,x1=x1,x2=x2,xlm=xlm,ylm=ylm,n=10000,resloc="./Results/jpeg_fig/",figformat = "jpeg",tagon = T)
Plotter_CMN4(f=f,km=10,kn=6,x1=x1,x2=x2,xlm=xlm,ylm=ylm,n=10000,resloc="./Results/eps_fig/",figformat = "eps",tagon = T)


#--------------------------------------------------------------------
# plotter function to plot variables (Ca,Cc or M,N) against time
Plotter_CMN_vs_t<-function(x1,axlim,nametag,taglegend,resloc,figformat,phi0){
  
  if(figformat=="eps"){
    setEPS()
    postscript(paste(resloc,nametag,"_vs_t.eps",sep=""),width=8,height=8)
  }else if(figformat=="jpeg"){
    jpeg(paste(resloc,nametag,"_vs_t.jpeg",sep=""),width=640,height=640)
  }else if(figformat=="pdf"){
    pdf(paste(resloc,nametag,"_vs_t.pdf",sep=""),width=8,height=8)
  }else{
    print("----------Error : figformat is not specified---------")
  }
  
  op<-par(mar=c(6,6,2,2))
  plot(x1[,1],x1[,2],xlab="time",ylab="",cex.lab=2.5,cex.axis=2,
       col="darkgreen",type="l",
       ylim=axlim,lwd=2)
  title(main=bquote(C[c]^0 == .(phi0)),cex.main=2.5,line=-8,adj=0.8)
  lines(x1[,1],x1[,3],col="darkmagenta",lty="dashed",lwd=2)
  legend("topright", taglegend, col = c("darkgreen", "darkmagenta"),
         cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.2,
         bty = "n")
  
  par(op)
  dev.off()
}

#------------------------------
# for f=0.2
x5C<-read.delim("CMN4_tCaCc_f_0.2_phi_5.0_km_10.0_kn_10.0.dat",sep="")
x5S<-read.delim("CMN4_tMN_f_0.2_phi_5.0_km_10.0_kn_10.0.dat",sep="")

axlim<-c(-5,10000)
Plotter_CMN_vs_t(x1=x5C,axlim=axlim,nametag="f_0.2_phi_5.0_km_10.0_kn_10.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=5)
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.2_phi_5.0_km_10.0_kn_10.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=5)

Plotter_CMN_vs_t(x1=x5C,axlim=axlim,nametag="f_0.2_phi_5.0_km_10.0_kn_10.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=5)
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.2_phi_5.0_km_10.0_kn_10.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=5)

Plotter_CMN_vs_t(x1=x5C,axlim=axlim,nametag="f_0.2_phi_5.0_km_10.0_kn_10.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=5)
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.2_phi_5.0_km_10.0_kn_10.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=5)
#-----------------------------

# for f=0.3, phi=0.5
x0.5C<-read.delim("CMN4_tCaCc_f_0.3_phi_0.5_km_10.0_kn_10.0.dat",sep="") # for allocated and construction C
x0.5S<-read.delim("CMN4_tMN_f_0.3_phi_0.5_km_10.0_kn_10.0.dat",sep="") # for both symbionts

# for f=0.3, phi=5.0
x5C<-read.delim("CMN4_tCaCc_f_0.3_phi_5.0_km_10.0_kn_10.0.dat",sep="")
x5S<-read.delim("CMN4_tMN_f_0.3_phi_5.0_km_10.0_kn_10.0.dat",sep="")

axlim<-c(0,max(x0.5C[,2],x0.5C[,3],x5C[,2],x5C[,3]))
Plotter_CMN_vs_t(x1=x0.5C,axlim=axlim,nametag="f_0.3_phi_0.5_km_10.0_kn_10.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=0.5)
Plotter_CMN_vs_t(x1=x5C,axlim=axlim,nametag="f_0.3_phi_5_km_10.0_kn_10.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=5)

Plotter_CMN_vs_t(x1=x0.5C,axlim=axlim,nametag="f_0.3_phi_0.5_km_10.0_kn_10.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=0.5)
Plotter_CMN_vs_t(x1=x5C,axlim=axlim,nametag="f_0.3_phi_5_km_10.0_kn_10.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=5)

Plotter_CMN_vs_t(x1=x0.5C,axlim=axlim,nametag="f_0.3_phi_0.5_km_10.0_kn_10.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=0.5)
Plotter_CMN_vs_t(x1=x5C,axlim=axlim,nametag="f_0.3_phi_5_km_10.0_kn_10.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=5)

axlim<-c(0,max(x5S[,2],x5S[,3],x5S[,2],x5S[,3]))
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.3_phi_0.5_km_10.0_kn_10.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=0.5)
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.3_phi_5_km_10.0_kn_10.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=5)

Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.3_phi_0.5_km_10.0_kn_10.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=0.5)
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.3_phi_5_km_10.0_kn_10.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=5)

Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.3_phi_0.5_km_10.0_kn_10.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=0.5)
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.3_phi_5_km_10.0_kn_10.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=5)
#-----------------------------------------------------------------------------
# for f=0.7, phi=0.5
x0.5C<-read.delim("CMN4_tCaCc_f_0.7_phi_0.5_km_10.0_kn_6.0.dat",sep="") # for allocated and construction C
x0.5S<-read.delim("CMN4_tMN_f_0.7_phi_0.5_km_10.0_kn_6.0.dat",sep="") # for both symbionts

# for f=0.7, phi=5.0
x5C<-read.delim("CMN4_tCaCc_f_0.7_phi_5.0_km_10.0_kn_6.0.dat",sep="")
x5S<-read.delim("CMN4_tMN_f_0.7_phi_5.0_km_10.0_kn_6.0.dat",sep="")

axlim<-c(0,max(x0.5C[,2],x0.5C[,3],x5C[,2],x5C[,3]))
Plotter_CMN_vs_t(x1=x0.5C,axlim=axlim,nametag="f_0.7_phi_0.5_km_10.0_kn_6.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=0.5)
Plotter_CMN_vs_t(x1=x5C,axlim=axlim,nametag="f_0.7_phi_5_km_10.0_kn_6.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=5)

Plotter_CMN_vs_t(x1=x0.5C,axlim=axlim,nametag="f_0.7_phi_0.5_km_10.0_kn_6.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=0.5)
Plotter_CMN_vs_t(x1=x5C,axlim=axlim,nametag="f_0.7_phi_5_km_10.0_kn_6.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=5)

Plotter_CMN_vs_t(x1=x0.5C,axlim=axlim,nametag="f_0.7_phi_0.5_km_10.0_kn_6.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=0.5)
Plotter_CMN_vs_t(x1=x5C,axlim=axlim,nametag="f_0.7_phi_5_km_10.0_kn_6.0_CaCc",taglegend=c("Allocated C","Construction C"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=5)

axlim<-c(0,max(x5S[,2],x5S[,3],x5S[,2],x5S[,3]))
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.7_phi_0.5_km_10.0_kn_6.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=0.5)
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.7_phi_5_km_10.0_kn_6.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/pdf_fig/",figformat = "pdf",phi0=5)

Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.7_phi_0.5_km_10.0_kn_6.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=0.5)
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.7_phi_5_km_10.0_kn_6.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/eps_fig/",figformat = "eps",phi0=5)

Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.7_phi_0.5_km_10.0_kn_6.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=0.5)
Plotter_CMN_vs_t(x1=x5S,axlim=axlim,nametag="f_0.7_phi_5_km_10.0_kn_6.0_MN",taglegend=c("Mutualist","Non-mutualist"),
                 resloc="./Results/jpeg_fig/",figformat = "jpeg",phi0=5)

#--------------------------------------
Plotter<-function(resloc,figname,figformat){
  if(figname=="CrSr_vs_phi"){
    if(figformat=="eps"){
      setEPS()
      postscript(paste(resloc,"CrSr_vs_phi.eps",sep=""),width=8,height=8)
    }else if(figformat=="jpeg"){
      jpeg(paste(resloc,"CrSr_vs_phi.jpeg",sep=""),width=640,height=640)
    }else if(figformat=="pdf"){
      pdf(paste(resloc,"CrSr_vs_phi.pdf",sep=""),width=8,height=8)
    }else{
      print("----------Error : figformat is not specified---------")
    }
    x1<-read.delim("CMN4_phiCrSr_f_0.3.dat",sep="")
    op<-par(mar=c(6,6,2,2),mgp=c(4,1,0))
    plot(x1[,1],x1[,2],xlab=expression("C"[c]^0),ylab="",cex.lab=2.5,cex.axis=2,col="darkgreen",type="l",
         ylim=c(0,2+max(x1[,2],x1[,3])),lwd=2)
    lines(x1[,1],x1[,3],col="darkmagenta",lty="dashed",lwd=2)
    abline(h=1,col="black",lty="dotted",lwd=2)
    legend("topleft", c("Construction C / Allocated C", "Non-mutualist / Mutualist"), col = c("darkgreen", "darkmagenta"),
           cex = 2.3, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1,x.intersp = 0.2,
           bty = "n") 
    par(op)
    dev.off()
  }else if(figname=="Puptake_vs_M_N"){
    # Plotting P uptake function by AMF
    # 3D figure
    if(figformat=="eps"){
      setEPS()
      postscript(paste(resloc,"Puptake_vs_M_N.eps",sep=""),width=8,height=8)
    }else if(figformat=="jpeg"){
      jpeg(paste(resloc,"Puptake_vs_M_N.jpeg",sep=""),width=640,height=640)
    }else if(figformat=="pdf"){
      pdf(paste(resloc,"Puptake_vs_M_N.pdf",sep=""),width=8,height=8)
    }else{
      print("----------Error : figformat is not specified---------")
    }
    op<-par(mar=c(6,6,2,2))
    f<-0.3
    u<-0.4
    kc<-5.0
    M    <- seq(from=0,to=100,by=2)
    N    <- seq(from=0,to=100,by=2)
    PUfun <- function(M,N){(M/(M+kc))*u*((M/(M+N))/(1-f+(f*(M/(M+N)))))}
    PU    <- outer(M,N, FUN="PUfun")
    
    persp(M,N,PU,theta = -45, phi = 25,col = "grey",xlab="Mutualist (M)",
          ylab="Non-mutualist (N)",
          zlab="P-uptake via AMF (F)",
          cex.lab=1.5,cex.axis=1.5)
    par(op)
    dev.off()
  }else if(figname=="schematic_diagram"){
    # schematic diagram
    if(figformat=="eps"){
      setEPS()
      postscript(paste(resloc,"schematic_diagram.eps",sep=""),width=8,height=8)
      linepos<- -35
      linepos2 <- -11.5
    }else if(figformat=="jpeg"){
      jpeg(paste(resloc,"schematic_diagram.jpeg",sep=""),width=640,height=640)
      linepos<- -40
      linepos2<- -13.5
    }else if(figformat=="pdf"){
      pdf(paste(resloc,"schematic_diagram.pdf",sep=""),width=8,height=8)
      linepos<- -34
      linepos2<- -11.5
    }else{
      print("----------Error : figformat is not specified---------")
    }
    op<-par(mar=c(6,6,2,2))
    plot(-1,-2,xlim=c(0,0.8),ylim=c(0,0.8),xlab="Construction carbon", ylab="Allocated carbon",
         xaxt="n",yaxt="n",cex.lab=1.8)
    abline(a=0.5,b=-1,lwd=2)
    abline(a=0.7,b=-2.5,lty="dashed",lwd=2)
    legend("topright", c("Mutualist : slope = -1","Non-mutualist : slope = -1/(1-f)"), 
           cex = 1.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1,x.intersp = 0.2,
           bty = "n") 
    mtext(adj=0.65,line=-15,(bquote("C"[cM]^"*"~"= C"[aM]^"*")),cex=1.5)
    mtext(adj=0.7,line=-18,(bquote("C"[cN]^"*"~"= (1-f)C"[aN]^"*")),cex=1.5)
    mtext(adj=0.7,line=linepos,(bquote("C"[cM]^"*")),cex=1.5)
    mtext(adj=0.4,line=linepos,(bquote("C"[cN]^"*")),cex=1.5)
    mtext(adj=-0.08,line=linepos2,(bquote("C"[aM]^"*")),cex=1.5)
    mtext(adj=-0.08,line=-2.6,(bquote("C"[aN]^"*")),cex=1.5)
    par(op)
    dev.off()
  }
}

#------------------------------------------
# Now call the plotter function
Plotter(resloc = "./Results/pdf_fig/", figname = "CrSr_vs_phi", figformat = "pdf")
Plotter(resloc = "./Results/eps_fig/", figname = "CrSr_vs_phi", figformat = "eps")
Plotter(resloc = "./Results/jpeg_fig/", figname = "CrSr_vs_phi", figformat = "jpeg")

Plotter(resloc = "./Results/pdf_fig/", figname = "Puptake_vs_M_N", figformat = "pdf")
Plotter(resloc = "./Results/eps_fig/", figname = "Puptake_vs_M_N", figformat = "eps")
Plotter(resloc = "./Results/jpeg_fig/", figname = "Puptake_vs_M_N", figformat = "jpeg")

Plotter(resloc = "./Results/pdf_fig/", figname = "schematic_diagram", figformat = "pdf")
Plotter(resloc = "./Results/eps_fig/", figname = "schematic_diagram", figformat = "eps")
Plotter(resloc = "./Results/jpeg_fig/", figname = "schematic_diagram", figformat = "jpeg")


#test lines













