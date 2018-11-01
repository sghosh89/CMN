# This plotter function genarates the plot for Ca vs. Cc for both symbionts

# Input
#     f = fidelity of plant C allocation to mutualist (M)
#     x1,x2 = data files generated from CMN4_tCaCc.dat for different phi
#     axlim = x and y axes-limit in the plot
Plotter_CMN4<-function(f,x1,x2,axlim){
  # specify the other parameters
  k_m<- 10# half saturation constant for mutualist 
  k_n<- 30 # half saturation constant for non-mutualist
  d<- 0.5 # death rate of mutualist and non-mutualist
  bmax<- 0.8 #maximum growth rate of symbionts
  s<- 0.3 # cost of mutualism
  
  a1 <-((k_m*d)/((bmax*(1-s))-d)) #intercept for mutualist
  a2<-((k_n*d)/((bmax-d)*((1-f)^2))) #intercept for non-mutualist
  
  plot(NA,xlim=axlim,ylim=axlim,xlab="C_construction",ylab="C_allocation")
  #lines(x = c(0,max(a1,a2)+20),y=c(0,max(a1,a2)+20),col="green")
  abline(a=a1,b=-1,col="red")
  abline(a =a2 , b=-1/(1-f),col="blue")
  legend("topright", c("Mutualist","Non-mutualist"), col = c("red", "blue"),
         cex = 1, lty = c(1, 1), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.5,
         bty = "n") 
  lines(x=axlim,y=c(0,0),col="black",lty="dotted")
  #grid()
  
  colnames(x1)<-c("time","C_alloc","C_cons")
  colnames(x2)<-colnames(x1)
  lines(x1$C_cons,x1$C_alloc,col="green",lty="dashed")
  lines(x2$C_cons,x2$C_alloc,col="green4",lty="dashed")
  
  n = 150
  
  arrows(x1$C_cons[which(1:nrow(x1) %% n == 0)-0.1], x1$C_alloc[which(1:nrow(x1) %% n == 0)-0.5], 
         x1$C_cons[1:nrow(x1) %% n == 0], x1$C_alloc[1:nrow(x1) %% n == 0] - 0.01, angle=45,
         length=0.1, col="green")
  
  arrows(x2$C_cons[which(1:nrow(x2) %% n == 0)-0.1], x2$C_alloc[which(1:nrow(x2) %% n == 0)-0.5], 
         x2$C_cons[1:nrow(x2) %% n == 0], x2$C_alloc[1:nrow(x2) %% n == 0] - 0.01, angle=45,
         length=0.1, col="green4")
  
}
#----------------------------------------
# call the function
f<-0.2
axlim<-c(0,160)
x1<-read.delim("CMN4_tCaCc_f_0.2_phi_5.dat",sep="")
x2<-read.delim("CMN4_tCaCc_f_0.2_phi_5.dat",sep="")
Plotter_CMN4(f=f,x1=x1,x2=x2,axlim=axlim)


# call the function
f<-0.3
axlim<-c(0,160)
x1<-read.delim("CMN4_tCaCc_f_0.3_phi_5.dat",sep="")
x2<-read.delim("CMN4_tCaCc_f_0.3_phi_25.dat",sep="")
Plotter_CMN4(f=f,x1=x1,x2=x2,axlim=axlim)

# call the function
f<-0.4
axlim<-c(0,160)
x1<-read.delim("CMN4_tCaCc_f_0.4_phi_5.dat",sep="")
x2<-read.delim("CMN4_tCaCc_f_0.4_phi_25.dat",sep="")
Plotter_CMN4(f=f,x1=x1,x2=x2,axlim=axlim)
