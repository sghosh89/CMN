# ALTERNATIVE MODEL
# It will solve 4 ODEs to give a dataframe with time and 4 state variables 

library(deSolve)

#----------------- Parameters ---------------------------------------
# ps  !concentration of soil phosphorous available for plant
# u      !phosphorous uptake per unit of A ??unit
# KA    !half saturation constant for A
# f     ! fidelity of plant allocation to mutualist
# phi     ! base root growth rate
# aM     ! colonization rate of new roots by M
# aN     ! colonization rate of new roots by N
# e     ! plant's efficiency to give resources to mutualist, non-mutualist
# bmax   ! max growth rate for symbionts
# s      ! cost of mutualism
# KM     ! half saturation constant for M
# d      ! death rate for symbionts
# KN     ! half saturation constant for N
#------------------------------------------------------------------------

ARMN_alternative<-function(t,states,params){
  with(as.list(c(states,params)),{
    
    #rate of change
    pm<-M/(M+N)
    cm<-(e*aM*R)+(A/((M+N)*(1-f+(f*pm))))
    cn<-(e*aN*R)+((A*(1-f))/((M+N)*(1-f+(f*pm))))
    
    dA<- 1-ps-(A*u*(M/(M+KA))*(f+((1-f)*pm)))
    dR<- phi-(aM*R*M)-(aN*R*N)
    dM<- (((bmax*(1-s)*cm)/(KM+cm))-d)*M
    dN<- (((bmax*cn)/(KN+cn))-d)*N
    
    #return the change
    list(c(dA,dR,dM,dN))
  }
  )
}


states<-c(A=0.5,R=0.5,M=0.1,N=0.1) # initial values of 4 state variables
t<-seq(from=0,to=1000,by=0.01) # times for simulation
params<-c(ps=0.3,
          u=0.4,
          KA=5,
          f=0.5, # change this f value as: 0.2, 0.3, 0.6 in three different simulations
          phi=5,
          aM=0.1,
          aN=0.2,
          e=0.5,
          bmax=0.8,
          s=0.1,
          KM=10,
          d=0.5,
          KN=10) #parameters

out2<-ode(y=states,times=t,func=ARMN_alternative,parms=params,method="rk4")
out2<-as.data.frame(out2)

# A,R vs. time plot
plot(out2$time,out2$A,xlim=c(0,1000),type="l",ylim=c(0,100),ylab="",xlab="time")
#grid()
lines(out2$time,out2$R,xlim=c(0,1000),lty=2,ylim=c(0,100))

# M,N vs. time plot
plot(out2$time,out2$M,xlim=c(0,1000),type="l",ylim=c(0,3),ylab="",xlab="time")
#grid()
lines(out2$time,out2$N,xlim=c(0,1000),lty=2,ylim=c(0,3))

#-----------------------
# compare between models
plot(out$time,out$N,xlim=c(0,1000),type="l",ylim=c(0,0.4),ylab="",xlab="time")
lines(out2$time,out2$N,xlim=c(0,1000),col="red")


pdf("./ARMN_Results/Puptake_vs_M_N_alternative.pdf",width=8,height=8)

op<-par(mar=c(2,2,2,2),mgp=c(3,0.5,0),pty="s",family="serif")
f<-0.3
u<-0.4
KA<-5
M    <- seq(from=0,to=100,by=2)
N    <- seq(from=0,to=100,by=2)
PUfun <- function(M,N){(M/(M+KA))*u*(f+((1-f)*(M/(M+N))))}
PU    <- outer(M,N, FUN="PUfun")

persp(M,N,PU,theta = -50, phi = 25,col = "grey",xlab="Mutualist (M)",
      ylab="Non-mutualist (N)",
      zlab="P-uptake via AMF (F)",ticktype = "detailed",
      cex.lab=1.5,cex.axis=1.4)
par(op)
dev.off()



