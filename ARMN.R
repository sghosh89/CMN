# This scripts does the same job as that of ARMN.f (when using Runge-Kutta 4th order for both FORTRAN and R code)
# But ofcourse R code will be slower than the FORTRAN code

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

ARMN<-function(t,states,params){
  with(as.list(c(states,params)),{
    
    #rate of change
    pm<-M/(M+N)
    cm<-(e*aM*R)+(A/((M+N)*(1-f+(f*pm))))
    cn<-(e*aN*R)+((A*(1-f))/((M+N)*(1-f+(f*pm))))
    
    dA<- 1-ps-(A*u*(M/(M+KA))*(pm/(1-f+(f*pm))))
    dR<- phi-(aM*R*M)-(aN*R*N)
    dM<- (((bmax*(1-s)*cm)/(KM+cm))-d)*M
    dN<- (((bmax*cn)/(KN+cn))-d)*N
    
    #return the change
    list(c(dA,dR,dM,dN))
  }
  )
}


states<-c(A=0.5,R=0.5,M=0.1,N=0.1) # initial values of 4 state variables
t<-seq(from=0,to=10000,by=0.01) # times for simulation
params<-c(ps=0.3,
          u=0.4,
          KA=5,
          f=0.6, # change this f value as: 0.2, 0.3, 0.6 in three different simulations
          phi=5,
          aM=0.1,
          aN=0.2,
          e=0.5,
          bmax=0.8,
          s=0.1,
          KM=10,
          d=0.5,
          KN=10) #parameters

out<-ode(y=states,times=t,func=ARMN,parms=params,method="rk4")
out<-as.data.frame(out)

# A,R vs. time plot
plot(out$time,out$A,xlim=c(0,1000),type="l",ylim=c(0,100),ylab="",xlab="time")
#grid()
lines(out$time,out$R,xlim=c(0,1000),lty=2,ylim=c(0,100))

# M,N vs. time plot
plot(out$time,out$M,xlim=c(0,1000),type="l",ylim=c(0,3),ylab="",xlab="time")
#grid()
lines(out$time,out$N,xlim=c(0,1000),lty=2,ylim=c(0,3))













