library(deSolve)

ARMN_model<-function(t,states,params){
  with(as.list(c(states,params)),{
   
    pM<- M/(M+N) #proportion of mutualist
    F_MN<- u*(M/(KA+M))*(pM/(1-f+(f*pM))) #saturating function
    CM<- (eM*aM*R) + (A/((M+N)*(1-f+(f*pM)))) #carbon consumed by M
    CN<- (eN*aN*R) + ((A*(1-f))/((M+N))*(1-f+(f*pM))) #carbon consumed by N
    
    #rate of change
    dA<- 1-ps-(A*F_MN)
    dR<- D-(aM*R*M)-(aN*R*N)
    dM<- (((bmax*(1-s)*CM)/(KM+CM))- d)*M
    dN<- (((bmax*CN)/(KN+CN))-d)*N
   
    #return the change
    list(c(dA,dR,dM,dN))
  }
  )
}

states<-c(A=0.5,R=0.5,M=0.1,N=0.1) # initial values of 4 state variables
                       # A = allocation rate of C, unit mass C/time
                       # R = length of uncolonized roots
                       # M = density of mutualist
                       # N = density of non-mutualist

t<-seq(from=0,to=800,by=0.01) # times for simulation

#parameters
params<-c(ps=0.3,  #concentration of soil phosphorous available for plant
          u=0.4,  #phosphorous uptake per unit of A ??unit
          KA=5, # half saturation constant for A
          f=0.6,  # fidelity of plant allocation to mutualist
          D=5,  # base root growth rate
          aM=0.5, # colonization rate of new roots by M
          aN=0.6, # colonization rate of new roots by N
          eM=0.5, # plant's efficiency to give resources to mutualist
          eN=0.5, # and non-mutualist eM should be equal to eN
          bmax=0.8, # max growth rate for symbionts
          s=0.1, # cost of mutualism
          KM=10, # half saturation constant for M
          d=0.5, # death rate for symbionts
          KN=6  # half saturation constant for N
)

out<-ode(y=states,times=t,func=ARMN_model,parms=params,method="rk4")
out<-as.data.frame(out)

op<-par(mfrow=c(2,2),pty="s")

plot(out$R,out$A,type="l",xlab="R",ylab="A")

ylm<-max(out$A,out$R)
plot(out$time,out$A,type="l",xlab="time",ylab="",ylim=c(0,ylm),main="black:A,red:R")
lines(out$time,out$R,col="red")

ylm<-max(out$M,out$N)
plot(out$time,out$M,type="l",xlab="time",ylab="",ylim=c(0,ylm),main="black:M,red:N")
lines(out$time,out$N,col="red")
par(op)


#==============================
















