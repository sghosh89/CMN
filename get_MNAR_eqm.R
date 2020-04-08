library(nleqslv)
#--------------------------------------------------------------------------
# This function returns the equlibrium values of M, N, A, R

get_MNAR_eqm<-function(x,KM,KN,f,ps,aM=0.1,aN=0.2,KA=5,d=0.5,bmax=0.8,s=0.1,phi=5,u=0.4,eM=0.5,eN=0.5){
  
  # x[1] = Meq
  # x[2] = Neq
  # x[3] = Aeq
  # x[4] = Req
  
  c0<-(bmax*(1-s))-d
  c1<-(KM*d*u*aM)*(x[1]^3)
  c2<-(((KM*d*u*aN)*x[2])-(u*eM*aM*phi*c0)-(c0*aM*(1-ps)))*(x[1]^2)
  c3<-((1-ps)*c0*(-(aN*x[2])-(KA*aM)))*x[1]
  
  neqn1<-c1+c2+c3-(KA*aN*(1-ps)*x[2]*c0) 
  
  c4<-(KN*d*u*aM)*(x[1]^3)
  c5<-(((KN*d*u*aN)*x[2])-(u*eN*aN*phi*(bmax-d))-((bmax-d)*aM*(1-ps)*(1-f)))*(x[1]^2)
  c6<-((1-ps)*(bmax-d)*(1-f)*(-(aN*x[2])-(KA*aM)))*x[1]
  neqn2<-c4+c5+c6-(KA*aN*(1-ps)*(1-f)*(bmax-d)*x[2])
  
  pm<-x[1]/(x[1]+x[2])
  FMN<- u*(x[1]/(x[1]+KA))*(pm/(1-f+(f*pm)))
  neqn3<-1-ps-(x[3]*FMN) 
  neqn4<-phi-(aM*x[1]*x[4])-(aN*x[2]*x[4])
  
  return(c(neqn1, # Meq
           neqn2, #Neq
           neqn3, #Aeq
           neqn4  #Req
  ))
}

# test the function 
#set.seed(101)
#ans<-nleqslv(x=c(1,1,1,1),fn=get_MNAR_eqm,KM=10,KN=10,f=fmin+0.01,ps=0.3,method="Newton")
#ans$x

#---- These values we get from solving 4 ODEs ----------
#Caeq=   42.8853715958310     
#Cceq=   11.6547464961649     
#Meq=  0.669353544566645     
#Neq=   1.81037069330730    




