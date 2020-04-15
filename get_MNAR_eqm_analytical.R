# This function returns the analytical equlibrium values of M, N, A, R

# Must be:
# KM=KN
# eM=eN

get_MNAR_eqm_analytical<-function(x,KM=10,KN=10,f,ps,aM=0.1,aN=0.2,KA=5,d=0.5,bmax=0.8,s=0.1,phi=5,u=0.4,eM=0.5,eN=0.5){
  
  if(KM!=KN | eM!=eN){
    stop("Error: KM, KN are not equal and/ or eM, eN are not equal",call. = T)
  }
  
  fmin<- (s*bmax)/(bmax-d)
  K<-KM<-KN
  e<-eM<-eN
  
  temp<-(aN-(aM*(1-f)))*(bmax-d)
    
  Req<- (K*d*(f-fmin))/(e*(1-fmin)*temp)
  A_by_alpha_eqm<- ((K*d*(bmax-d)*(aN-aM))+(s*K*d*bmax*aM))/((bmax*(1-s)-d)*temp)
  Meq<-((1-ps)+sqrt(((1-ps)^2)+(4*u*KA*A_by_alpha_eqm*(1-ps))))/(2*u*A_by_alpha_eqm)
  Neq<-((phi/Req)-(aM*Meq))/aN
  Aeq<-A_by_alpha_eqm*(Meq+((1-f)*Neq))
  
  return(list(Meq=Meq,
           Neq=Neq,
           Aeq=Aeq,
           Req=Req))
}

