# This function returns the analytical equlibrium values of M, N, A, R

# Must be:
# KM=KN
# eM=eN

get_MNAR_eqm_analytical<-function(f,ps=0.3,s=0.1,aM=0.1,aN=0.2,phi=5,getalleqmval=F){
  
  # defining parameters
  #phi<-5
  #aM<-0.1 
  #aN<-0.2
  d<-0.5
  bmax<-0.8
  #s<-0.1
  e<-0.5
  u<-0.4
  K<-10
  KA<-5
  
  # get Req
  Req1<-((1-f)*K*d)/(e*(aN-aM*(1-f)))
  Req2<-(1/((1-f)*(bmax-d))) - (1/((bmax*(1-s))-d))
  Req<-Req1*Req2
  
  #get Meq
  yeq<- (K*d)/((bmax*(1-s))-d)-(e*aM*Req) # this is A by alpha at eqm from mutualist ZNGI eqn.
  tempo<-(1-ps)/(u*yeq)
  
  Meq<-0.5*(tempo+sqrt((tempo^2)+(4*tempo*KA)))
  
  # get Neq
  Neq<-((phi/Req) -(aM*Meq))/aN
  
  #get Aeq
  PMeq<-Meq/(Meq+Neq)
  alphaeq<-(Meq+Neq)*(1-f+(f*PMeq))
  Aeq<-yeq*alphaeq
  
  if(getalleqmval==T){
    return(list(Aeq=Aeq,
                Req=Req,
                Meq=Meq,
                Neq=Neq))
  }
  
  if(getalleqmval==F){
    return(Neq) # to get uniroot solution
  }
  
}



# test case