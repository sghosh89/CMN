# This function returns the equlibrium value for mutualist and non-mutualist
CaCcMN_eqm<-function(km=10,kn=30,kc=5,d,bmax,s,f,phi,g,ps,u){
  
  beta1<- -(phi*f)/(g*(1-f))
  beta2<- (km/((bmax*(1-s))-d)) - (kn/((bmax-d)*((1-f)^2)))
  beta<-beta1/(d*beta2)
  
  alp1<- (km*d*u*beta)/((bmax*(1-s))-d)
  alp1<-alp1 - (phi*u/g) -(f*(1-ps)*beta)
  
  alp2<- -((1-ps)*kc*beta*f)-((1-f)*(1-ps)*(beta^2))
  alp3<- -(kc*(1-f)*(1-ps)*(beta^2))
  Meq<- (-alp2 + sqrt((alp2^2)-(4*alp1*alp3)))/(2*alp1)
  Neq<- beta-Meq
  
  Caeq<- ((Meq+(Neq*(1-f)))*(Meq+kc)*(1-ps))/(u*(Meq^2))
  Cceq<- phi/(g*(Meq+Neq))
  
  res<-list(Caeq=Caeq,Cceq=Cceq,Meq=Meq,Neq=Neq)
  
  return(res)
}

# call the function 
km<- 10 # half saturation constant for mutualist 
kn<- 30 # half saturation constant for non-mutualist
kc<- 5 # half saturation constant for allocated carbon
d<- 0.5 # death rate of mutualist and non-mutualist
bmax<- 0.8 # maximum growth rate of symbionts
s<- 0.3 # cost of mutualism
f<-0.3 # fidelity
phi <- 5 # constant resource value for construction carbon
g<- 0.2 # Rate at which construction carbon is allocated to both symbionts
ps <- 0.3 # P-availability in the soil
#u<-0.4 # Phosphorous uptake per unit of preferentially allocated carbon received by mutualists

#CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u) # you can check the numerical M, N equilibria value with this analytical finding

ulist<-seq(from=0,to=1,by=0.01)
df<-matrix(NA,nrow=length(ulist),ncol=5)
colnames(df)<-c("u","Caeq","Cceq","Meq","Neq")
df<-as.data.frame(df)
for(i in seq_along(ulist)){
  u<-ulist[i]
  res<-CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u)
  df$u[i]<-u
  df$Caeq[i]<-res$Caeq
  df$Cceq[i]<-res$Cceq
  df$Meq[i]<-res$Meq
  df$Neq[i]<-res$Neq
}

plot(df$u,df$Meq,type="b",xlab="u",ylab="equilibrium value of mutualist")

plot(df$u,(df$u*df$Caeq),type="b",col="red",xlab="u",ylab="uCa_eqm")
