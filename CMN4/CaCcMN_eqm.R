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
#f<-0.3 # fidelity
phi <- 5 # constant resource value for construction carbon
g<- 0.2 # Rate at which construction carbon is allocated to both symbionts
ps <- 0.3 # P-availability in the soil
u<-0.4 # Phosphorous uptake per unit of preferentially allocated carbon received by mutualists

#CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u) # you can check the numerical M, N equilibria value with this analytical finding

mylist<-seq(from=0,to=1,by=0.01)
df<-matrix(NA,nrow=length(mylist),ncol=5)
colnames(df)<-c("f","Caeq","Cceq","Meq","Neq")
df<-as.data.frame(df)
for(i in seq_along(mylist)){
  f<-mylist[i]
  res<-CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u)
  df$f[i]<-f
  df$Caeq[i]<-res$Caeq
  df$Cceq[i]<-res$Cceq
  df$Meq[i]<-res$Meq
  df$Neq[i]<-res$Neq
}

op<-par(mfrow=c(2,2))
plot(df$f,df$Caeq,type="b",xlab="f",ylab="Ca_eqm")
abline(h=0,col="grey")
plot(df$f,df$Cceq,type="b",xlab="f",ylab="Cc_eqm")
abline(h=0,col="grey")
plot(df$f,df$Meq,type="b",xlab="f",ylab="M_eqm")
abline(h=0,col="grey")
plot(df$f,df$Neq,type="b",xlab="f",ylab="N_eqm")
abline(h=0,col="grey")
par(op)

#----------------------------------
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
#ps <- 0.3 # P-availability in the soil
u<-0.4 # Phosphorous uptake per unit of preferentially allocated carbon received by mutualists

mylist<-seq(from=0,to=1,by=0.01)
df<-matrix(NA,nrow=length(mylist),ncol=5)
colnames(df)<-c("ps","Caeq","Cceq","Meq","Neq")
df<-as.data.frame(df)
for(i in seq_along(mylist)){
  ps<-mylist[i]
  res<-CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u)
  df$ps[i]<-ps
  df$Caeq[i]<-res$Caeq
  df$Cceq[i]<-res$Cceq
  df$Meq[i]<-res$Meq
  df$Neq[i]<-res$Neq
}

op<-par(mfrow=c(2,2))
plot(df$ps,df$Caeq,type="b",xlab="Ps",ylab="Ca_eqm")
abline(h=0,col="grey")
plot(df$ps,df$Cceq,type="b",xlab="Ps",ylab="Cc_eqm")
abline(h=0,col="grey")
plot(df$ps,df$Meq,type="b",xlab="Ps",ylab="M_eqm")
abline(h=0,col="grey")
plot(df$ps,df$Neq,type="b",xlab="Ps",ylab="N_eqm")
abline(h=0,col="grey")
par(op)





