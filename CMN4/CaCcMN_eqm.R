# This function returns the equlibrium value for mutualist and non-mutualist
CaCcMN_eqm<-function(km,kn,kc,d,bmax,s,f,phi,g,ps,u){
  
  beta1<- -(phi*f)/(g*(1-f))
  beta2<- (km/((bmax*(1-s))-d)) - (kn/((bmax-d)*((1-f))))
  beta<-beta1/(d*beta2)
  
  alp1<- (km*d*u*beta)/((bmax*(1-s))-d)
  alp1<-alp1 - (phi*u/g) -(f*(1-ps)*beta)
  
  alp2<- -((1-ps)*kc*beta*f)-((1-f)*(1-ps)*(beta^2))
  alp3<- -(kc*(1-f)*(1-ps)*(beta^2))
  Meq<- (-alp2 + sqrt((alp2^2)-(4*alp1*alp3)))/(2*alp1)
  Neq<- beta-Meq
  
 # if((is.finite(Meq)&(Meq<0))==T){
 #   Meq<-0
 #   Neq<- beta-Meq
 # }
  
  if((is.finite(Neq)&(Neq<0))==T){
    Neq<-0
  }
  
  
  Caeq<- ((Meq+(Neq*(1-f)))*(Meq+kc)*(1-ps))/(u*(Meq^2))
  Cceq<- phi/(g*(Meq+Neq))
  
  res<-list(Caeq=Caeq,Cceq=Cceq,Meq=Meq,Neq=Neq)
  
  return(res)
}
#---------------------------------------------------------------------------------------------
# Plot equilibrium values of Ca, Cc, M, N against fidelity (f)

# call the function 
km<- 10 # half saturation constant for mutualist 
kn<- 10 # half saturation constant for non-mutualist
kc<- 5 # half saturation constant for allocated carbon
d<- 0.5 # death rate of mutualist and non-mutualist
bmax<- 0.8 # maximum growth rate of symbionts
s<- 0.1 # cost of mutualism
#f<-0.3 # fidelity
phi <- 5 # constant resource value for construction carbon
g<- 0.2 # Rate at which construction carbon is allocated to both symbionts
ps <- 0.3 # P-availability in the soil
u<- 0.4 # Phosphorous uptake per unit of preferentially allocated carbon received by mutualists

#CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u) # you can check the numerical M, N equilibria value with this analytical finding

mylist<-seq(from=0,to=1,by=0.01)
df<-matrix(NA,nrow=length(mylist),ncol=5)
colnames(df)<-c("f","Caeq","Cceq","Meq","Neq")
df<-as.data.frame(df)
for(i in seq_along(mylist)){
  f<-mylist[i]
  cat("i=",i,"\n")
  res<-CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u)
  df$f[i]<-f
  df$Caeq[i]<-res$Caeq
  df$Cceq[i]<-res$Cceq
  df$Meq[i]<-res$Meq
  df$Neq[i]<-res$Neq
}

#df<-na.omit(df)
df<-df[which(df$Caeq>0),]
range(df$f)

#-----------------
pdf("./Results/pdf_fig/Ca_eqm_vs_f.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2))
plot(df$f,df$Caeq,type="b",xlab="f",ylab=expression(hat(C[a])),cex.lab=2.5,cex.axis=2,ylim=c(0,max(df$Caeq)))
abline(h=0,col="grey")
par(op)
dev.off()

pdf("./Results/pdf_fig/Cc_eqm_vs_f.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2))
plot(df$f,df$Cceq,type="b",xlab="f",ylab=expression(hat(C[c])),cex.lab=2.5,cex.axis=2,ylim=c(0,max(df$Cceq)))
abline(h=0,col="grey")
par(op)
dev.off()

pdf("./Results/pdf_fig/M_eqm_vs_f.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2))
plot(df$f,df$Meq,type="b",xlab="f",ylab=expression(hat(M)),cex.lab=2.5,cex.axis=2,ylim=c(0,max(df$Meq)))
abline(h=0,col="grey")
par(op)
dev.off()

pdf("./Results/pdf_fig/N_eqm_vs_f.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2))
plot(df$f,df$Neq,type="b",xlab="f",ylab=expression(hat(N)),cex.lab=2.5,cex.axis=2,ylim=c(0,max(df$Neq)))
abline(h=0,col="grey")
par(op)
dev.off()


#----------------------------------------------------------------------------------------------
# Plot equilibrium values of Ca, Cc, M, N against soil Phosphorous availability (ps)

# call the function 
km<- 10 # half saturation constant for mutualist 
kn<- 10 # half saturation constant for non-mutualist
kc<- 5 # half saturation constant for allocated carbon
d<- 0.5 # death rate of mutualist and non-mutualist
bmax<- 0.8 # maximum growth rate of symbionts
s<- 0.1 # cost of mutualism
f<-0.3 # fidelity
phi <- 5 # constant resource value for construction carbon
g<- 0.2 # Rate at which construction carbon is allocated to both symbionts
#ps <- 0.3 # P-availability in the soil
u<-0.4 # Phosphorous uptake per unit of preferentially allocated carbon received by mutualists

mylist<-seq(from=0,to=0.99,by=0.01)
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

pdf("./Results/pdf_fig/CaCc_eqm_vs_ps.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2))
ylm<-max(df$Caeq,df$Cceq,na.rm = T)
plot(df$ps,df$Caeq,type="l",xlab=expression(P[s]),ylab="",
     cex.lab=2.5,cex.axis=2,lwd=2,
     ylim=c(0,ylm+7))
lines(df$ps,df$Cceq,cex.lab=2.5,cex.axis=2,lty="dashed",lwd=2)
legend("topleft", c(expression(hat(C[a])),expression(hat(C[c]))), 
       cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1,x.intersp = 0.2,
       bty = "n")

legend("topright", expression(paste("C"[c]^0, " = ", 5)),
        cex = 2.5, horiz = F, inset = c(0,0),
              bty = "n") 
par(op)
dev.off()

pdf("./Results/pdf_fig/MN_eqm_vs_ps.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2))
ylm<-max(df$Meq,df$Neq,na.rm = T)
plot(df$ps,df$Meq,type="l",xlab=expression(P[s]),ylab="",
     cex.lab=2.5,cex.axis=2,lwd=2,
     ylim=c(0,ylm+1))
lines(df$ps,df$Neq,cex.lab=2.5,cex.axis=2,lty="dashed",lwd=2)
legend("topleft", c(expression(hat(M)),expression(hat(N))), 
       cex = 2.5, lty = c(1, 2), lwd=c(2,2), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 1.2,x.intersp = 0.2,
       bty = "n")

legend("topright", expression(paste("C"[c]^0, " = ", 5)),
       cex = 2.5, horiz = F, inset = c(0,0),
       bty = "n") 
par(op)
dev.off()







