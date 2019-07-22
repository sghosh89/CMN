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
#ans<-nleqslv(x=c(0.5,0.1,0.1,0.5),fn=get_MN_eqm,KM=10,KN=10,f=0.6,ps=0.3)
#ans$x

#---- These values we get from solving 4 ODEs ----------
#Caeq=   42.8853715958310     
#Cceq=   11.6547464961649     
#Meq=  0.669353544566645     
#Neq=   1.81037069330730    

#--------------------------------- Function to get eigen values ------------------------

# NOTE : the input eqmvals should be a vector containing M,N,A,R eqm values in order as we get same output from get_MNAR_eqm fn.
get_eigenvalues<-function(KM,KN,f,ps,eqmvals,aM=0.1,aN=0.2,KA=5,d=0.5,bmax=0.8,s=0.1,phi=5,u=0.4,eM=0.5,eN=0.5){
  
  jA<-expression(1-ps-(A*(u*(M/(KA+M))*((M/(M+N))/(1-f+(f*(M/(M+N))))))))
  
  #CM<- ((eM*aM*R)+(A/((M+N)*(1-f+(f*(M/(M+N))))))) 
  jM<-expression((((bmax*(1-s)*((eM*aM*R)+(A/((M+N)*(1-f+(f*(M/(M+N))))))))/(KM+((eM*aM*R)+(A/((M+N)*(1-f+(f*(M/(M+N)))))))))-d)*M) 
  
  #CN<-((eN*aN*R)+((A*(1-f))/((M+N)*(1-f+(f*(M/(M+N)))))))
  jN<-expression((((bmax*((eN*aN*R)+((A*(1-f))/((M+N)*(1-f+(f*(M/(M+N))))))))/(KN+((eN*aN*R)+((A*(1-f))/((M+N)*(1-f+(f*(M/(M+N)))))))))-d)*N) 
  
  jR<-expression(phi-(aM*R*M)-(aN*R*N)) 
  
  j11<-D(jA,"A")
  j12<-D(jA,"M")
  j13<-D(jA,"N")
  j14<-D(jA,"R")
  
  j21<-D(jM,"A")
  j22<-D(jM,"M")
  j23<-D(jM,"N")
  j24<-D(jM,"R")
  
  j31<-D(jN,"A")
  j32<-D(jN,"M")
  j33<-D(jN,"N")
  j34<-D(jN,"R")
  
  j41<-D(jR,"A")
  j42<-D(jR,"M")
  j43<-D(jR,"N")
  j44<-D(jR,"R")
  
  #----------------evaluate at eqm------------
  
  
  M<-eqmvals[1]
  N<-eqmvals[2]
  A<-eqmvals[3]
  R<-eqmvals[4]
  
  j11_eqm<-eval(j11)
  j12_eqm<-eval(j12)
  j13_eqm<-eval(j13)
  j14_eqm<-eval(j14)
  
  j21_eqm<-eval(j21)
  j22_eqm<-eval(j22)
  j23_eqm<-eval(j23)
  j24_eqm<-eval(j24)
  
  j31_eqm<-eval(j31)
  j32_eqm<-eval(j32)
  j33_eqm<-eval(j33)
  j34_eqm<-eval(j34)
  
  j41_eqm<-eval(j41)
  j42_eqm<-eval(j42)
  j43_eqm<-eval(j43)
  j44_eqm<-eval(j44)
  
  
  #-----------------------------
  J_mat<-matrix(c(j11_eqm,j12_eqm,j13_eqm,j14_eqm,
                  j21_eqm,j22_eqm,j23_eqm,j24_eqm,
                  j31_eqm,j32_eqm,j33_eqm,j34_eqm,
                  j41_eqm,j42_eqm,j43_eqm,j44_eqm),nrow=4,ncol=4,byrow = T)
  
  tJ<-sum(diag(J_mat)) #trace of J_mat
  dJ<-det(J_mat) #det of J_Mat
  
  xx<-sqrt((tJ^2)-(4*dJ))
  
  if(tJ<0){
    print("stable")
  }else if (tJ>0){
    print("unstable")
  }else{
    print("imaginary eigen values")
  }
  
  if(xx>0){
    print("node")
  }else if(xx<0){
    print("spiral")
  }else{
    print("star")
  }
  
  lamda1<-0.5*(tJ+xx)
  lamda2<-0.5*(tJ-xx)
  
  eigenvalues<-c(lamda1,lamda2)
  
  return(eigenvalues)
}

# test the function 
ans<-nleqslv(x=c(0.5,0.1,0.1,0.5),fn=get_MNAR_eqm,KM=10,KN=10,f=0.3,ps=0.3)
ans$x
if(all(ans$x>0)){
  eg<-get_eigenvalues(KM=10,KN=10,f=0.6,ps=0.3,eqmvals=ans$x)
  maxeg<-max(eg)
  maxeg
}else{
  cat("no co-existence with given parameters")
}

#-----------------------------------------------------------------------------------------
# Now call the functions to see eigen value variation against fidelity with KM=KN
KM<- 10 # half saturation constant for mutualist 
KN<- 10 # half saturation constant for non-mutualist

bmax<-0.8
s<-0.1
d<-0.5
tempo<-((bmax*(1-s))-d)/(bmax-d)
fmin<-1-((KN/KM)*(tempo))

mylist<-seq(from=fmin,to=1,by=0.01)
df<-matrix(NA,nrow=length(mylist),ncol=6)
colnames(df)<-c("f","Aeq","Req","Meq","Neq","maxeg")
df<-as.data.frame(df)

df$f<-mylist

for(i in c(1:length(mylist))){
  f<-df$f[i]
  ans<-nleqslv(x=c(1,1,1,1),fn=get_MNAR_eqm,KM=KM,KN=KN,f=f,ps=0.3)
  eqmvals<-ans$x
  
  df$Meq[i]<-eqmvals[1]
  df$Neq[i]<-eqmvals[2]
  df$Aeq[i]<-eqmvals[3]
  df$Req[i]<-eqmvals[4]
  
  if(all(eqmvals>0)){
    eg<-get_eigenvalues(KM=KM,KN=KN,f=f,ps=0.3,eqmvals=eqmvals)
    maxeg<-max(eg)
    df$maxeg[i]<-maxeg
  }else{
    cat("no co-existence with given parameters \n")
  }
}

df<-na.omit(df)

pdf("./ARMN_Results/max_eigenval_vs_f_with_KM_10_KN_10.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(df$f,df$maxeg,xlab="f",ylab="max(eigenvalues)",
     cex.lab=2.5,cex.axis=2,lwd=2,type="l",ylim=c(range(df$maxeg,0)))
abline(h=0,col="grey",lwd=2)

par(op)
dev.off()

# Now call the functions to see eigen value variation against soil P availability 

mylist<-seq(from=0,to=1,by=0.01)
df<-matrix(NA,nrow=length(mylist),ncol=6)
colnames(df)<-c("ps","Aeq","Req","Meq","Neq","maxeg")
df<-as.data.frame(df)

df$ps<-mylist

for(i in c(1:length(mylist))){
  ps<-df$ps[i]
  ans<-nleqslv(x=c(1,1,1,1),fn=get_MNAR_eqm,KM=KM,KN=KN,f=0.3,ps=ps)
  eqmvals<-ans$x
  
  df$Meq[i]<-eqmvals[1]
  df$Neq[i]<-eqmvals[2]
  df$Aeq[i]<-eqmvals[3]
  df$Req[i]<-eqmvals[4]
  
  if(all(eqmvals>0)){
    eg<-get_eigenvalues(KM=KM,KN=KN,f=0.3,ps=ps,eqmvals=eqmvals)
    maxeg<-max(eg)
    df$maxeg[i]<-maxeg
  }else{
    cat("no co-existence with given parameters \n")
  }
}

df<-na.omit(df)

pdf("./ARMN_Results/max_eigenval_vs_ps_with_KM_10_KN_10.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(df$ps,df$maxeg,xlab=expression(P[s]),ylab="max(eigenvalues)",
     cex.lab=2.5,cex.axis=2,lwd=2,type="l",ylim=c(range(df$maxeg,0)))
abline(h=0,col="grey",lwd=2)

par(op)
dev.off()

#=======================================================================================

# Now call the functions to see eigen value variation against fidelity with KM is not equal to KN
KM<- 10 # half saturation constant for mutualist 
KN<- 6 # half saturation constant for non-mutualist

bmax<-0.8
s<-0.1
d<-0.5
tempo<-((bmax*(1-s))-d)/(bmax-d)
fmin<-1-((KN/KM)*(tempo))

mylist<-seq(from=fmin,to=1,by=0.01)
df<-matrix(NA,nrow=length(mylist),ncol=6)
colnames(df)<-c("f","Aeq","Req","Meq","Neq","maxeg")
df<-as.data.frame(df)

df$f<-mylist

for(i in c(1:length(mylist))){
  f<-df$f[i]
  ans<-nleqslv(x=c(1,1,1,1),fn=get_MNAR_eqm,KM=KM,KN=KN,f=f,ps=0.3)
  eqmvals<-ans$x
  
  df$Meq[i]<-eqmvals[1]
  df$Neq[i]<-eqmvals[2]
  df$Aeq[i]<-eqmvals[3]
  df$Req[i]<-eqmvals[4]
  
  if(all(eqmvals>0)){
    eg<-get_eigenvalues(KM=KM,KN=KN,f=f,ps=0.3,eqmvals=eqmvals)
    maxeg<-max(eg)
    df$maxeg[i]<-maxeg
  }else{
    cat("no co-existence with given parameters \n")
  }
}

df<-na.omit(df)

pdf("./ARMN_Results/max_eigenval_vs_f_with_KM_10_KN_6.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(df$f,df$maxeg,xlab="f",ylab="max(eigenvalues)",
     cex.lab=2.5,cex.axis=2,lwd=2,type="l",ylim=c(range(df$maxeg,0)))
abline(h=0,col="grey",lwd=2)

par(op)
dev.off()

# Now call the functions to see eigen value variation against soil P availability 

mylist<-seq(from=0,to=1,by=0.01)
df<-matrix(NA,nrow=length(mylist),ncol=6)
colnames(df)<-c("ps","Aeq","Req","Meq","Neq","maxeg")
df<-as.data.frame(df)

df$ps<-mylist

for(i in c(1:length(mylist))){
  ps<-df$ps[i]
  ans<-nleqslv(x=c(1,1,1,1),fn=get_MNAR_eqm,KM=KM,KN=KN,f=0.6,ps=ps)
  eqmvals<-ans$x
  
  df$Meq[i]<-eqmvals[1]
  df$Neq[i]<-eqmvals[2]
  df$Aeq[i]<-eqmvals[3]
  df$Req[i]<-eqmvals[4]
  
  if(all(eqmvals>0)){
    eg<-get_eigenvalues(KM=KM,KN=KN,f=0.6,ps=ps,eqmvals=eqmvals)
    maxeg<-max(eg)
    df$maxeg[i]<-maxeg
  }else{
    cat("no co-existence with given parameters \n")
  }
}

df<-na.omit(df)

pdf("./ARMN_Results/max_eigenval_vs_ps_with_KM_10_KN_6.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(df$ps,df$maxeg,xlab=expression(P[s]),ylab="max(eigenvalues)",
     cex.lab=2.5,cex.axis=2,lwd=2,type="l",ylim=c(range(df$maxeg,0)))
abline(h=0,col="grey",lwd=2)

par(op)
dev.off()


