source("./get_MNAR_eqm.R")
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
#ans<-nleqslv(x=c(0.5,0.1,0.1,0.5),fn=get_MNAR_eqm,KM=10,KN=10,f=0.3,ps=0.3)
#ans$x
#if(all(ans$x>0)){
#  eg<-get_eigenvalues(KM=10,KN=10,f=0.6,ps=0.3,eqmvals=ans$x)
#  maxeg<-max(eg)
#  maxeg
#}else{
#  cat("no co-existence with given parameters")
#}

#-----------------------------------------------------------------------------------------
# Now call the functions to see eigen value variation against fidelity with KM=KN
KM<- 10 # half saturation constant for mutualist 
KN<- 10 # half saturation constant for non-mutualist

bmax<-0.8
s<-0.1
d<-0.5
tempo<-((bmax*(1-s))-d)/(bmax-d)
fmin<-1-((KN/KM)*(tempo))

# read data file from fortran code output
fAR<-read.delim("./ARMN_Results/ARMN_dat/ARMN_fAR_ps_0.3_km_10_kn_10_phi_5.dat",sep="",header = F) # read f, A_eqm, R_eqm for ps=0.3
fMN<-read.delim("./ARMN_Results/ARMN_dat/ARMN_fMN_ps_0.3_km_10_kn_10_phi_5.dat",sep="",header=F) # read f, M_eqm, N_eqm for ps=0.3
fMNAR<-cbind(fMN,fAR[,c(2:3)])
colnames(fMNAR)<-c("f","Meq","Neq","Aeq","Req")
fMNAR<-as.data.frame(fMNAR)
fMNAR<-subset(fMNAR,f<1)
#fMNAR<-as.matrix(fMNAR)
fMNAR$maxeg<-NA

for(i in c(1:nrow(fMNAR))){
  f<-fMNAR$f[i]
  eqmvals<-fMNAR[i,c(2:5)]
  eqmvals<-as.matrix(unname(eqmvals))
  eg<-get_eigenvalues(KM=KM,KN=KN,ps=0.3,f=f,eqmvals=eqmvals)
  fMNAR$maxeg[i]<-max(eg)
}

#=============================
pdf("./ARMN_Results/max_eigenval_vs_f_with_KM_10_KN_10.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(fMNAR$f,fMNAR$maxeg,xlab="f",ylab="max(eigenvalues)",
     cex.lab=2.5,cex.axis=2,lwd=2,type="b",
     xlim=c(fMNAR$f[1],0.496),ylim=c(range(fMNAR$maxeg,0)))
abline(h=0,col="grey",lwd=2)

par(op)
dev.off()

# Now call the functions to see eigen value variation against soil P availability 

# read data file from fortran code output
psAR<-read.delim("./ARMN_Results/ARMN_dat/ARMN_psAR_f_0.3_km_10_kn_10_phi_5.dat",sep="",header = F) # read ps, A_eqm, R_eqm for f=0.3
psMN<-read.delim("./ARMN_Results/ARMN_dat/ARMN_psMN_f_0.3_km_10_kn_10_phi_5.dat",sep="",header=F) # read ps, M_eqm, N_eqm for f=0.3
psMNAR<-cbind(psMN,psAR[,c(2:3)])
colnames(psMNAR)<-c("ps","Meq","Neq","Aeq","Req")
psMNAR<-as.data.frame(psMNAR)

psMNAR$maxeg<-NA

for(i in c(1:nrow(psMNAR))){
  ps<-psMNAR$ps[i]
  eqmvals<-psMNAR[i,c(2:5)]
  eqmvals<-as.matrix(unname(eqmvals))
  eg<-get_eigenvalues(KM=KM,KN=KN,f=0.3,ps=ps,eqmvals=eqmvals)
  psMNAR$maxeg[i]<-max(eg)
}

#=============================
pdf("./ARMN_Results/max_eigenval_vs_ps_with_KM_10_KN_10.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(psMNAR$ps,psMNAR$maxeg,xlab=expression(P[s]),ylab="max(eigenvalues)",
     cex.lab=2.5,cex.axis=2,lwd=2,type="l",
     xlim=c(0,1),ylim=c(range(psMNAR$maxeg,0)))
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

# read data file from fortran code output
fAR<-read.delim("./ARMN_Results/ARMN_dat/ARMN_fAR_ps_0.3_km_10_kn_6_phi_5.dat",sep="",header = F) # read f, A_eqm, R_eqm for ps=0.3
fMN<-read.delim("./ARMN_Results/ARMN_dat/ARMN_fMN_ps_0.3_km_10_kn_6_phi_5.dat",sep="",header=F) # read f, M_eqm, N_eqm for ps=0.3
fMNAR<-cbind(fMN,fAR[,c(2:3)])
colnames(fMNAR)<-c("f","Meq","Neq","Aeq","Req")
fMNAR<-as.data.frame(fMNAR)
fMNAR<-subset(fMNAR,f<1)
#fMNAR<-as.matrix(fMNAR)
fMNAR$maxeg<-NA

for(i in c(1:nrow(fMNAR))){
  f<-fMNAR$f[i]
  eqmvals<-fMNAR[i,c(2:5)]
  eqmvals<-as.matrix(unname(eqmvals))
  eg<-get_eigenvalues(KM=KM,KN=KN,ps=0.3,f=f,eqmvals=eqmvals)
  fMNAR$maxeg[i]<-max(eg)
}

#=============================
pdf("./ARMN_Results/max_eigenval_vs_f_with_KM_10_KN_6.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(fMNAR$f,fMNAR$maxeg,xlab="f",ylab="max(eigenvalues)",
     cex.lab=2.5,cex.axis=2,lwd=2,type="b",
     xlim=c(fMNAR$f[1],0.84),ylim=c(range(fMNAR$maxeg,0)))
abline(h=0,col="grey",lwd=2)

par(op)
dev.off()

# Now call the functions to see eigen value variation against soil P availability 

# read data file from fortran code output
psAR<-read.delim("./ARMN_Results/ARMN_dat/ARMN_psAR_f_0.6_km_10_kn_6_phi_5.dat",sep="",header = F) # read ps, A_eqm, R_eqm for f=0.3
psMN<-read.delim("./ARMN_Results/ARMN_dat/ARMN_psMN_f_0.6_km_10_kn_6_phi_5.dat",sep="",header=F) # read ps, M_eqm, N_eqm for f=0.3
psMNAR<-cbind(psMN,psAR[,c(2:3)])
colnames(psMNAR)<-c("ps","Meq","Neq","Aeq","Req")
psMNAR<-as.data.frame(psMNAR)

psMNAR$maxeg<-NA

for(i in c(1:nrow(psMNAR))){
  ps<-psMNAR$ps[i]
  eqmvals<-psMNAR[i,c(2:5)]
  eqmvals<-as.matrix(unname(eqmvals))
  eg<-get_eigenvalues(KM=KM,KN=KN,f=0.6,ps=ps,eqmvals=eqmvals)
  psMNAR$maxeg[i]<-max(eg)
}

#=============================
pdf("./ARMN_Results/max_eigenval_vs_ps_with_KM_10_KN_6.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(psMNAR$ps,psMNAR$maxeg,xlab=expression(P[s]),ylab="max(eigenvalues)",
     cex.lab=2.5,cex.axis=2,lwd=2,type="l",
     xlim=c(0,1),ylim=c(range(psMNAR$maxeg,0)))
abline(h=0,col="grey",lwd=2)

par(op)
dev.off()

