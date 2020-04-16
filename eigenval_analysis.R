#source("./get_MNAR_eqm_analytical.R")
#--------------------------------- Function to get eigen values ------------------------

# NOTE : the input eqmvals should be a vector containing M,N,A,R eqm values in order as we get same output from get_MNAR_eqm fn.
get_eigenvalues<-function(KM=10,KN=10,f,ps,eqmvals,aM=0.1,aN=0.2,KA=5,d=0.5,bmax=0.8,s=0.1,phi=5,u=0.4,eM=0.5,eN=0.5){
  
  # note: KM=KN=K, eM=eN=e
  
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
  
  
  M<-eqmvals$Meq
  N<-eqmvals$Neq
  A<-eqmvals$Aeq
  R<-eqmvals$Req
  
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
  
  #xximg<-(tJ^2) - (4*dJ)
  
#  xx<-sqrt((tJ^2)-(4*dJ))

  if(tJ<0){
    print("stable")
  }else if (tJ>0){
    print("unstable")
  }else{
    print("imaginary eigen values")
  }
  
# This is valid only for 2 by 2 matrix for two eigen values  
#  lamda1<-0.5*(tJ+xx)
#  lamda2<-0.5*(tJ-xx)
  
 # eigenvalues<-c(lamda1,lamda2)
  
  max_re_part<-max(Re(eigen(J_mat,symmetric=FALSE,only.values=TRUE)$values))
  
  return(list(J_mat=J_mat,
              max_re_part=max_re_part))
  #return(J_mat)
}

# test the function 
#ans<-nleqslv(x=c(0.5,0.1,0.1,0.5),fn=get_MNAR_eqm,KM=10,KN=10,f=0.3,ps=0.3)
#ans$x
#if(all(ans$x>0)){
#  eg<-get_eigenvalues(KM=10,KN=10,f=0.6,ps=0.3,eqmvals=ans$x)$max_re_part
#  eg
#}else{
#  cat("no co-existence with given parameters")
#}

#-----------------------------------------------------------------------------------------


#######################################################################################################
#
################### PLOT for eigen val. vs Ps (soil phosphorous) #########################
#
#######################################################################################################

ps_range<-seq(from=0,to=1,by=0.01)
ps_maxeg<-data.frame(ps=ps_range,maxeg=NA)

for(i in c(1:nrow(ps_maxeg))){
  ps<-ps_maxeg$ps[i]
  eqmvals<-get_MNAR_eqm_analytical(f=0.3,ps=ps,s=0.1,aM=0.1,aN=0.2,phi=5,getalleqmval=T)
  eg<-get_eigenvalues(KM=10,KN=10,f=0.3,ps=ps,eqmvals=eqmvals)$max_re_part
  ps_maxeg$maxeg[i]<-eg
}

pdf("./ARMN_Results/max_eigenval_vs_ps.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(ps_maxeg$ps,ps_maxeg$maxeg,xlab=expression(P[s]),ylab="Max[Re(eigenvalues)]",
     cex.lab=2.5,cex.axis=2,lwd=2,type="l",
     xlim=c(0,1),ylim=c(range(ps_maxeg$maxeg,0)))
abline(h=0,col="grey",lwd=2)

par(op)
dev.off()

#=======================================================================================

#######################################################################################################
#
################### PLOT for eigen val. vs fidelity, f (range= fmin to fmax) #########################
#
#######################################################################################################
s<-0.1
bmax<-0.8
d<-0.5
fmin<-(s*bmax)/(bmax-d)
myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,0.9),ps=0.3,s=0.1,aM=0.1,aN=0.2,phi=5,getalleqmval=F) 
fmax<-myfmax$root
f_init<-fmin+0.00000001 # at f=fmin Req =0, Neq=infinity, so we consider a limiting case
f_range<-seq(from=f_init,to=fmax,by=(fmax-f_init)/30)
f_maxeg<-data.frame(f=f_range,maxeg=NA)

for(i in c(1:nrow(f_maxeg))){
  
  f<-f_maxeg$f[i]
  eqmvals<-get_MNAR_eqm_analytical(f=f,ps=0.3,s=0.1,aM=0.1,aN=0.2,phi=5,getalleqmval=T)
  eg<-get_eigenvalues(KM=10,KN=10,f=f,ps=0.3,eqmvals=eqmvals)$max_re_part
  f_maxeg$maxeg[i]<-eg
  cat("f=",f,"\n")
}

pdf("./ARMN_Results/max_eigenval_vs_f.pdf",width=8,height=8)
op<-par(mar=c(6,6.2,2,2),pty="s")

plot(f_maxeg$f,f_maxeg$maxeg,xlab="f",ylab="Max[Re(eigenvalues)]",
     cex.lab=2.5,cex.axis=2,lwd=2,type="l",
     xlim=range(f_range),ylim=c(range(f_maxeg$maxeg,0)))
abline(h=0,col="grey",lwd=2)
#abline(v=fMNAR$f[1],col="grey",lwd=2)
#abline(v=fMNAR$f[which((fMNAR$Meq/fMNAR$Neq)>100)[1]],col="grey",lwd=2,lty=2)

par(op)
dev.off()

