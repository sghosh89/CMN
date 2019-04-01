#----------------------------------------------------
# This function returns the equlibrium values 
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
#------------------------------------------------------------------

get_eigenvalues<-function(km,kn,kc,d,bmax,s,f,phi,g,ps,u){
  jCa<-expression(1-ps-(Ca*(u*(M/(kc+M))*((M/(M+N))/(1-f+(f*(M/(M+N))))))))
  jM<-expression((((bmax*(1-s)*(Ca+Cc))/(km+Ca+Cc))-d)*M)
  jN<-expression((((bmax*(((1-f)*Ca)+Cc))/(kn+Cc+((1-f)*Ca)) - d)*N))
  jCc<-expression(phi-(g*(M+N)*Cc))
  
  j11<-D(jCa,"Ca")
  j12<-D(jCa,"M")
  j13<-D(jCa,"N")
  j14<-D(jCa,"Cc")
  
  j21<-D(jM,"Ca")
  j22<-D(jM,"M")
  j23<-D(jM,"N")
  j24<-D(jM,"Cc")
  
  j31<-D(jN,"Ca")
  j32<-D(jN,"M")
  j33<-D(jN,"N")
  j34<-D(jN,"Cc")
  
  j41<-D(jCc,"Ca")
  j42<-D(jCc,"M")
  j43<-D(jCc,"N")
  j44<-D(jCc,"Cc")
  
  #----------------evaluate at eqm------------
  
  eqmvals<-CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u)
  Ca<-eqmvals$Caeq
  Cc<-eqmvals$Cceq
  M<-eqmvals$Meq
  N<-eqmvals$Neq
  
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


#--------------------------------------------------------
# Now call the functions
km<- 10 # half saturation constant for mutualist 
kn<- 6 # half saturation constant for non-mutualist
kc<- 5 # half saturation constant for allocated carbon
d<- 0.5 # death rate of mutualist and non-mutualist
bmax<- 0.8 # maximum growth rate of symbionts
s<- 0.1 # cost of mutualism
f<-0.7 # fidelity
phi <- 5 # constant resource value for construction carbon
g<- 0.2 # Rate at which construction carbon is allocated to both symbionts
ps <- 0.3 # P-availability in the soil
u<- 0.4 # Phosphorous uptake per unit of preferentially allocated carbon received by mutualists

get_eigenvalues(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u)



#fn<-expression(a*(x^2)+(2*x*y))
#fn_d1<-D(fn,"x")
#fn_d1

#a<-2
#x<-4
#y<-3

#eval(fn_d1)





















