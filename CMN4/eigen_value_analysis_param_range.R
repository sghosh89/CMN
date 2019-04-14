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
  
  #if((is.finite(Neq)&(Neq<0))==T){
  #  Neq<-0
  #}
  
  
  Caeq<- ((Meq+(Neq*(1-f)))*(Meq+kc)*(1-ps))/(u*(Meq^2))
  Cceq<- phi/(g*(Meq+Neq))
  
  res<-list(Caeq=Caeq,Cceq=Cceq,Meq=Meq,Neq=Neq)
  
  return(res)
}
#------------------------------------------------------------------

get_eigenvalues<-function(km,kn,kc,d,bmax,s,f,phi,g,ps,u,eqmvals,printstatus){
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
  
  if(printstatus==T){
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
  }
  
  
  lamda1<-0.5*(tJ+xx)
  lamda2<-0.5*(tJ-xx)
  
  eigenvalues<-c(lamda1,lamda2)
  
  return(eigenvalues)
}

#-----------------------------------------------------------------------------------------
# Function to judge if the eqm is stable or not with variation of f and ps

eigen_f_ps_vary<-function(fvary=fvary,psvary=psvary,km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,phi=phi,g=g,u=u,printstatus){
  
  ca_eqm_mat<-matrix(NA,nrow=length(fvary),ncol=length(psvary))
  cc_eqm_mat<-ca_eqm_mat
  m_eqm_mat<-ca_eqm_mat
  n_eqm_mat<-ca_eqm_mat
  coex_mat<-ca_eqm_mat
  
  for(i_f in c(1:length(fvary))){
    for(i_ps in c(1:length(psvary))){
      f<-fvary[i_f]
      ps<-psvary[i_ps]
      eqms<-CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u)
      ca_eqm_mat[i_f,i_ps]<-eqms$Caeq
      cc_eqm_mat[i_f,i_ps]<-eqms$Cceq
      m_eqm_mat[i_f,i_ps]<-eqms$Meq
      n_eqm_mat[i_f,i_ps]<-eqms$Neq
      
      can_by_cam_star<-(kn*((bmax*(1-s))-d))/(km*(bmax-d)*(1-f))
      if(can_by_cam_star>1 & can_by_cam_star<(1/(1-f))){
       coex_mat[i_f,i_ps]<-"COEX"
      }
    }
  }
  
  id_m_n_pos<-which(n_eqm_mat>0 & m_eqm_mat>0,arr.ind=T)
  id_m_n_pos<-as.data.frame(id_m_n_pos)
  
  if(nrow(id_m_n_pos)>0){
    # This f and ps combinations are good for co-exs of both mutualists
    f_ps_combo<-split(id_m_n_pos,id_m_n_pos$row)
    
    res<-c()
    for(i_f in c(1:length(f_ps_combo))){
      
      f<-fvary[unique(f_ps_combo[[i_f]]$row)]
      lenps<-length(f_ps_combo[[i_f]]$row)
      # initialize
      mat_max_eigen<-matrix(NA,nrow=lenps,ncol=7)
      colnames(mat_max_eigen)<-c("good_f","good_ps","ca_eqm","cc_eqm","m_eqm","n_eqm","max_eigen_val")
      
      mat_max_eigen[,1]<-f
      psrange<-psvary[f_ps_combo[[i_f]]$col]
      
      for(i_ps in c(1:length(psrange))){
        ps<-psrange[i_ps]
        mat_max_eigen[i_ps,2]<-ps
        
        eqmvals<-CaCcMN_eqm(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u)
        mat_max_eigen[i_ps,3]<-eqmvals$Caeq
        mat_max_eigen[i_ps,4]<-eqmvals$Cceq
        mat_max_eigen[i_ps,5]<-eqmvals$Meq
        mat_max_eigen[i_ps,6]<-eqmvals$Neq
        
        ans<-get_eigenvalues(km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,f=f,phi=phi,g=g,ps=ps,u=u,eqmvals=eqmvals,printstatus=printstatus)
        mat_max_eigen[i_ps,7]<-max(ans)
      }
      
      res<-rbind(res,mat_max_eigen)
      cat("=========i_f = ",i_f,"=========\n")
    }
    
    is_stable_eqm<-all(res[,7]<=0) # if TRUE then stable eqm as no positive eigen value
    return(list(res=res,
                is_stable_eqm=is_stable_eqm))
  }else{
    return(list(res="no coext.",
                is_stable_eqm="no coext."))
  }
}

#------------------------------------------------------------------------------------------------
# Now call the functions to see eigen value variation against fidelity with km=10, kn=6
#km<- 5 # half saturation constant for mutualist 
#kn<- 7 # half saturation constant for non-mutualist
#kc<- 15 # half saturation constant for allocated carbon
#d<- 0.5 # death rate of mutualist and non-mutualist
#bmax<- 0.8 # maximum growth rate of symbionts
#s<- 0.1 # cost of mutualism
#phi <- 5 # constant resource value for construction carbon
#g<- 0.2 # Rate at which construction carbon is allocated to both symbionts
#u<- 0.4 # Phosphorous uptake per unit of preferentially allocated carbon received by mutualists

# initialize
#fvary<-seq(from=0,to=1,by=0.01)
#psvary<-seq(from=0,to=1,by=0.01)

#ans<-eigen_f_ps_vary(fvary=fvary,psvary=psvary,km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,phi=phi,g=g,u=u,printstatus=T)
#ans$is_stable_eqm
#---------------------------------------------------------------------------------------------------

d<- 0.5 # death rate of mutualist and non-mutualist
bmax<- 0.8 # maximum growth rate of symbionts
s<- 0.1 # cost of mutualism
phi <- 5 # constant resource value for construction carbon
g<- 0.2 # Rate at which construction carbon is allocated to both symbionts
u<- 0.4 # Phosphorous uptake per unit of preferentially allocated carbon received by mutualists

# initialize
kmvary<-seq(from=1,to=20,by=1)
knvary<-seq(from=1,to=20,by=1)
kcvary<-seq(from=1,to=20,by=1)

fvary<-seq(from=0,to=1,by=0.01)
psvary<-seq(from=0,to=1,by=0.01)

sink("./Results/pdf_fig/eigen_res/km_kn_kc_f_ps_vary.txt", append=TRUE, split=TRUE)
stable_eqm<-c()  
for(i_km in c(1:length(kmvary))){
  for(i_kn in c(1:length(knvary))){
    for(i_kc in c(1:length(kcvary))){
      km<-kmvary[i_km]
      kn<-knvary[i_kn]
      kc<-kcvary[i_kc]
      ans<-eigen_f_ps_vary(fvary=fvary,psvary=psvary,km=km,kn=kn,kc=kc,d=d,bmax=bmax,s=s,phi=phi,g=g,u=u,printstatus=F)
      cat("------------km = ",km,"----------kn = ",kn,"--------kc = ",kc,"------stable eqm = ",ans$is_stable_eqm,"-----\n")
      stable_eqm<-c(stable_eqm,ans$is_stable_eqm)
    }
  }
}
sink()

saveRDS(stable_eqm,"./Results/pdf_fig/eigen_res/km_kn_kc_1_20_f_ps_vary.RDS")

any(stable_eqm==FALSE) # it should be False







