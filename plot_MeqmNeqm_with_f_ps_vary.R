source("./get_MNAR_eqm.R")
set.seed(seed=101)
PUfun<-function(f,ps,KM,KN,s,d,bmax,vareqm){
  ans<-nleqslv(x=c(1,1,10,40),fn=get_MNAR_eqm,KM=KM,KN=KN,f=f,ps=ps,aM=0.1,aN=0.2,KA=5,d=d,bmax=bmax,s=s,phi=5,u=0.4,eM=0.5,eN=0.5)
  #cat(" f = ",f," ------ ps = ",ps," -------","\n")
  if(vareqm=="M"){
    res<-ans$x[1]
  }else if(vareqm=="N"){
    res<-ans$x[2]
  }else{
    stop("error: vareqm should be 'M'or 'N'")
  }
  return(res)
}


KM<- 10 # half saturation constant for mutualist 
KN<- 10 # half saturation constant for non-mutualist

bmax<-0.8
s<-0.1
d<-0.5
tempo<-((bmax*(1-s))-d)/(bmax-d)
fmin<-1-((KN/KM)*(tempo))

f_seq<-seq(from=round(fmin,2),to=1,by=0.01)
lenf<-length(f_seq)
ps_seq<-seq(from=0,to=1,by=0.01)
lenps<-length(ps_seq)
PUM<-matrix(NA,nrow=lenf,ncol=lenps)
PUN<-PUM
PU_PM<-PUM # to store proportion of mutualist at eqm.

for(i in c(1:lenf)){ #loop for f
  
  for(j in 1:lenps){ #loop for ps
    f<-f_seq[i]
    ps<-ps_seq[j]
    resM<-PUfun(f=f,ps=ps,KM=KM,KN=KN,s=s,d=d,bmax=bmax,vareqm="M")
    resN<-PUfun(f=f,ps=ps,KM=KM,KN=KN,s=s,d=d,bmax=bmax,vareqm="N")
    
    if(resN<0){
      cat(" f = ",f," ------ ps = ",ps," -------Meqm = ",resM,"-------Neqm = ",resN,"\n")
    }
    
    if(resM<0){
      resM<-NA
    }
    if(resN<0){
      resN<-NA
    }
    
    PUM[i,j]<-resM
    PUN[i,j]<-resN
    PU_PM[i,j]<-resM/(sum(resM,resN,na.rm=T))
  }
}

range(PUM,na.rm=T)
range(PUN,na.rm=T)
range(PU_PM,na.rm=T)


pdf("./ARMN_Results/Meqm_vs_f_vs_ps_3dplot_KM_10_KN_10.pdf",height=8,width=8)
op<-par(mar=c(2,2,2,2),pty="s",mgp=c(4,0.5,0))
persp(f_seq,ps_seq,PUM,theta = -60, phi = 30,col = "lightblue",xlab="fidelity",
      ylab="soil P",
      zlab="M eqm.",expand=0.5,ltheta = 150, shade = 0.5, ticktype = "detailed",
      cex.lab=1.5,cex.axis=1.5,box=T,axes=T,xlim=range(f_seq),ylim=range(ps_seq))
par(op)
dev.off()

pdf("./ARMN_Results/Neqm_vs_f_vs_ps_3dplot_KM_10_KN_10.pdf",height=8,width=8)
op<-par(mar=c(2,2,2,2),pty="s",mgp=c(4,1,0))
persp(f_seq,ps_seq,PUN,theta = -30, phi = 30,col = "lightgreen",xlab="fidelity",
      ylab="soil P",
      zlab="N eqm.",expand=0.5,ltheta = 150, shade = 0.75, ticktype = "detailed",
      cex.lab=1.5,cex.axis=1.5,box=T,axes=T,xlim=range(f_seq),ylim=range(ps_seq))
par(op)
dev.off()

pdf("./ARMN_Results/Proportion_of_M_eqm_vs_f_vs_ps_3dplot_KM_10_KN_10.pdf",height=8,width=8)
op<-par(mar=c(2,2,2,2),pty="s",mgp=c(4,1,0))
persp(f_seq,ps_seq,PU_PM,theta = -30, phi = 30,col = "lightgreen",xlab="fidelity",
      ylab="soil P",
      zlab="PM eqm.",expand=0.5,ltheta = 150, shade = 0.75, ticktype = "detailed",
      cex.lab=1.5,cex.axis=1.5,box=T,axes=T,xlim=range(f_seq),ylim=range(ps_seq),zlim=c(0,1))
par(op)
dev.off()

id<-which(!is.finite(PUN),arr.ind=T)
