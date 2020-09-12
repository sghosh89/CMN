
source("get_MNAR_eqm_analytical.R")
#===========================================

#==========================================================================
# Now, call the function
#d<-0.5
#bmax<-0.8

#s<-0.1
#fmin<-(s*bmax)/(bmax-d)
#fmin
#curve(get_MNAR_eqm_analytical, xlim=c(fmin,1), lwd=2, lty=2, ylab='Neq',xlab='f',ylim=c(-0.1,1))
#abline(h=0)
#myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,0.9),ps=0.3,s=0.1,aM=0.1,aN=0.2,phi=5,getalleqmval=F) 
#fmax<-myfmax$root #
#cat("range of fidelity for coexistence: (fmin,fmax)=(",fmin,",",fmax,")\n")
#(fmax-fmin)

#------------------ when aN=aM=0.1, s=0.1
#s<-0.1
#fmin<-(s*bmax)/(bmax-d)
#myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,0.9), ps=0.3,s=0.1,aM=0.1,aN=0.1,phi=5,getalleqmval=F) #fmax=0.3146655
#fmax<-myfmax$root
#cat("range of fidelity for coexistence: (fmin,fmax)=(",fmin,",",fmax,")\n")
#(fmax-fmin)

#----------------- when s=0, aN>aM 
#s<-0
#fmin<-(s*bmax)/(bmax-d)
#myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,0.9), ps=0.3, s=s, aM=0.1, aN=0.2, phi=5, getalleqmval=F) #fmax=0.2103166
#fmax<-myfmax$root
#cat("range of fidelity for coexistence: (fmin,fmax)=(",fmin,",",fmax,")\n")
#(fmax-fmin)

#==========================================================================

#========================== generate plot with variation of s, aM<aN, phi=5 and phi=1 ==================

bmax<-0.8
d<-0.5
#smax<-1-(d/bmax)
smax<-0.31
s<-0
fmins<-c() 
fmaxs_phi5_aM_0.1_aN_0.2<-c()
fmaxs_phi1_aM_0.1_aN_0.2<-c()

s_vec<-c()
while(s<=smax){
  cat("s=",s,"\n")
  s_vec<-c(s_vec,s)
  fmin<-(s*bmax)/(bmax-d) # indep. of phi, aM or aN
  fmins<-c(fmins,fmin)
  
  #---------- aM=0.1, aN=0.2 ---------------
  
  myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin+0.0001,0.999), ps=0.3, s=s, aM=0.1, aN=0.2, phi=5, getalleqmval=F) 
  fmax_phi5<-myfmax$root
  fmaxs_phi5_aM_0.1_aN_0.2<-c(fmaxs_phi5_aM_0.1_aN_0.2,fmax_phi5)
  
  myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin+0.0001,0.999), ps=0.3, s=s, aM=0.1, aN=0.2, phi=1, getalleqmval=F) 
  fmax_phi1<-myfmax$root
  fmaxs_phi1_aM_0.1_aN_0.2<-c(fmaxs_phi1_aM_0.1_aN_0.2,fmax_phi1)
  
  s<-s+0.01
}

resloc<-"./ARMN_Results/"
pdf(paste(resloc,"range_of_fidelity_vs_s_aM_0.1_aN_0.2.pdf",sep=""),width=8,height=8)

op<-par(mar=c(6,6,2,2),pty="s",family="serif")

plot(s_vec,fmins,ylim=c(0,1),type="l",lwd=2,ylab="",xlab="s",cex.lab=2.5,cex.axis=2)
lines(s_vec,fmaxs_phi1_aM_0.1_aN_0.2,col="red",lwd=2)
lines(s_vec,fmaxs_phi5_aM_0.1_aN_0.2,col="blue",lwd=2)
legend("topleft", 
       c(expression(paste("f"["max"], " , D=5")),
         expression(paste("f"["max"], " , D=1")),
         expression("f"["min"])), 
       col = c("blue", "red", "black"),
       cex = 2.5, lty = c(1, 1, 1), lwd=c(2,2,2), xpd = TRUE, 
       horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
       bty = "n") 
legend("bottomright", c(expression(paste("a"["M"],"=0.1, ")),expression(paste("a"["N"],"=0.2"))),
       cex = 2.5,
       horiz = T, x.intersp = -0.6,
       bty = "n")

par(op)
dev.off()

#---------------------------------------------------------------

#---------- aM=0.1, aN=0.1 ---------------

bmax<-0.8
d<-0.5
fmins<-c() 
fmaxs_phi5_aM_0.1_aN_0.1<-c()
fmaxs_phi1_aM_0.1_aN_0.1<-c()
smax<-0.31
s<-0
s_vec<-c()
while(s<=smax){
#  if(s>0){
    cat("s=",s,"\n")
    s_vec<-c(s_vec,s)
    
    fmin<-(s*bmax)/(bmax-d) # indep. of phi, aM or aN
    fmins<-c(fmins,fmin)
    
    myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin+0.001,0.9), ps=0.3, s=s, aM=0.1, aN=0.1, phi=5, getalleqmval=F) 
    fmax_phi5<-myfmax$root
    fmaxs_phi5_aM_0.1_aN_0.1<-c(fmaxs_phi5_aM_0.1_aN_0.1,fmax_phi5)
    
    myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin+0.001,0.9), ps=0.3, s=s, aM=0.1, aN=0.1, phi=1, getalleqmval=F) 
    fmax_phi1<-myfmax$root
    fmaxs_phi1_aM_0.1_aN_0.1<-c(fmaxs_phi1_aM_0.1_aN_0.1,fmax_phi1)
#  }
  s<-s+0.01
}

resloc<-"./ARMN_Results/"
pdf(paste(resloc,"range_of_fidelity_vs_s_aM_0.1_aN_0.1.pdf",sep=""),width=8,height=8)

op<-par(mar=c(6,6,2,2),pty="s",family="serif")

plot(s_vec,fmins,ylim=c(0,1),type="l",lwd=2,ylab="",xlab="s",cex.lab=2.5,cex.axis=2)
lines(s_vec,fmaxs_phi1_aM_0.1_aN_0.1,col="red",lwd=2)
lines(s_vec,fmaxs_phi5_aM_0.1_aN_0.1,col="blue",lwd=2)
legend("topleft", 
       c(expression(paste("f"["max"], " , D=5")),
         expression(paste("f"["max"], " , D=1")),
         expression("f"["min"])), 
       col = c("blue", "red", "black"),
       cex = 2.5, lty = c(1, 1, 1), lwd=c(2,2,2), xpd = TRUE, 
       horiz = F, inset = c(0,0),y.intersp = 0.8,x.intersp = 0.1,
       bty = "n") 
legend("bottomright", c(expression(paste("a"["M"],"=0.1,")),expression(paste("a"["N"],"=0.1"))),
       cex = 2.5,
       horiz = T, x.intersp = -0.6,
       bty = "n")


par(op)
dev.off()



#======================================================================================================

#------------------ Contour plots suggested by Am. Nat.'s Erol Akcay -----------------------

# contour plot of Meq/(Meq+Neq) vs. fidelity,f vs. constant new root growth rate, phi

get_phi_f_PM_data<-function(aM,aN,s){
  
  # start of a function here
  ps<-0.3
  bmax<-0.8
  d<-0.5
  fmin<-(s*bmax)/(bmax-d)
  
  phi_range<-seq(from=1,to=5,by=1/50)
  
  phi_and_fmax<-data.frame(phi=phi_range,fmax=NA)
  
  for(i in 1:length(phi_range)){
    phi<-phi_range[i]
    cat("phi=",phi,"\n")
    
    myfmax<-uniroot(get_MNAR_eqm_analytical, interval=c(fmin,200), # ideally the interval we should look (fmin,1) 
                    ps=ps,s=s,aM=aM,aN=aN,phi=phi,getalleqmval=F) 
    fmax<-myfmax$root
    
    # if(fmax>1){
    #  fmax<-1 # force set fmax = 1 
    #stop("Stop code: fmax is not less that 1 for a given ps",call.=T)
    #}
    
    
    phi_and_fmax$fmax[i]<-fmax
  }
  
  
  # Now make a table for all possible combo of phi, f for co-existence 
  # We varied f in between >= fmin to < fmax for a given phi
  
  len<-15000 # give it a big number
  
  phi_f_PM<-data.frame(phi=NA*numeric(len),
                       f=NA*numeric(len),
                       PM=NA*numeric(len),
                       Aeq=NA*numeric(len),
                       Req=NA*numeric(len),
                       Meq=NA*numeric(len),
                       Neq=NA*numeric(len))
  
  k<-1
  finit<-fmin#0
  for(i in c(1:nrow(phi_and_fmax))){
    phi<-phi_and_fmax$phi[i]
    #fmax<-1
    fmax<-phi_and_fmax$fmax[i] 
    
    #cat("ps=",ps,"fmax=",fmax,"\n")
    
    f_range<-seq(from=finit,to=fmax,by=(fmax-finit)/50)
    
    for(j in c(1:length(f_range))){
      
      f<-f_range[j]
      
      ans<-get_MNAR_eqm_analytical(f=f,ps=ps,s=s,aM=aM,aN=aN,phi=phi,getalleqmval=T)
      PM<-ans$Meq/(ans$Meq+ans$Neq)
      
      phi_f_PM$f[k]<-f
      phi_f_PM$phi[k]<-phi
      phi_f_PM$PM[k]<-PM
      phi_f_PM$Aeq[k]<-ans$Aeq
      phi_f_PM$Req[k]<-ans$Req
      phi_f_PM$Meq[k]<-ans$Meq
      phi_f_PM$Neq[k]<-ans$Neq
      cat("k=",k,"f=",f,"phi=",phi,"PM=",PM,"\n")
      k<-k+1
    }
  }
  
  (z<-phi_f_PM[which(phi_f_PM$PM>1),]) # this table shows at f=fmax, PM should be 1 but instead it's slightly >1
  # I think, this is because uniroot function just finds the root (that could be 0,>0,<0) and these are numerical precession error
  # but for meaningful biological variable Neq can't be negative, so we can consider 
  # Neq goes to zero and PM = 1 for given phi, f combination
  
  ind<-which(phi_f_PM$PM>1)
  phi_f_PM[ind,]$PM<-1
  
  phi_f_PM<-na.omit(phi_f_PM) # just to omit unnecessary blank lines
  
  return(phi_f_PM)
}

# call the function now

#----------------------------- aM = aN -----------------------------------------------------

# for aM = aN, s=0.025
aM<-0.1
aN<-0.1
s<-0.025
ans<-get_phi_f_PM_data(aM=aM,aN=aN,s=s)
write.csv(ans,paste("./ARMN_Results/phi_f_PM_aM_",aM,"_aN_",aN,"_s_",s,".csv",sep=""),row.names = F)


# for aM = aN, s=0.1
#aM<-0.1
#aN<-0.1
#s<-0.1
#ans<-get_phi_f_PM_data(aM=aM,aN=aN,s=s)
#write.csv(ans,paste("./ARMN_Results/phi_f_PM_aM_",aM,"_aN_",aN,"_s_",s,".csv",sep=""),row.names = F)

# for aM = aN, s=0.2
#aM<-0.1
#aN<-0.1
#s<-0.2
#ans<-get_phi_f_PM_data(aM=aM,aN=aN,s=s)
#write.csv(ans,paste("./ARMN_Results/phi_f_PM_aM_",aM,"_aN_",aN,"_s_",s,".csv",sep=""),row.names = F)

# for aM = aN, s=0.25
aM<-0.1
aN<-0.1
s<-0.25
ans<-get_phi_f_PM_data(aM=aM,aN=aN,s=s)
write.csv(ans,paste("./ARMN_Results/phi_f_PM_aM_",aM,"_aN_",aN,"_s_",s,".csv",sep=""),row.names = F)

#----------------------------- aM < aN -----------------------------------------------------

# for aM < aN, s=0.025
aM<-0.1
aN<-0.2
s<-0.025
ans<-get_phi_f_PM_data(aM=aM,aN=aN,s=s)
write.csv(ans,paste("./ARMN_Results/phi_f_PM_aM_",aM,"_aN_",aN,"_s_",s,".csv",sep=""),row.names = F)

# for aM < aN, s=0.1
#aM<-0.1
#aN<-0.2
#s<-0.1
#ans<-get_phi_f_PM_data(aM=aM,aN=aN,s=s)
#write.csv(ans,paste("./ARMN_Results/phi_f_PM_aM_",aM,"_aN_",aN,"_s_",s,".csv",sep=""),row.names = F)

# for aM < aN, s=0.2
#aM<-0.1
#aN<-0.2
#s<-0.2
#ans<-get_phi_f_PM_data(aM=aM,aN=aN,s=s)
#write.csv(ans,paste("./ARMN_Results/phi_f_PM_aM_",aM,"_aN_",aN,"_s_",s,".csv",sep=""),row.names = F)

# for aM < aN, s=0.25
aM<-0.1
aN<-0.2
s<-0.25
ans<-get_phi_f_PM_data(aM=aM,aN=aN,s=s)
write.csv(ans,paste("./ARMN_Results/phi_f_PM_aM_",aM,"_aN_",aN,"_s_",s,".csv",sep=""),row.names = F)
















