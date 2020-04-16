rm(list=ls())
graphics.off()

#***
#The functions
#***

#Given all the params except f, this function gives the acceptable range of f. The value
#given for fmax can be greater than 1, which means, when used, it should typically be 
#replaced with 1.
#
#Args
#PS, u, KA, D, aM, aN, bmax, s, K, d, e         Model parameters. Each a single number.
#
#Output
#A vector of length 3 containing fmin and fmax and the estimated precision of fmax
RangeOff<-function(PS, u, KA, D, aM, aN, bmax, s, K, d, e)
{
  fmin<-bmax*s/(bmax-d)
  
  Rhat<-function(f){return(K*d*(f-fmin)/(e*(aN-aM*(1-f))*(bmax-d)*(1-fmin)))}
  Ahat_alphahat<-function(f){return((K*d*(bmax-d)*(aN-aM)+s*K*d*bmax*aM)/((bmax-d)*(bmax*(1-s)-d)*(aN-aM*(1-f))))}
  Mhat<-function(f){return(((1-PS)+sqrt((1-PS)^2+4*u*Ahat_alphahat(f)*(1-PS)*KA))/(2*u*Ahat_alphahat(f)))}
  Nhat<-function(f){return(D/(aN*Rhat(f))-aM*Mhat(f)/aN)}
  
  if (Nhat(1)<0)
  {
    h<-uniroot(Nhat,c(fmin,1),maxiter=10000)
    fmax<-h$root 
    estim.prec<-h$estim.prec
    return(c(fmin=fmin,fmax=fmax,estim.prec=estim.prec))    
  }
  if (Nhat(10)<0)
  {
    h<-uniroot(Nhat,c(fmin,10),maxiter=10000)
    fmax<-h$root 
    estim.prec<-h$estim.prec
    return(c(fmin=fmin,fmax=fmax,estim.prec=estim.prec))    
  }
  if (Nhat(100)<0)
  {
    h<-uniroot(Nhat,c(fmin,100),maxiter=10000)
    fmax<-h$root 
    estim.prec<-h$estim.prec
    return(c(fmin=fmin,fmax=fmax,estim.prec=estim.prec))    
  }
  if (Nhat(1000)<0)
  {
    h<-uniroot(Nhat,c(fmin,1000),maxiter=10000)
    fmax<-h$root 
    estim.prec<-h$estim.prec
    return(c(fmin=fmin,fmax=fmax,estim.prec=estim.prec))    
  }
  if (Nhat(10000)<0)
  {
    h<-uniroot(Nhat,c(fmin,10000),maxiter=10000)
    fmax<-h$root 
    estim.prec<-h$estim.prec
    return(c(fmin=fmin,fmax=fmax,estim.prec=estim.prec))    
  }
  if (Nhat(100000)<0)
  {
    h<-uniroot(Nhat,c(fmin,100000),maxiter=10000)
    fmax<-h$root 
    estim.prec<-h$estim.prec
    return(c(fmin=fmin,fmax=fmax,estim.prec=estim.prec))    
  }
  stop("Error in RangeOff: 100000 was not big enough")
}

#Given all the params including f, write a function that gives the equilibrium of the model.
#
#Args
#PS, u, KA, D, aM, aN, bmax, s, K, d, e         Model parameters. Each a single number.
#f      Another model parameter, again a single number
#
#Output
#A vector of Ahat, Rhat, Mhat, Nhat
GetEquil<-function(PS, u, KA, D, aM, aN, bmax, s, K, d, e, f)
{
  fRange<-unname(RangeOff(PS, u, KA, D, aM, aN, bmax, s, K, d, e))
  #if (f<=fRange[1] || f>=fRange[2])
  #{
  #  stop("Error in GetEquil: f not in acceptable range")
  #}
  
  fmin<-fRange[1]
  Rhat<-K*d*(f-fmin)/(e*(aN-aM*(1-f))*(bmax-d)*(1-fmin))
  Ahat_alphahat<-(K*d*(bmax-d)*(aN-aM)+s*K*d*bmax*aM)/((bmax-d)*(bmax*(1-s)-d)*(aN-aM*(1-f)))
  Mhat<-((1-PS)+sqrt((1-PS)^2+4*u*Ahat_alphahat*(1-PS)*KA))/(2*u*Ahat_alphahat)
  Nhat<-D/(aN*Rhat)-aM*Mhat/aN
  Ahat<-Ahat_alphahat*(Mhat+(1-f)*Nhat)
  
  return(c(Ahat=Ahat,Rhat=Rhat,Mhat=Mhat,Nhat=Nhat,Ahat_alphahat=Ahat_alphahat))
}

#Given all the params including f, computes the Jacobian at the equilibrium. 
#
#Args
#PS, u, KA, D, aM, aN, bmax, s, K, d, e         Model parameters. Each a single number.
#
#Output - a list containing the following named elements
#equil            The equilibrium - Ahat, Rhat, Mhat, Nhat
#jac              The Jacobian at the equilibrium
#max_re_part      The maximal real part of all Jacobian eigenvales
GetJac<-function(PS, u, KA, D, aM, aN, bmax, s, K, d, e, f)
{
  #Get the equilirium
  equil<-GetEquil(PS, u, KA, D, aM, aN, bmax, s, K, d, e, f)
  A<-unname(equil["Ahat"])
  R<-unname(equil["Rhat"])
  M<-unname(equil["Mhat"])
  N<-unname(equil["Nhat"])
  equil<-equil[1:4]
  alph<-M+(1-f)*N
  
  #Get the Jacobian, derivs computed by hand
  dfAdA<-(-u*M^2)/((KA+M)*alph)
  dfAdR<-0
  dfAdM<-(-u*A/alph)*(((KA+M)*2*M-M^2)/((KA+M)^2))+(M^2*u*A)/((KA+M)*alph^2)
  dfAdN<-(u*M^2*A*(1-f))/((KA+M)*alph^2)
    
  dfRdA<-0
  dfRdR<-(-(aM*M+aN*N))
  dfRdM<-(-aM*R)
  dfRdN<-(-aN*R)

  CM<-e*aM*R+A/alph
  dfMdA<-(K*bmax*(1-s)*M)/(alph*(K+CM)^2)
  dfMdR<-(K*bmax*(1-s)*e*aM*M)/((K+CM)^2)
  dfMdM<-(bmax*(1-s)*(K*CM-K*A*M/alph^2+CM^2))/((K+CM)^2)-d
  dfMdN<-(-bmax*(1-s)*M*A*(1-f)*K)/(alph^2*(K+CM)^2)

  CN<-e*aN*R+(1-f)*A/alph
  dfNdA<-(bmax*N*(1-f)*K)/(alph*(K+CN)^2)
  dfNdR<-(bmax*N*e*aN*K)/((K+CN)^2)
  dfNdM<-(-bmax*N*(1-f)*A*K)/(alph^2*(K+CN)^2)
  dfNdN<-((K+CN)*bmax*CN-K*bmax*N*(1-f)^2*A/alph^2)/((K+CN)^2)-d
  
  jac<-matrix(c(dfAdA,dfAdR,dfAdM,dfAdN,
                dfRdA,dfRdR,dfRdM,dfRdN,
                dfMdA,dfMdR,dfMdM,dfMdN,
                dfNdA,dfNdR,dfNdM,dfNdN),4,4,byrow = TRUE)
  
  #Now take the maximum of the real parts of the eigenvalues
  max_re_part<-max(Re(eigen(jac,symmetric=FALSE,only.values=TRUE)$values))
  
  #Assemble results and return
  return(list(equil=equil,jac=jac,max_re_part=max_re_part))
}

#Another way of getting the Jacobian, for comparison
fA<-expression(1-PS-u*(M^2/(KA+M))*A/(M+(1-f)*N))
fR<-expression(D-(aM*M+aN*N)*R)
fM<-expression(((bmax*(1-s)*(e*aM*R+A/(M+(1-f)*N)))/(K+e*aM*R+A/(M+(1-f)*N))-d)*M)
fN<-expression(((bmax*(e*aN*R+(1-f)*A/(M+(1-f)*N)))/(K+e*aN*R+(1-f)*A/(M+(1-f)*N))-d)*N)

dfAdA<-stats::D(fA,"A")
dfAdR<-stats::D(fA,"R")
dfAdM<-stats::D(fA,"M")
dfAdN<-stats::D(fA,"N")

dfRdA<-stats::D(fR,"A")
dfRdR<-stats::D(fR,"R")
dfRdM<-stats::D(fR,"M")
dfRdN<-stats::D(fR,"N")

dfMdA<-stats::D(fM,"A")
dfMdR<-stats::D(fM,"R")
dfMdM<-stats::D(fM,"M")
dfMdN<-stats::D(fM,"N")

dfNdA<-stats::D(fN,"A")
dfNdR<-stats::D(fN,"R")
dfNdM<-stats::D(fN,"M")
dfNdN<-stats::D(fN,"N")

GetJac2<-function(PS, u, KA, D, aM, aN, bmax, s, K, d, e, f)
{
  #Get the equilirium
  equil<-GetEquil(PS, u, KA, D, aM, aN, bmax, s, K, d, e, f)
  A<-unname(equil["Ahat"])
  R<-unname(equil["Rhat"])
  M<-unname(equil["Mhat"])
  N<-unname(equil["Nhat"])
  equil<-equil[1:4]

  jac<-matrix(c(eval(dfAdA),eval(dfAdR),eval(dfAdM),eval(dfAdN),
                eval(dfRdA),eval(dfRdR),eval(dfRdM),eval(dfRdN),
                eval(dfMdA),eval(dfMdR),eval(dfMdM),eval(dfMdN),
                eval(dfNdA),eval(dfNdR),eval(dfNdM),eval(dfNdN)),4,4,byrow = TRUE)
  
  #Now take the maximum of the real parts of the eigenvalues
  max_re_part<-max(Re(eigen(jac,symmetric=FALSE,only.values=TRUE)$values))
  
  #Assemble results and return
  return(list(equil=equil,jac=jac,max_re_part=max_re_part))
}

#***
#Tests
#***

#Params from Shya, she wrote them in the paper
u<-.4 #
bmax<-.8 #
d<-.5 #
s<-.1 #
e<-.5 #
aM<-.1 #
aN<-.2 #
D<-5 #
KA<-5 #
K<-10 #
PS<-.3 #

#get the allowed range for f and check with Shya results
rg<-RangeOff(PS, u, KA, D, aM, aN, bmax, s, K, d, e)
rg #Agrees with Shya

#now make plots of how the equilibrium depends on f
inputf<-seq(from=rg[1],to=rg[2],length.out=102)
out<-matrix(NA,5,length(inputf))
for (counter in 1:length(inputf))
{
  out[,counter]<-GetEquil(PS, u, KA, D, aM, aN, bmax, s, K, d, e, inputf[counter])
}

out[is.infinite(out)]<-NA
plot(inputf,out[1,],type="l",ylim=range(out[1:2,],na.rm=TRUE),xlab="f",ylab="Ahat or Rhat")
lines(inputf,out[2,],type="l",lty="dashed")
legend(x="topright",legend=c("Ahat","Rhat"),lty=c("solid","dashed")) #Agrees with Shya plot

plot(inputf,out[3,],type="l",ylim=range(out[3:4,],na.rm=TRUE),xlab="f",ylab="Mhat or Nhat")
lines(inputf,out[4,],type="l",lty="dashed")
lines(inputf,rep(0,length(inputf)),type="l",col="red")
legend(x="topright",legend=c("Mhat","Nhat","zero"),lty=c("solid","dashed","solid"),col=c("black","black","red")) #Agrees with Shya plot

plot(inputf,out[3,],type="l",ylim=range(out[3,],na.rm=TRUE),xlab="f",ylab="Mhat")

#for inspection
out

#Jacobian
h1<-GetJac(PS, u, KA, D, aM, aN, bmax, s, K, d, e, .32)
h2<-GetJac2(PS, u, KA, D, aM, aN, bmax, s, K, d, e, .32)
h1$equil
h2$equil
max(abs(h1$equil-h2$equil))
h1$jac
h2$jac
max(abs(h1$jac-h2$jac))
h1<-GetJac(PS, u, KA, D, aM, aN, bmax, s, K, d, e, .42)
h2<-GetJac2(PS, u, KA, D, aM, aN, bmax, s, K, d, e, .42)
h1$equil
h2$equil
max(abs(h1$equil-h2$equil))
h1$jac
h2$jac
max(abs(h1$jac-h2$jac))

out<-NA*numeric(length(inputf))
for (counter in 2:length(inputf))
{
  out[counter]<-GetJac(PS, u, KA, D, aM, aN, bmax, s, K, d, e, inputf[counter])$max_re_part
}
plot(inputf,out,type="p",pch=20,cex=.5,xlim=range(inputf))
lines(inputf,out,type="l")
lines(range(inputf),rep(0,2),type="l",col="red")
lines(rep(rg[1],2),range(out,na.rm=TRUE),type="l",lty="dashed")
lines(rep(rg[2],2),range(out,na.rm=TRUE),type="l",lty="dashed")

#compare to another way of computing the Jacobian
out<-NA*numeric(length(inputf))
for (counter in 2:length(inputf))
{
  out[counter]<-GetJac2(PS, u, KA, D, aM, aN, bmax, s, K, d, e, inputf[counter])$max_re_part
}
plot(inputf,out,type="p",pch=20,cex=.5,xlim=range(inputf))
lines(inputf,out,type="l")
lines(range(inputf),rep(0,2),type="l",col="red")
lines(rep(rg[1],2),range(out,na.rm=TRUE),type="l",lty="dashed")
lines(rep(rg[2],2),range(out,na.rm=TRUE),type="l",lty="dashed")

#***
#Now apply the above functions to examine stability
#***

#"soft" range of parameters to consider
rg_u<-c(.001,1)
rg_bmax<-c(.001,5) 
rg_d<-c(.001,5) 
rg_s<-c(.001,.999) 
rg_e<-c(.001,.999) 
rg_aM<-c(.001,5)
rg_aN<-c(.001,5) 
rg_D<-c(.001,5) 
rg_KA<-c(.01,15) 
rg_K<-c(.01,15) 
rg_PS<-c(.001,.999) 

lowbd<-c(rg_u[1],
         rg_bmax[1],
         rg_d[1],
         rg_s[1],
         rg_e[1],
         rg_aM[1],
         rg_aN[1],
         rg_D[1],
         rg_KA[1],
         rg_K[1],
         rg_PS[1])
pdiff<-c(diff(rg_u),
         diff(rg_bmax),
         diff(rg_d),
         diff(rg_s),
         diff(rg_e),
         diff(rg_aM),
         diff(rg_aN),
         diff(rg_D),
         diff(rg_KA),
         diff(rg_K),
         diff(rg_PS))

numruns<-1000000 #number of sets of non-f parameters to try
numf<-25 #number of attempted f values per set of non-f parameters
set.seed(101)
randnums<-matrix(runif(numruns*length(lowbd)),numruns,length(lowbd)) #generate all random nums at outset
results<-matrix(NA,numruns*numf,length(lowbd)+9)
colnames(results)<-c("u","bmax","d","s","e","aM","aN","D","KA","K","PS", #parameters
                     "fmin","fmax","fmax.precision","f", #fmin and fmax
                     "Ahat","Rhat","Mhat","Nhat", #equilibrium
                     "max_re_evals") #max real part of Jacobian eigenvalues at equil
resrowcounter<-1
for (counter in 1:numruns)
{
  #Get the non-f parameters to use for this time through this outer loop
  parms<-randnums[counter,]
  parms<-parms*pdiff+lowbd
  u<-parms[1]
  bmax<-parms[2]
  d<-parms[3]
  s<-parms[4]
  e<-parms[5]
  aM<-parms[6]
  aN<-parms[7]
  D<-parms[8]
  KA<-parms[9]
  K<-parms[10]
  PS<-parms[11]
  
  #See if all the hard bounds are satisfied, if not don't proceed. There is some redundancy here but 
  #I don't care. I am keeping it for robustness.
  if (aN<aM || aM<=0)
  {
    next
  }
  if (s<=0)
  {
    next
  }
  if (bmax-d<=0)
  {
    next
  }
  if (bmax*(1-s)-d<=0)
  {
    next
  }
  if (PS<=0 || PS>=1)
  {
    next
  }
  if (D<=0)
  {
    next
  }
  if (e<=0 || e>=1)
  {
    next
  }
  if (u<=0 || u>=1)
  {
    next
  }
  if (KA<=0)
  {
    next
  }
  if (K<=0)
  {
    next
  }
  
  #Get the range of f and the specific f you will use
  rgf<-RangeOff(PS, u, KA, D, aM, aN, bmax, s, K, d, e)
  if (rgf[2]>1)
  {
    rgf[2]<-1
  }
  if (unname((rgf[2]-rgf[1])/rgf[3]<2*numf))
  {
    #if the gap between rgf[1] and rgf[2] is not at least 2*(numf+1) times the precision rgf[3], don't proceed
    next
  }
  these_fs<-seq(from=rgf[1],to=rgf[2],length.out=numf+2)
  these_fs<-these_fs[-c(1,length(these_fs))]

  #If you get to here, you have valid parameter combinations, so save the values in the next
  #several rows of results...
  results[resrowcounter:(resrowcounter+length(these_fs)-1),1:11]<-rep(parms,each=length(these_fs))
  results[resrowcounter:(resrowcounter+length(these_fs)-1),12:14]<-rep(unname(rgf),each=length(these_fs))
  
  #...and put the f values to be used into the correct column of results
  results[resrowcounter:(resrowcounter+length(these_fs)-1),15]<-these_fs
  
  #Now get the equilibrium, and the Jacobian, and evaluate stability for each f
  for (fcount in 1:numf)
  {
    h<-GetEquil(PS, u, KA, D, aM, aN, bmax, s, K, d, e, these_fs[fcount])
    results[resrowcounter+fcount-1,16:19]<-h[1:4]
    
    h<-GetJac(PS, u, KA, D, aM, aN, bmax, s, K, d, e, these_fs[fcount])
    results[resrowcounter+fcount-1,20]<-h$max_re_part
  }
  
  resrowcounter<-resrowcounter+numf
}

#now trim the left-over rows which are full of NAs
results<-results[!is.na(results[,1]),]
dim(results)

#now look at stability results and capture the ones that are supposedly unstable
hist(results[,"max_re_evals"],50)
max(results[,"max_re_evals"])
results<-results[order(results[,"max_re_evals"],decreasing=TRUE),]
head(results)
