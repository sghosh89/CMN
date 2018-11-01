f<-0.3 # fidelity of plant C allocation to mutualist (M)
kc<-5 # half saturation constant 
ps<-0.3 # P availability in soil
u<- 0.4 # rate of P return by M

k_m<- 10# half saturation constant for mutualist 
k_n<- 20 # half saturation constant for non-mutualist
d<- 0.5 # death rate of mutualist and non-mutualist
bmax<- 0.8 #maximum growth rate of symbionts
s<- 0.3 # cost of mutualism
phi<- 5# initial resource for construction C : taken as constant and same for both symbionts

g<-0.2

#op<-par(mfrow=c(1,2))
plot(NA,xlim=c(0,5),ylim=c(0,150),xlab="M",ylab="C")

for(jn in 0:2){
  N<-jn*1
  #print(N)
  Ca<-c() # to strore the values for allocated carbon
  Cc<-c() # to strore the values for construction carbon
  Ca_m<-c()
  Ca_n<-c()
  M<-c()
  for(im in 1:50){
    m<-im*0.1
    temp<-(m+(N*(1-f)))*(m+kc)*(1-ps)
    temp<-temp/(u*(m^2))
    tempc<-phi/(g*(m+N))
    tempm<- -tempc + ((k_m*d)/((bmax*(1-s))-d))
    tempn<- -(tempc/(1-f)) + ((k_n*d)/((bmax-d)*((1-f)^2)))
    Ca<-c(Ca,temp)
    Cc<-c(Cc,tempc)
    M<-c(M,m)
    Ca_m<-c(Ca_m,tempm)
    Ca_n<-c(Ca_n,tempn)
    #abline(v=tempc,col="grey")
  } 
  
  lines(M,Ca,lty=jn+1)
  lines(M,Ca_m,col="red",lty=jn+1)
  lines(M,Ca_n,col="blue",lty=jn+1)
  lines(M,Cc,col="green4",lty=jn+1)
  #rect(xleft=min(Cc),ybottom=-100,xright=max(Cc),ytop=150,col=rgb(0,0.5,0,0.1),border = "white",lty=NULL)
}


s<(1-(d/bmax))                              #condtion for coexistence 


legend("topright",c("Mutualist", "Non-mutualist","Ca","Cc"), col = c("red","blue","black","green4"),bty="n",
       cex = 0.8, lty=c(1,1,1,1), xpd = TRUE, horiz = F, inset = c(0,-0.15),y.intersp = 0.8) 

legend("topright",c("N=0", "N=1","N=2"), bty="n",
       cex = 0.8, lty=c(1,2,3), xpd = TRUE, horiz = F, inset = c(0,0),y.intersp = 0.8) 


#----------------------------------------------------------------------
# 3D figure
#M    <- seq(from=1,to=30,by=2)
#N    <- seq(from=0,to=30,by=5)
#Cfun <- function(M,N){((M+(N*(1-f)))*(M+kc)*(1-ps))/((M^2)*u)}
#C    <- outer(M,N, FUN="Cfun")

#persp(M,N,C,theta = 20, phi = 20,col = "yellow", shade = 0.5,zlim=c(0,500))
#------------------------------------------------------------------------

op<-par(mfrow=c(1,3))
for(jn in 0:2){
  N<-jn*1
  #print(N)
  Cc<-c() # to strore the values for construction carbon
  M<-c()
  for(im in 1:50){
    m<-im*0.1
    tempc<-phi/(g*(m+N))
    Cc<-c(Cc,tempc)
    M<-c(M,m)
  } 
  
  a1 <-((k_m*d)/((bmax*(1-s))-d))
  a2<-((k_n*d)/((bmax-d)*((1-f)^2)))
  plot(NA,xlim=c(0,max(a1,a2)+20),ylim=c(0,max(a1,a2)+20),xlab="C_construction",ylab="C_allocation")
  lines(x = c(0,max(a1,a2)+20),y=c(0,max(a1,a2)+20),col="green")
  abline(a=a1,b=-1,col="red")
  abline(a =a2 , b=-1/(1-f),col="blue")
  mtext(paste0("N=",N),line=1,cex=1)
  grid()
  
  
  rect(xleft=min(Cc),ybottom=-100,xright=max(Cc),ytop=150,col=rgb(0,0.5,0,0.1),border = "white",lty=NULL)
}
par(op)
op<-par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 1, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("top", c("Mutualist","Non-mutualist"), col = c("red", "blue"),
       cex = 0.8, lty = c(1, 1), xpd = TRUE, horiz = T, inset = c(0,0),
       bty = "n") 
par(op)



