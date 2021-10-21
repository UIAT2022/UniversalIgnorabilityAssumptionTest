####################################################################################################### ########################################
# Simulation of scenario 5 (Compliance Class Model)
# This code is the modified version of the original code of Guo et al. (2014)
#####################################################################################################################################
#install.packages("AER")
#install.packages("MASS")
#install.packages("mlogit")
library(AER)
library(MASS)
library(mlogit)

#Simulation Input
library(AER)
library(MASS)
library(mlogit)
tol<-10^(-6)      #maximum tolerated distance between two iterations for EM 
maxiter<-1000     #maximum iteration number


#EM algorithm maximize over whole parameter space
EM<-function(data,delta_nt_0,tau_nt_0,delta_at_0,tau_at_0,lambda_at_0,lambda_co_0,lambda_nt_0,logodd_0,sigma_0,alpha_ra_0,beta_ra_0,tol,maxiter){   
  y<-data[,1]                         # input the simulated data
  z<-data[,2]
  a<-data[,3]
  x<-data[,5]
  cova=cbind(1,x)
  n<-length(y)
  # insert the initial value of the EM algorithm
  
  delta_nt<-delta_nt_0               # initial parameters for compliance classes model
  tau_nt<-tau_nt_0
  delta_at<-delta_at_0
  tau_at<-tau_at_0
  lambda_at<-lambda_at_0             # initial parameters for outcome model and it is inserted with vectors              
  lambda_co<-lambda_co_0
  lambda_nt<-lambda_nt_0
  logodd<-logodd_0
  sigma<-sigma_0
  alpha_ra<-alpha_ra_0               # initial parameters for random assignment model
  beta_ra<-beta_ra_0
  
  flag<-T
  ittime=0
  
  while(flag){
    #E step using the old parameter to compute the expectation
    denominator=1+exp(delta_nt+tau_nt*x)+exp(delta_at+tau_at*x)     # the denominator matrix of P(Ci|Xi)
    phi_n=exp(delta_nt+tau_nt*x)/denominator                        #P(Ci=nt|xi) which is a matrix recording the prob for subjects from 1 to n 
    phi_a=exp(delta_at+tau_at*x)/denominator                        #P(Ci=at|xi)
    phi_c = 1 /denominator 
    w1 <- matrix(0,n,1)
    w2 <- matrix(0,n,1)
    pc1=(1/sigma)*exp(-(y-cova%*%t(lambda_co)-cova%*%t(logodd))^2/(sigma^2))                #p(Yi=1|Zi=1,Ci=co,xi)
    pat=(1/sigma)*exp(-(y-cova%*%t(lambda_at))^2/(sigma^2))                                 #p(Yi=1|Ci=at,xi)
    pnt=(1/sigma)*exp(-(y-cova%*%t(lambda_nt))^2/(sigma^2))                                 #p(Yi=1|Ci=nt,xi)
    pc0=(1/sigma)*exp(-(y-cova%*%t(lambda_co))^2/(sigma^2))                               #p(Yi=1|Zi=0,Ci=co,xi)
    w1<-phi_c*pc1/(phi_c*pc1+phi_a*pat)
    w2<-phi_c*pc0/(phi_c*pc0+phi_n*pnt)
    
    
    #M step 
    obj <- function (para){                                                   # define the expected likelihood function 
      delta_nt<-para[1]                                                         # parameters for compliance classes model
      tau_nt<-para[2] 
      delta_at<-para[3] 
      tau_at<-para[4] 
      lambda_at<-cbind(para[5],para[6])                                                        #  parameters for outcome model              
      lambda_co<-cbind(para[7],para[8])  
      lambda_nt<-cbind(para[9],para[10])  
      logodd<-cbind(para[11],para[12])  
      sigma<-para[13]
      alpha_ra<-para[14]                                                        # parameters for random assignment model
      beta_ra<-para[15] 
      
      denominator=1+exp(delta_nt+tau_nt*x)+exp(delta_at+tau_at*x)      # the denominator matrix of P(Ci|Xi) of para
      phi_n=exp(delta_nt+tau_nt*x)/denominator                        #P(Ci=nt|xi) which is a matrix recording the prob for subjects from 1 to n 
      phi_a=exp(delta_at+tau_at*x)/denominator                        #P(Ci=at|xi)
      phi_c = 1 /denominator 
      pra=exp(alpha_ra+beta_ra*x)/(1+exp(alpha_ra+beta_ra*x))                     #P(Zi=1|xi) in the form of para
      
      pc1=(1/sigma)*exp(-(y-cova%*%t(lambda_co)-cova%*%t(logodd))^2/(sigma^2))                #p(Yi=1|Zi=1,Ci=co,xi) in the form of para
      pat=(1/sigma)*exp(-(y-cova%*%t(lambda_at))^2/(sigma^2))                                 #p(Yi=1|Ci=at,xi)
      pnt=(1/sigma)*exp(-(y-cova%*%t(lambda_nt))^2/(sigma^2))                                 #p(Yi=1|Ci=nt,xi)
      pc0=(1/sigma)*exp(-(y-cova%*%t(lambda_co))^2/(sigma^2))                               #p(Yi=1|Zi=0,Ci=co,xi)
      
      f1 <- sum(z*a*(log(pra)+w1*log(phi_c)+(1-w1)*log(phi_a)))
      f2 <- sum(z*a*(w1*log(pc1)+(1-w1)*log(pat)))
      f3 <- sum(z*(1-a)*(log(pra)+log(phi_n)+log(pnt)))
      f4 <- sum((1-z)*(1-a)*(log(1-pra)+w2*log(phi_c)+(1-w2)*log(phi_n)))
      f5 <- sum((1-z)*(1-a)*(w2*log(pc0)+(1-w2)*log(pnt)))
      f6 <- sum((1-z)*a*(log(1-pra)+log(phi_a)+log(pat)))
      obj<-  -(f1+f2+f3+f4+f5+f6)
    }
    
    p<- c(delta_nt,tau_nt,delta_at,tau_at,lambda_at,lambda_co,lambda_nt,logodd,sigma,alpha_ra,beta_ra)
    pnew<-optim(p,obj)$par
    distance=sqrt(sum(p-pnew)^2);
    ittime=ittime+1
    if ((distance<tol)|(ittime>=maxiter)){
      flag<-FALSE
    }else{
      delta_nt<-pnew[1]
      tau_nt<-pnew[2]
      delta_at<-pnew[3]
      tau_at<-pnew[4]
      lambda_at<-cbind(pnew[5],pnew[6])                                                                                                            #  update parameters for outcome model              
      lambda_co<-cbind(pnew[7],pnew[8]) 
      lambda_nt<-cbind(pnew[9],pnew[10])
      logodd<-cbind(pnew[11],pnew[12])
      sigma<-pnew[13]
      alpha_ra<-pnew[14]                                                                                                            #  update parameters for random assignment model
      beta_ra<-pnew[15] 
    }
  }
  cbind(ittime,distance,delta_nt,tau_nt,delta_at,tau_at,lambda_at,lambda_co,lambda_nt,logodd,sigma,alpha_ra,beta_ra)
}

#Constrained EM algorithm maximize over the constrained subspace lambda_at=lambda_co+logodd, lambda_nt=lambda_co
CEM<-function(data,delta_nt_0,tau_nt_0,delta_at_0,tau_at_0,lambda_co_0,logodd_0,sigma_0,alpha_ra_0,beta_ra_0,tol,maxiter){   
  y<-data[,1]                         # input the simulated data
  z<-data[,2]
  a<-data[,3]
  x<-data[,5]
  cova=cbind(1,x)
  n<-length(y)
  # insert the initial value of the EM algorithm
  
  delta_nt<-delta_nt_0               # initial parameters for compliance classes model
  tau_nt<-tau_nt_0
  delta_at<-delta_at_0
  tau_at<-tau_at_0
  lambda_co<-lambda_co_0
  logodd<-logodd_0
  lambda_nt<-lambda_co
  lambda_at<-lambda_co+logodd
  sigma<-sigma_0
  alpha_ra<-alpha_ra_0               # initial parameters for random assignment model
  beta_ra<-beta_ra_0
  
  flag<-T
  ittime=0
  
  while(flag){
    #E step using the old parameter to compute the expectation
    denominator=1+exp(delta_nt+tau_nt*x)+exp(delta_at+tau_at*x)     # the denominator matrix of P(Ci|Xi)
    phi_n=exp(delta_nt+tau_nt*x)/denominator                        #P(Ci=nt|xi) which is a matrix recording the prob for subjects from 1 to n 
    phi_a=exp(delta_at+tau_at*x)/denominator                        #P(Ci=at|xi)
    phi_c = 1 /denominator 
    w1 <- matrix(0,n,1)
    w2 <- matrix(0,n,1)
    pc1=(1/sigma)*exp(-(y-cova%*%t(lambda_co)-cova%*%t(logodd))^2/(sigma^2))                #p(Yi=1|Zi=1,Ci=co,xi)
    pat=(1/sigma)*exp(-(y-cova%*%t(lambda_at))^2/(sigma^2))                                 #p(Yi=1|Ci=at,xi)
    pnt=(1/sigma)*exp(-(y-cova%*%t(lambda_nt))^2/(sigma^2))                                 #p(Yi=1|Ci=nt,xi)
    pc0=(1/sigma)*exp(-(y-cova%*%t(lambda_co))^2/(sigma^2))                               #p(Yi=1|Zi=0,Ci=co,xi)
    w1<-phi_c*pc1/(phi_c*pc1+phi_a*pat)
    w2<-phi_c*pc0/(phi_c*pc0+phi_n*pnt)
    
    
    #M step 
    obj <- function (para){                                                   # define the expected likelihood function 
      delta_nt<-para[1]                                                         # parameters for compliance classes model
      tau_nt<-para[2] 
      delta_at<-para[3] 
      tau_at<-para[4] 
      lambda_co<-cbind(para[5],para[6])                                                        #  parameters for outcome model              
      logodd<-cbind(para[7],para[8])
      lambda_at<-lambda_co+logodd
      lambda_nt<-lambda_co
      sigma<-para[9]
      alpha_ra<-para[10]                                                        # parameters for random assignment model
      beta_ra<-para[11] 
      
      denominator=1+exp(delta_nt+tau_nt*x)+exp(delta_at+tau_at*x)      # the denominator matrix of P(Ci|Xi) of para
      phi_n=exp(delta_nt+tau_nt*x)/denominator                        #P(Ci=nt|xi) which is a matrix recording the prob for subjects from 1 to n 
      phi_a=exp(delta_at+tau_at*x)/denominator                        #P(Ci=at|xi)
      phi_c = 1 /denominator 
      pra=exp(alpha_ra+beta_ra*x)/(1+exp(alpha_ra+beta_ra*x))                     #P(Zi=1|xi) in the form of para
      
      pc1=(1/sigma)*exp(-(y-cova%*%t(lambda_co)-cova%*%t(logodd))^2/(sigma^2))                #p(Yi=1|Zi=1,Ci=co,xi) in the form of para
      pat=(1/sigma)*exp(-(y-cova%*%t(lambda_at))^2/(sigma^2))                                 #p(Yi=1|Ci=at,xi)
      pnt=(1/sigma)*exp(-(y-cova%*%t(lambda_nt))^2/(sigma^2))                                 #p(Yi=1|Ci=nt,xi)
      pc0=(1/sigma)*exp(-(y-cova%*%t(lambda_co))^2/(sigma^2))                               #p(Yi=1|Zi=0,Ci=co,xi)
      
      f1 <- sum(z*a*(log(pra)+w1*log(phi_c)+(1-w1)*log(phi_a)))
      f2 <- sum(z*a*(w1*log(pc1)+(1-w1)*log(pat)))
      f3 <- sum(z*(1-a)*(log(pra)+log(phi_n)+log(pnt)))
      f4 <- sum((1-z)*(1-a)*(log(1-pra)+w2*log(phi_c)+(1-w2)*log(phi_n)))
      f5 <- sum((1-z)*(1-a)*(w2*log(pc0)+(1-w2)*log(pnt)))
      f6 <- sum((1-z)*a*(log(1-pra)+log(phi_a)+log(pat)))
      obj<-  -(f1+f2+f3+f4+f5+f6)
    }
    
    p<- c(delta_nt,tau_nt,delta_at,tau_at,lambda_co,logodd,sigma,alpha_ra,beta_ra)
    pnew<-optim(p,obj)$par
    distance=sqrt(sum(p-pnew)^2);
    ittime=ittime+1
    if ((distance<tol)|(ittime>=maxiter)){
      flag<-FALSE
    }else{
      delta_nt<-pnew[1]
      tau_nt<-pnew[2]
      delta_at<-pnew[3]
      tau_at<-pnew[4]
      lambda_co<-cbind(pnew[5],pnew[6])                                                                                                            #  update parameters for outcome model              
      logodd<-cbind(pnew[7],pnew[8])
      lambda_at<-lambda_co+logodd
      lambda_nt<-lambda_co
      sigma<-pnew[9]
      alpha_ra<-pnew[10]                                                                                                            #  update parameters for random assignment model
      beta_ra<-pnew[11] 
    }
  }
  cbind(ittime,distance,delta_nt,tau_nt,delta_at,tau_at,lambda_at,lambda_co,lambda_nt,logodd,sigma,alpha_ra,beta_ra)
}
# observed data likelihood for EM estimator
#EMcoef[i,]=EM(data,delta_nt_0,tau_nt_0,delta_at_0,tau_at_0,lambda_at_0,lambda_co_0,lambda_nt_0,logodd_0,sigma_0,alpha_ra_0,beta_ra_0,10^(-6),maxiter)
observed.data.likelihood.EM<-function(para){
  ntp=rep(0,2)
  atp=rep(0,2)
  oc_atp=rep(0,2)
  oc_cop=rep(0,2)
  oc_ntp=rep(0,2)
  logoddp=rep(0,2)
  sigmap<-NA
  rap=rep(0,2)
  for(i in 1:2){
    ntp[i]=para[i]
  }
  for(i in 1:2){
    atp[i]=para[2+i]
  }
  for(i in 1:2){
    oc_atp[i]=para[4+i]
  }
  for(i in 1:2){
    oc_cop[i]=para[6+i]
  }
  for(i in 1:2){
    oc_ntp[i]=para[8+i]
  }
  for(i in 1:2){
    logoddp[i]=para[10+i]
  }
  logsigmap<-para[13]
  sigmap<-exp(logsigmap)
  for(i in 1:2){
    rap[i]=para[13+i]
  }
  
  #if(sigma <= 0 ) cat("Negative sigma=",sigma,"\n") 
  denominator=1+exp(XC%*%ntp)+exp(XC%*%atp)                # the expression of denominator matrix of P(Ci|Xi)
  phi_n=exp(XC%*%ntp)/denominator                          # P(Ci=nt|xi)
  phi_a=exp(XC%*%atp)/denominator                          # P(Ci=at|xi)
  phi_c=1/denominator
  pra=exp(XC%*%rap)/(1+exp(XC%*%rap))                        #P(Zi=1|xi)
  
  pc1=dnorm(y,mean=XC%*%oc_cop+XC%*%logoddp,sd=sigmap)        #p(Yi|Zi=1,Ci=co,xi)
  pat=dnorm(y,mean=XC%*%oc_atp,sd=sigmap)                       #p(Yi=1|Ci=at,xi)
  pnt=dnorm(y,mean=XC%*%oc_ntp,sd=sigmap)                          #p(Yi=1|Ci=nt,xi)
  pc0=dnorm(y,mean=XC%*%oc_cop,sd=sigmap)                            #p(Yi=1|Zi=0,Ci=co,xi)
  
  EM1 <- sum(z*a*log(pra))
  EM2 <- sum(z*a*log(phi_c*pc1 +phi_a*pat))
  EM3 <- sum(z*(1-a)*(log(pra)+log(phi_n)+log(pnt)))
  EM4 <- sum((1-z)*(1-a)*(log(1-pra)))
  EM5 <- sum((1-z)*(1-a)*log(phi_c*(pc0)+phi_n*(pnt)))
  EM6 <- sum((1-z)*a*(log(1-pra)+log(phi_a)+log(pat)))
  LEM=-(EM1+EM2+EM3+EM4+EM5+EM6)
  out<-LEM
}

observed.data.likelihood.CEM<-function(para){
  ntp=rep(0,2)
  atp=rep(0,2)
  oc_atp=rep(0,2)
  oc_cop=rep(0,2)
  oc_ntp=rep(0,2)
  logoddp=rep(0,2)
  sigmap<-0
  rap=rep(0,2)
  for(i in 1:2){
    ntp[i]=para[i]
  }
  for(i in 1:2){
    atp[i]=para[2+i]
  }
  for(i in 1:2){
    oc_atp[i]=para[4+i]
  }
  for(i in 1:2){
    oc_cop[i]=para[6+i]
  }
  for(i in 1:2){
    oc_ntp[i]=para[8+i]
  }
  for(i in 1:2){
    logoddp[i]=para[10+i]
  }
  logsigmap<-para[13]
  sigmap<-exp(logsigmap)
  for(i in 1:2){
    rap[i]=para[13+i]
  }
  oc_atp=oc_cop+logoddp
  oc_ntp=oc_cop
  # if (sigmap<0){
  #sigmap<-para[13]
  #}
  #if(sigma <= 0 ) cat("Negative sigma=",sigma,"\n") 
  denominator=1+exp(XC%*%ntp)+exp(XC%*%atp)                # the expression of denominator matrix of P(Ci|Xi)
  phi_n=exp(XC%*%ntp)/denominator                          # P(Ci=nt|xi)
  phi_a=exp(XC%*%atp)/denominator                          # P(Ci=at|xi)
  phi_c=1/denominator
  pra=exp(XC%*%rap)/(1+exp(XC%*%rap))                        #P(Zi=1|xi)
  
  pc1=dnorm(y,mean=XC%*%oc_cop+XC%*%logoddp,sd=sigmap)        #p(Yi|Zi=1,Ci=co,xi)
  pat=dnorm(y,mean=XC%*%oc_atp,sd=sigmap)                       #p(Yi=1|Ci=at,xi)
  pnt=dnorm(y,mean=XC%*%oc_ntp,sd=sigmap)                          #p(Yi=1|Ci=nt,xi)
  pc0=dnorm(y,mean=XC%*%oc_cop,sd=sigmap)                            #p(Yi=1|Zi=0,Ci=co,xi)
  
  EM1 <- sum(z*a*log(pra))
  EM2 <- sum(z*a*log(phi_c*pc1 +phi_a*pat))
  EM3 <- sum(z*(1-a)*(log(pra)+log(phi_n)+log(pnt)))
  EM4 <- sum((1-z)*(1-a)*(log(1-pra)))
  EM5 <- sum((1-z)*(1-a)*log(phi_c*(pc0)+phi_n*(pnt)))
  EM6 <- sum((1-z)*a*(log(1-pra)+log(phi_a)+log(pat)))
  LCEM=-(EM1+EM2+EM3+EM4+EM5+EM6)
  out<-LCEM
}


###simulation start!####


sims=1000 #number of iteration

datastore<-array(0,c(n,5,sims))
EMcoef=matrix(0,sims,17)
EMcoef.BFGS<-matrix(0,sims,15)
CEMcoef=matrix(0,sims,17)
CEMcoef.BFGS<-matrix(0,sims,15)
LRtest=matrix(0,sims,8)
Htest=matrix(0,sims,2)

set.seed(1)
for(i in 1:sims){

  #scenario 5
  n<-1000
  k0_co<- 0.8
  k1_co <- 1
  k0_at <- 0.8
  k1_at <- 1
  k0_nt <- 0.8
  k1_nt <- 1
  k0_de <- 0.8
  k1_de <- 1
  
  a0 <- 0 #changeable
  a1 <- 0 
  
  X <- rbinom(n,1,0.5)
  
  mu_Z <- 1/(1+exp(1-2*X))
  Z <- rbinom(n,1,mu_Z)
  
  p_at <- exp(-2.5+3.5*X)/(1+exp(-2.5+3.5*X)+exp(-2.5+3.5*X))
  p_nt <- exp(-2.5+3.5*X)/(1+exp(-2.5+3.5*X)+exp(-2.5+3.5*X))
  p_co <- 0.75/(1+exp(-2.5+3.5*X)+exp(-2.5+3.5*X))
  p_de <- 0.25/(1+exp(-2.5+3.5*X)+exp(-2.5+3.5*X)) #Note that p_de is 1/3 of p_co.
  
  C <- c()
  for(j in 1:1000){
    C<-c(C,sample(c("co","at","nt","de"),1,replace=TRUE,prob=c(p_at[j],p_nt[j],p_co[j],p_de[j])))
  }
  
  
  A<-rep(0,1000)
  A[C=="at"]<-1
  A[C=="co"&Z==1]<-1
  A[C=="de"&Z==0]<-1
  
  Y<-rep(-1,1000)
  for(k in 1:1000){
    
    if(C[k]=="co"){
      mu_y<- k0_co+k1_co*X[k]+(a0+a1*X[k])*A[k]
      Y[k] <-rnorm(1,mu_y,1)
      
    }else if(C[k]=="at"){
      mu_y<- k0_at+k1_at*X[k]
      Y[k] <- rnorm(1,mu_y,1)   
    }else if(C[k]=="nt"){
      mu_y<- k0_nt+k1_nt*X[k]
      Y[k] <- rnorm(1,mu_y,1)
    }else{
      mu_y<- k0_de+k1_de*X[k]
      Y[k] <- rnorm(1,mu_y,1) 
    }
  }
  
  
  data<-data.frame(col1=Y,col2=Z,col3=A,col4=C,col5=X)
  data[data$col4=="nt",4]<-0
  data[data$col4=="co",4]<-1
  data[data$col4=="at",4]<-2
  
  
  comp=matrix(1,n,1)
  y=data[,1]
  z=data[,2]
  a=data[,3]
  x=data[,5]

  ols.model=lm(y~x+a)
  cov.ols=vcov(ols.model)
  ivreg.model=ivreg(y~x+a|x+z)
  cov.ivreg=vcov(ivreg.model)
  Htest[i,1]=t(coef(ivreg.model)-coef(ols.model))%*%solve(cov.ivreg-cov.ols)%*%(coef(ivreg.model)-coef(ols.model))
  Htest[i,2]=pchisq(Htest[i,1], df = 1, lower.tail = FALSE) #DWHT

  cova=cbind(1,x)
  XC<-cova
  ## estimate the coefficients for P(Ci|Xi) 
  comp[(a*(1-z))==1]=2
  comp[((1-a))*z==1]=0
  mydata<-data.frame(comp,x)
  mydata$comp=as.factor(mydata$comp)
  mldata<-mlogit.data(mydata, varying=NULL, choice="comp", shape="wide")
  mlogit.model<- mlogit(comp~1|x, data = mldata, reflevel="1")
  delta_nt_0=coef(mlogit.model)[1]
  delta_at_0=coef(mlogit.model)[2]
  tau_nt_0=coef(mlogit.model)[3]
  tau_at_0=coef(mlogit.model)[4]
  ## estimate the coefficients for the random assignment 
  randomag.model<-glm(z~x,family=binomial("logit"))
  alpha_ra_0=coef(randomag.model)[1]
  beta_ra_0=coef(randomag.model)[2]
  ## estimate the coefficients for the outcome model
  wholedata=data.frame(y,z,a,comp,x)
  complier.data<-subset(wholedata,comp==1)
  nevertaker.data<-subset(wholedata,comp==0)
  alwaystaker.data<-subset(wholedata,comp==2)
  complier.model<-glm(y~x+a+I(a*x),family=gaussian,data=complier.data)
  nevertaker.model<-glm(y~x,family=gaussian,data=nevertaker.data)
  alwaystaker.model<-glm(y~x,family=gaussian,data=alwaystaker.data)
  ##summary(complier.model)
  ##summary(nevertaker.model)
  ##summary(alwaystaker.model)
  lambda_co_0=cbind(coef(complier.model)[1],coef(complier.model)[2])
  logodd_0=cbind(coef(complier.model)[3],coef(complier.model)[4])
  lambda_at_0=cbind(coef(alwaystaker.model)[1],coef(alwaystaker.model)[2])
  lambda_nt_0=cbind(coef(nevertaker.model)[1],coef(nevertaker.model)[2])
  sigma_0=sqrt(var(alwaystaker.data[,1]))
  
  EMcoef[i,]=EM(data,delta_nt_0,tau_nt_0,delta_at_0,tau_at_0,lambda_at_0,lambda_co_0,lambda_nt_0,logodd_0,sigma_0,alpha_ra_0,beta_ra_0,10^(-6),maxiter)
  EM.est.start<-EMcoef[i,3:17]
  EM.est.start[13]<-log(EM.est.start[13])
  EM.est <- optim(EM.est.start,observed.data.likelihood.EM,method = "BFGS")
  EMcoef.BFGS[i,]=EM.est$par
  LEM.BFGS=-EM.est$value
  
  CEMcoef[i,]=CEM(data,delta_nt_0,tau_nt_0,delta_at_0,tau_at_0,lambda_co_0,logodd_0,sigma_0,alpha_ra_0,beta_ra_0,10^(-6),maxiter)
  CEM.est.start<-CEMcoef[i,3:17]
  CEM.est.start[13]<-log(CEM.est.start[13])
  CEM.est<-optim(CEM.est.start,observed.data.likelihood.CEM,method = "BFGS")
  CEMcoef.BFGS[i,]=CEM.est$par
  LCEM.BFGS=-CEM.est$value
  
  ## assign the EM to compute the log likelihood 
  delta_nt_EM<-EMcoef[i,3]                                                # obtain the EM estiamte
  tau_nt_EM<-EMcoef[i,4]
  delta_at_EM<-EMcoef[i,5]
  tau_at_EM<-EMcoef[i,6] 
  lambda_at_EM<-EMcoef[i,7:8]
  lambda_co_EM<-EMcoef[i,9:10]
  lambda_nt_EM<-EMcoef[i,11:12] 
  logodd_EM<-EMcoef[i,13:14]
  sigma_EM<-EMcoef[i,15]
  alpha_ra_EM<-EMcoef[i,16]                                                                                                         
  beta_ra_EM<-EMcoef[i,17]

  
  denominator_EM=1+exp(delta_nt_EM+tau_nt_EM*x)+exp(delta_at_EM+tau_at_EM*x)    # the expression of denominator matrix of P(Ci|Xi) of EM estimate 
  phi_n_EM=exp(delta_nt_EM+tau_nt_EM*x)/denominator_EM                                   #P(Ci=nt|xi) of EM estimate
  phi_a_EM=exp(delta_at_EM+tau_at_EM*x)/denominator_EM                                  #P(Ci=at|xi) of EM estimate
  phi_c_EM= 1 - phi_a_EM - phi_n_EM
  pra_EM=exp(alpha_ra_EM+beta_ra_EM*x)/(1+exp(alpha_ra_EM+beta_ra_EM*x))              #P(Zi=1|xi) of EM estimate
  
  
  
  pc1_EM=dnorm(y,mean=cova%*%(lambda_co_EM)+cova%*%(logodd_EM),sd=sigma_EM)                        #p(Yi=1|Zi=1,Ci=co,xi) in the form of para
  pat_EM=dnorm(y,mean=cova%*%(lambda_at_EM),sd=sigma_EM)                                           #p(Yi=1|Ci=at,xi)
  pnt_EM=dnorm(y,mean=cova%*%(lambda_nt_EM),sd=sigma_EM)                                           #p(Yi=1|Ci=nt,xi)
  pc0_EM=dnorm(y,mean=cova%*%(lambda_co_EM),sd=sigma_EM)                                           #p(Yi=1|Zi=0,Ci=co,xi)
  EM1 <- sum(z*a*log(pra_EM))
  EM2 <- sum(z*a*log(phi_c_EM*pc1_EM+phi_a_EM*pat_EM))
  EM3 <- sum(z*(1-a)*(log(pra_EM)+log(phi_n_EM)+log(pnt_EM)))
  EM4 <- sum((1-z)*(1-a)*(log(1-pra_EM)))
  EM5 <- sum((1-z)*(1-a)*log(phi_c_EM*pc0_EM+phi_n_EM*pnt_EM))
  EM6 <- sum((1-z)*a*(log(1-pra_EM)+log(phi_a_EM)+log(pat_EM)))
  LEM.check=EM1+EM2+EM3+EM4+EM5+EM6 
  EM.assigned<-EMcoef[i,3:17]
  EM.assigned[13]<-log(EM.assigned[13])
  LEM=-observed.data.likelihood.EM(EM.assigned)
  
  ## assign the CEM to compute the log likelihood 
  delta_nt_CEM<-CEMcoef[i,3]                                                # obtain the EM estiamte
  tau_nt_CEM<-CEMcoef[i,4]
  delta_at_CEM<-CEMcoef[i,5]
  tau_at_CEM<-CEMcoef[i,6] 
  lambda_at_CEM<-CEMcoef[i,7:8]
  lambda_co_CEM<-CEMcoef[i,9:10]
  lambda_nt_CEM<-CEMcoef[i,11:12] 
  logodd_CEM<-CEMcoef[i,13:14]
  sigma_CEM<-CEMcoef[i,15]
  alpha_ra_CEM<-CEMcoef[i,16]                                                                                                         
  beta_ra_CEM<-CEMcoef[i,17]
  
  denominator_CEM=1+exp(delta_nt_CEM+tau_nt_CEM*x)+exp(delta_at_CEM+tau_at_CEM*x)    # the expression of denominator matrix of P(Ci|Xi) of EM estimate 
  phi_n_CEM=exp(delta_nt_CEM+tau_nt_CEM*x)/denominator_CEM                                   #P(Ci=nt|xi) of EM estimate
  phi_a_CEM=exp(delta_at_CEM+tau_at_CEM*x)/denominator_CEM                                  #P(Ci=at|xi) of EM estimate
  phi_c_CEM= 1 - phi_a_CEM - phi_n_CEM
  pra_CEM=exp(alpha_ra_CEM+beta_ra_CEM*x)/(1+exp(alpha_ra_CEM+beta_ra_CEM*x))              #P(Zi=1|xi) of EM estimate
  
  pc1_CEM=dnorm(y,mean=cova%*%(lambda_co_CEM)+cova%*%(logodd_CEM),sd=sigma_CEM)                        #p(Yi=1|Zi=1,Ci=co,xi) in the form of para
  pat_CEM=dnorm(y,mean=cova%*%(lambda_at_CEM),sd=sigma_CEM)                                           #p(Yi=1|Ci=at,xi)
  pnt_CEM=dnorm(y,mean=cova%*%(lambda_nt_CEM),sd=sigma_CEM)                                           #p(Yi=1|Ci=nt,xi)
  pc0_CEM=dnorm(y,mean=cova%*%(lambda_co_CEM),sd=sigma_CEM)                                           #p(Yi=1|Zi=0,Ci=co,xi)
  
  
  CEM1 <- sum(z*a*log(pra_CEM))
  CEM2 <- sum(z*a*log(phi_c_CEM*pc1_CEM+phi_a_CEM*pat_CEM))
  CEM3 <- sum(z*(1-a)*(log(pra_CEM)+log(phi_n_CEM)+log(pnt_CEM)))
  CEM4 <- sum((1-z)*(1-a)*(log(1-pra_CEM)))
  CEM5 <- sum((1-z)*(1-a)*log(phi_c_CEM*pc0_CEM+phi_n_CEM*pnt_CEM))
  CEM6 <- sum((1-z)*a*(log(1-pra_CEM)+log(phi_a_CEM)+log(pat_CEM)))
  LCEM.check=CEM1+CEM2+CEM3+CEM4+CEM5+CEM6 
  CEM.assigned<-CEMcoef[i,3:17]
  CEM.assigned[13]<-log(CEM.assigned[13])
  LCEM=-observed.data.likelihood.CEM(CEM.assigned)
  
  LRtest[i,1]=LEM
  LRtest[i,2]=LEM.BFGS
  LRtest[i,3]=LCEM
  LRtest[i,4]=LCEM.BFGS
  LRtest[i,5]=-2*(LCEM-LEM)
  LRtest[i,6]=pchisq(LRtest[i,5], df = 4, lower.tail = FALSE) 
  LRtest[i,7]=-2*(LCEM.BFGS-LEM.BFGS)
  LRtest[i,8]=pchisq(LRtest[i,7], df = 4, lower.tail = FALSE) 

}  



sum(Htest[,2]<0.05)/sims #the result of DWHT
sum(Htest[,2]<0.01)/sims #the result of DWHT
sum(LRtest[,6]<0.05)/sims
sum(LRtest[,6]<0.01)/sims
sum(LRtest[,8]<0.05)/sims #the result of CLRT
sum(LRtest[,8]<0.01)/sims #the result of CLRT

