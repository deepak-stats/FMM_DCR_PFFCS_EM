
result<-function(n,k,m,R,pi,alpha_1, alpha_2,lambda_1,lambda_2, PQ, NRep)
{
  
  MLL<-function(x)
  {
    p<-x[1]
    a<-x[2]
    b<-x[3]
    c<-x[4]
    d<-x[5]
    
    m1<-length(Data1)
    m2<-length(Data2)
    sum1<-0
    sum2<-0
    sum3<-0
    for(i in 1:length(Data1))
    {
      sum1=sum1+log(b)-b*Data1[i]+(a-1)*log(1-exp(-b*Data1[i]))
    }
    for(i in 1:length(Data2))
    {
      sum2=sum2+log(d)-d*Data2[i]+(c-1)*log(1-exp(-d*Data2[i]))
    }
    
    for(i in 1: length(Data))
    {
      ss<- max(p*(1-(1-exp(-b*Data[i]))^(a))+(1-p)*(1-(1-exp(-d*Data[i]))^(c)), 0.0000000000001)
      #print(ss)
      sum3=sum3+(k*(R[i]+1)-1)*log(ss)
    }
    f<-m1*log(p)+m2*log(1-p)+m1*log(a)+m2*log(c)+sum1+sum2+sum3
  }
  
  ##########################################################################
  check<-function(x,a,b)
  {
    if((x<=b)&&(x>=a))
    {
      m=1
    }
    else
    {
      m=0
    }
    return(m)
  }
  diff<-function(p,q)
  {
    r=q-p
    return(r)
  }
  #############################################
  ind_asymp<-matrix(0,3,5)
  length<-matrix(0,3,5)
  gam_3<-cbind(1.6448,1.96,2.576)
  V<-numeric()
  U<-numeric()
  MLE_pi=rep(0,Nrep)
  MLE_Alpha_1=rep(0,Nrep)
  MLE_Alpha_2=rep(0,Nrep)
  MLE_Lambda_1=rep(0,Nrep)
  MLE_Lambda_2=rep(0,Nrep)
  MLE_pi_MSE=rep(0,Nrep)
  MLE_Alpha1_MSE=rep(0,Nrep)
  MLE_Alpha2_MSE=rep(0,Nrep)
  MLE_Lambda_1_MSE=rep(0,Nrep)
  MLE_Lambda_2_MSE=rep(0,Nrep)
  PthQ=rep(0,Nrep)
  RbiasPQ=rep(0,Nrep)
  MSE_PthQ = rep(0, Nrep)
  h=1
  while(h<Nrep+1)
  {
    W<-runif(m,0,1)
    
    for(i in 1:m)
    {
      sum<-0
      for(j in 1:i)
      {
        sum<-sum+R[m-j+1]
      }
      ss<-i+sum
      V[i]<-W[i]^(1/(ss))
    }
    
    for(i in 1:m)
    {
      pp<-1
      for(j in 1:i)
      {
        pp<-pp*V[m-j+1]
      }
      U[i]<-1-pp
    }
    
    Data<-rep(0,m)
    Indicator=rep(0,m)
    
    for(i in 1:m)
    {
      if(W[i]<=pi)
      {
        Data[i]= ((1-(1-(1-U[i])^(1/k))^(1/alpha_1))^(-1/lambda_1))-1
        Indicator[i]=1
      }else
      {
        Data[i]=((1-(1-(1-U[i])^(1/k))^(1/alpha_2))^(-1/lambda_2))-1
        Indicator[i]=2
      }
      
    }
    #print(Data)
    Data1=Data[which(Indicator==1)]
    Data2=Data[which(Indicator==2)]
    counter_1=which(Indicator==1)
    counter_2=which(Indicator==2)
    count_1=length(which(Indicator==1))
    count_2=length(which(Indicator==2))
    ####Estimation
    MLE<-optim(c(0.5,1.5,1.5,0.5,1),MLL,lower=c(0.1,0.01,0.01,0.01,0.01), upper=c(0.99,6,6,6,6), method="L-BFGS-B", control=list(fnscale=-1))$par
    M_pi=MLE[1]
    M_Alpha_1=MLE[2]
    M_Alpha_2=MLE[4]
    M_Lambda_1=MLE[3]
    M_Lambda_2=MLE[5]
    #print(MLE)
    #X=Fisher(M_pi,M_Alpha_1,M_Lambda_1,M_Alpha_2,M_Lambda_2)
    
    if((M_Alpha_1<6.5)&& (M_Alpha_2<6.5)&&(M_Lambda_1<6.5)&&(M_Lambda_2<6.5) && (M_Alpha_1>0.1)&& (M_Alpha_2>0.1)&&(M_Lambda_1>0.1)&&(M_Lambda_2>0.1))
    {
      MLE_pi[h]=M_pi
      MLE_Alpha_1[h]=M_Alpha_1
      MLE_Alpha_2[h]=M_Alpha_2
      MLE_Lambda_1[h]= M_Lambda_1
      MLE_Lambda_2[h]=M_Lambda_2
      MLE_pi_MSE[h]=(MLE_pi[h]-pi)^2
      MLE_Alpha1_MSE[h]=(MLE_Alpha_1[h]-alpha_1)^2
      MLE_Alpha2_MSE[h]=(MLE_Alpha_2[h]-alpha_2)^2
      MLE_Lambda_1_MSE[h]=(MLE_Lambda_1[h]-lambda_1)^2
      MLE_Lambda_2_MSE[h]=(MLE_Lambda_2[h]-lambda_2)^2
      
      
      PthQ[h]= MLE_pi[h]*(-(1/MLE_Lambda_1[h])*log(1-(1-(1-0.5)^(1/k))^(1/MLE_Alpha_1[h]))) + (1-MLE_pi[h])*(-(1/MLE_Lambda_2[h])*log(1-(1-(1-0.5)^(1/k))^(1/MLE_Alpha_2[h])))
      MSE_PthQ[h]= (PthQ[h]-PQ)^2 
      RbiasPQ[h] = (PthQ[h]-PQ)/PQ
      
      
    }else
    {
      h=h-1
    }
    print(h)
    h=h+1
  }
  #IND=ind_asymp/Nrep
  #AL=length/Nrep
  
  return(list(ML=c(mean(MLE_pi),mean(MLE_pi_MSE),mean(MLE_Alpha_1),mean(MLE_Alpha1_MSE),mean(MLE_Alpha_2),mean(MLE_Alpha2_MSE),mean(MLE_Lambda_1),mean(MLE_Lambda_1_MSE),mean(MLE_Lambda_2),mean(MLE_Lambda_2_MSE), mean(PthQ),mean(MSE_PthQ),mean(RbiasPQ))))
}
##############Estimation #################
alpha_1<-3.0
lambda_1<-1.5
alpha_2<-4.0
lambda_2<-1.75
pi=0.35
k<-3
n<-60
m<-60
PQ<-(((1-(1-(1-0.5)^(1/k))^(1/alpha_1))^(-1/lambda_1))-1)*pi+(((1-(1-(1-0.5)^(1/k))^(1/alpha_2))^(-1/lambda_2))-1)*(1-pi)
#R<-c(n-m,rep(0,m-1))
#R<-c(rep(0,m-1),n-m)
R=c(rep(0,m/2-1),(n-m),rep(0,m/2))
Nrep=1000
result(n,k,m,R,pi,alpha_1, alpha_2,lambda_1,lambda_2,PQ, NRep)