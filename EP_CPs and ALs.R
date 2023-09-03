##############Estimation #################
alpha_1<-4.75
lambda_1<-2.25
alpha_2<-1.75
lambda_2<-1.20
pi=0.40
k<-3
n<-60
m<-60
#R<-c(n-m,rep(0,m-1))
#R<-c(rep(0,m-1),n-m)
R=c(rep(0,m/2-1),(n-m),rep(0,m/2))
Nrep=1000

MLL_EP<-function(x)
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
    sum1=sum1+log(b)-(b+1)*log(1+Data1[i])+(a-1)*log(1-(1+Data1[i])^(-b))
  }
  for(i in 1:length(Data2))
  {
    sum2=sum2+log(d)-(d+1)*log(1+Data2[i])+(c-1)*log(1-(1+Data2[i])^(-d))
  }
  for(i in 1: length(Data))
  {
    sum3=sum3+(k*(R[i]+1)-1)*log(p*(1-(1-(1+Data[i])^(-b))^a)+(1-p)*(1-(1-(1+Data[i])^(-d))^c))
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
ind_asymp<-matrix(0,1,5)
length<-matrix(0,1,5)
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
  MLE<-optim(c(pi,alpha_1,lambda_1,alpha_2,lambda_2),MLL_EP,lower=c(0.01, 0.01,0.01,0.01,0.01), upper=c(0.99,5,5,5,5), method="L-BFGS-B", control=list(fnscale=-1),hessian = TRUE)
  M_pi=MLE$par[1]
  M_Alpha_1=MLE$par[2]
  M_Alpha_2=MLE$par[4]
  M_Lambda_1=MLE$par[3]
  M_Lambda_2=MLE$par[5]
  hupps=-MLE$hessian
  X=solve(hupps)
  sd_one=sqrt(X[1,1])
  sd_two=sqrt(X[2,2])
  sd_three=sqrt(X[3,3])
  
  
  sd_four=sqrt(X[4,4])
  sd_five=sqrt(X[5,5])
  #X=Fisher(M_pi,M_Alpha_1,M_Lambda_1,M_Alpha_2,M_Lambda_2)
  if((M_pi<0.6)&&(M_Alpha_1<5.5)&& (M_Alpha_2<5.5)&&(M_Lambda_1<5.5)&&(M_Lambda_2<5.5) && (M_Alpha_1>0.8)&& (M_Alpha_2>0.8)&&(M_Lambda_1>0.2)&&(M_Lambda_2>0.2)&&(sd_one>0)&&(sd_two>0)&&(sd_three>0)&&(sd_four>0)&&(sd_five>0)&& (is.finite(sd_one)))
  {
    MLE_pi[h]=M_pi
    MLE_Alpha_1[h]=M_Alpha_1
    MLE_Lambda_1[h]=M_Lambda_1
    MLE_Alpha_2[h]=M_Alpha_2
    MLE_Lambda_2[h]=M_Lambda_2
    MLE_pi_MSE[h]=(MLE_pi[h]-pi)^2
    MLE_Alpha1_MSE[h]=(MLE_Alpha_1[h]-alpha_1)^2
    MLE_Alpha2_MSE[h]=(MLE_Alpha_2[h]-alpha_2)^2
    MLE_Lambda_1_MSE[h]=(MLE_Lambda_1[h]-lambda_1)^2
    MLE_Lambda_2_MSE[h]=(MLE_Lambda_2[h]-lambda_2)^2
    
    
    ind_asymp[1,1]=ind_asymp[1,1]+check(pi,(M_pi-gam_3[2]*sd_one),(M_pi+gam_3[2]*sd_one))
    ind_asymp[1,2]=ind_asymp[1,2]+check(alpha_1,(M_Alpha_1-gam_3[2]*sd_two),(M_Alpha_1+gam_3[2]*sd_two))
    ind_asymp[1,3]=ind_asymp[1,3]+check(lambda_1,(M_Lambda_1-gam_3[2]*sd_three),(M_Lambda_1+gam_3[2]*sd_three))
    ind_asymp[1,4]=ind_asymp[1,4]+check(alpha_2,(M_Alpha_2-gam_3[2]*sd_four),(M_Alpha_2+gam_3[2]*sd_four))
    ind_asymp[1,5]=ind_asymp[1,5]+check(lambda_2,(M_Lambda_2-gam_3[2]*sd_five),(M_Lambda_2+gam_3[2]*sd_five))
    
    
    
    
    length[1,]=length[1,]+cbind(2*gam_3[2]*sd_one,2*gam_3[2]*sd_two,2*gam_3[2]*sd_three,2*gam_3[2]*sd_four,2*gam_3[2]*sd_five)
    
    
    
    
    
  }else
  {
    h=h-1
  }
  print(h)
  h=h+1
}
IND=ind_asymp/Nrep
AL=length/Nrep

c(mean(MLE_pi),mean(MLE_pi_MSE),mean(MLE_Alpha_1),mean(MLE_Alpha1_MSE),mean(MLE_Alpha_2),mean(MLE_Alpha2_MSE),mean(MLE_Lambda_1),mean(MLE_Lambda_1_MSE),mean(MLE_Lambda_2),mean(MLE_Lambda_2_MSE))

aa<- c(IND[1,1],AL[1,1],IND[1,2],AL[1,2],IND[1,3],AL[1,3],IND[1,4],AL[1,4],IND[1,5],AL[1,5])
round(aa, digits = 4)

