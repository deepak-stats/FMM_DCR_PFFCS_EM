rm(list=ls())
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
L_1=function(xx)
{
  x<-xx[1]
  y<-xx[2]
  sum_1=0
  sum_2=0
  sum_3=0
  for(i in 1:length(Data))
  {
    sum_1=sum_1+(k*(R[i]+1)-1)*(updated_pi*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1))/((updated_pi*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((1-updated_pi)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2)))
    
    sum_2=sum_2+log(1-(1-exp(-y*Data[i]))^x)*(k*(R[i]+1)-1)*(updated_pi*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1))/((updated_pi*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((1-updated_pi)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2)))
  }
  for(i in 1:length(Data1))
  {
    sum_3=sum_3+log(x)+log(y)-y*Data1[i]+(x-1)*log(1-exp(-y*Data1[i]))
  }
  return((length(Data1)+sum_1)*log(updated_pi)+sum_3+sum_2)
}

L_2=function(yy)
{
  x<-yy[1]
  y<-yy[2]
  sum_4=0
  sum_5=0
  sum_6=0
  for(i in 1:length(Data))
  {
    sum_4=sum_4+(k*(R[i]+1)-1)*((1-updated_pi)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))/((updated_pi*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((1-updated_pi)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2)))
    
    sum_5=sum_5+log(1-(1-exp(-y*Data[i]))^x)*(k*(R[i]+1)-1)*((1-updated_pi)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))/((updated_pi*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((1-updated_pi)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2)))
  }
  for(i in 1:length(Data2))
  {
    sum_6=sum_6+log(x)+log(y)-y*Data2[i]+(x-1)*log(1-exp(-y*Data2[i]))
  }
  return((length(Data2)+sum_4)*log(1-updated_pi)+sum_6+sum_5)
}
M_pi<-numeric()
M_alpha_1<-numeric()
M_alpha_2<-numeric()
M_lambda_1<-numeric()
M_lambda_2<-numeric()
MSE_pi<-numeric()
MSE_alpha_1<-numeric()
MSE_alpha_2<-numeric()
MSE_lambda_1<-numeric()
MSE_lambda_2<-numeric()
cc1<-0
for(ll in 1:1000)
{
  cc1<-cc1+1
V<-numeric()
U<-numeric()

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
    Data[i]= -(1/lambda_1)*log(1-(1-(1-U[i])^(1/k))^(1/alpha_1))
    Indicator[i]=1
  }else
  {
    Data[i]=-(1/lambda_2)*log(1-(1-(1-U[i])^(1/k))^(1/alpha_2))
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







########## EM Algorithm ########
updated_pi=pi
updated_alpha_1=alpha_1
updated_lambda_1=lambda_1
updated_alpha_2=alpha_2
updated_lambda_2=lambda_2
pi_mle=c(0)
alpha_1_mle=c(0)
lambda_1_mle=c(0)
alpha_2_mle=c(0)
lambda_2_mle=c(0)
count<-0
indicatorr<-0
while(indicatorr<1)
{
  count<-count+1
  if(count>5)
  {
    indicatorr<-1
  }
  old_pi<-updated_pi
  old_alpha_1<-updated_alpha_1
  old_lambda_1<-updated_lambda_1
  old_alpha_2<-updated_alpha_2
  old_lambda_2<-updated_lambda_2
  
  
  sum_6=0
  
  for(i in 1:length(Data))
  {
    sum_6=sum_6+(k*(R[i]+1)-1)*(updated_pi*(1-   (1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1  ))/((updated_pi*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((1-updated_pi)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2)))
    
  }
  
  updated_pi=(count_1+sum_6)/(k*n)
  updated_1<-optim(c(alpha_1,lambda_1),L_1,lower=2, upper=5, method="L-BFGS-B", control=list(fnscale=-1))$par
  updated_alpha_1<-updated_1[1]
  updated_lambda_1<-updated_1[2]
  updated_2<-optim(c(alpha_2,lambda_2),L_2,lower=1, upper=2.5, method="L-BFGS-B", control=list(fnscale=-1))$par
  updated_alpha_2<-updated_2[1]
  updated_lambda_2<-updated_2[2]
 
  
  if((updated_alpha_1>0.5)&&(updated_lambda_1>0.5)&&(updated_lambda_2>0.5)&&(updated_alpha_2>0.5)&&(updated_pi>0.2))
  {
    updated_alpha_1<-updated_alpha_1
    updated_lambda_1<-updated_lambda_1
    updated_lambda_2<-updated_lambda_2
    updated_alpha_2<-updated_alpha_2
    updated_pi=updated_pi
  }  
  pi_mle=c(old_pi,updated_pi)
  alpha_1_mle=c(old_alpha_1,updated_alpha_1)
  lambda_1_mle=c(old_lambda_1,updated_lambda_1)
  alpha_2_mle=c(old_alpha_2,updated_alpha_2)
  lambda_2_mle=c(old_lambda_2,updated_lambda_2)
  #print(c(count,updated_pi,updated_alpha_1,updated_lambda_1,updated_alpha_2,updated_lambda_2))
  if((abs(pi_mle[2]-pi_mle[1])<0.001)&&(abs(alpha_1_mle[2]-alpha_1_mle[1])<0.01)&&(abs(lambda_1_mle[2]-lambda_1_mle[1])<0.001)&&(abs(lambda_2_mle[2]-lambda_2_mle[1])<0.001)&&(abs(alpha_2_mle[2]-alpha_2_mle[1])<0.001))
  {
    indicatorr<-1
  }
}
print(cc1)
M_pi[ll]<-updated_pi
M_alpha_1[ll]<-updated_alpha_1
M_lambda_1[ll]<-updated_lambda_1
M_alpha_2[ll]<-updated_alpha_2
M_lambda_2[ll]<-updated_lambda_2

MSE_pi[ll]<-(updated_pi-pi)^2
MSE_alpha_1[ll]<-(updated_alpha_1-alpha_1)^2
MSE_alpha_2[ll]<-(updated_alpha_2-alpha_2)^2
MSE_lambda_1[ll]<-(updated_lambda_1-lambda_1)^2
MSE_lambda_2[ll]<-(updated_lambda_2-lambda_2)^2
}
aa<-c(mean(M_pi), mean(MSE_pi), mean(M_alpha_1), mean(MSE_alpha_1), mean(M_lambda_1), mean(MSE_lambda_1), mean(M_alpha_2), mean(MSE_alpha_2), mean(M_lambda_2), mean(MSE_lambda_2))
round(aa, digits = 4)
