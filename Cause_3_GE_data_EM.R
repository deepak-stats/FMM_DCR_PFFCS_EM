ptm <- proc.time()
#R<-c(n-m,rep(0,m-1))
#R<-c(rep(0,m-1),n-m)
#R=c(rep(0,m/2-1),(n-m),rep(0,m/2))
L_1=function(xx)
{
  x<-xx[1]
  y<-xx[2]
  sum_4=0
  sum_5=0
  sum_6=0
  for(i in 1:length(Data))
  {
    sum_4=sum_4+(k*(R[i]+1)-1)*((updated_pi1)*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1))/((updated_pi1*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))+ ((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3)))
    
    sum_5=sum_5+log(1-(1-exp(-y*Data[i]))^x)*(k*(R[i]+1)-1)*((updated_pi1)*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1))/((updated_pi1*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))+ ((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3)))
  }
  for(i in 1:length(Data3))
  {
    sum_6=sum_6+log(x)+log(y)-y*Data1[i]+(x-1)*log(1-exp(-y*Data1[i]))
  }
  return((length(Data1)+sum_4)*log(updated_pi1)+sum_6+sum_5)
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
    sum_4=sum_4+(k*(R[i]+1)-1)*((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))/((updated_pi1*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))+ ((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3)))
    
    sum_5=sum_5+log(1-(1-exp(-y*Data[i]))^x)*(k*(R[i]+1)-1)*((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))/((updated_pi1*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))+ ((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3)))
  }
  for(i in 1:length(Data2))
  {
    sum_6=sum_6+log(x)+log(y)-y*Data2[i]+(x-1)*log(1-exp(-y*Data2[i]))
  }
  return((length(Data2)+sum_4)*log(updated_pi2)+sum_6+sum_5)
}

L_3=function(zz)
{
  x<-zz[1]
  y<-zz[2]
  sum_4=0
  sum_5=0
  sum_6=0
  for(i in 1:length(Data))
  {
    sum_4=sum_4+(k*(R[i]+1)-1)*((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3))/((updated_pi1*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))+ ((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3)))
    
    sum_5=sum_5+log(1-(1-exp(-y*Data[i]))^x)*(k*(R[i]+1)-1)*((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3))/((updated_pi1*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))+ ((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3)))
  }
  for(i in 1:length(Data3))
  {
    sum_6=sum_6+log(x)+log(y)-y*Data3[i]+(x-1)*log(1-exp(-y*Data3[i]))
  }
  return((length(Data3)+sum_4)*log(1-updated_pi1-updated_pi2)+sum_6+sum_5)
}


n<-186
k<-2
m<-166
alpha_1<-0.62956
lambda_1<-4.11247
alpha_2<-0.71322
lambda_2<-4.24875
alpha_3<-1.04885
lambda_3<-4.58185 
pi1=0.34337
pi2=0.34939
da <-read.csv("Cause_3.csv")

R=da[,3]
Data=da[,1]/3000
Indicator=da[,2]

Data1=Data[which(Indicator==1)]
Data2=Data[which(Indicator==2)]
Data3=Data[which(Indicator==3)]
counter_1=which(Indicator==1)
counter_2=which(Indicator==2)
counter_3=which(Indicator==3)
count_1=length(which(Indicator==1))
count_2=length(which(Indicator==2))
count_3=length(which(Indicator==3))

########## EM Algorithm ########
updated_pi1=pi1
updated_pi2=pi2
updated_alpha_1=alpha_1
updated_lambda_1=lambda_1
updated_alpha_2=alpha_2
updated_lambda_2=lambda_2
updated_alpha_3=alpha_3
updated_lambda_3=lambda_3
pi1_mle=c(0)
pi2_mle=c(0)
alpha_1_mle=c(0)
lambda_1_mle=c(0)
alpha_2_mle=c(0)
lambda_2_mle=c(0)
alpha_3_mle=c(0)
lambda_3_mle=c(0)
count<-0
indicatorr<-0
while(indicatorr<1)
{
  count<-count+1
  if(count>200)
  {
    indicatorr<-1
  }
  old_pi1<-updated_pi1
  old_pi2<-updated_pi2
  old_alpha_1<-updated_alpha_1
  old_lambda_1<-updated_lambda_1
  old_alpha_2<-updated_alpha_2
  old_lambda_2<-updated_lambda_2
  old_alpha_3<-updated_alpha_3
  old_lambda_3<-updated_lambda_3
  
  
  sum_6=0
  sum_66=0
  
  for(i in 1:length(Data))
  {
    sum_6=sum_6+(k*(R[i]+1)-1)*(updated_pi1*(1-   (1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1  ))/((updated_pi1*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))+ ((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3)))
    
  }
  
  for(i in 1:length(Data))
  {
    sum_66=sum_66+(k*(R[i]+1)-1)*(updated_pi2*(1-   (1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2  ))/((updated_pi1*(1-(1-exp(-updated_lambda_1*Data[i]))^updated_alpha_1)) +((updated_pi2)*(1-(1-exp(-updated_lambda_2*Data[i]))^updated_alpha_2))+ ((1-updated_pi1-updated_pi2)*(1-(1-exp(-updated_lambda_3*Data[i]))^updated_alpha_3)))
    
  }
  
  updated_pi1=(count_1+sum_6)/(k*n)
  updated_pi2=(count_2+sum_66)/(k*n)
  updated_1<-optim(c(alpha_1,lambda_1),L_1,lower=0.01, upper=5, method="L-BFGS-B", control=list(fnscale=-1))$par
  updated_alpha_1<-updated_1[1]
  updated_lambda_1<-updated_1[2]
  updated_2<-optim(c(alpha_2,lambda_2),L_2,lower=0.01, upper=5, method="L-BFGS-B", control=list(fnscale=-1))$par
  updated_alpha_2<-updated_2[1]
  updated_lambda_2<-updated_2[2]
  updated_3<-optim(c(alpha_3,lambda_3),L_3,lower=0.01, upper=5, method="L-BFGS-B", control=list(fnscale=-1))$par
  updated_alpha_3<-updated_3[1]
  updated_lambda_3<-updated_3[2]
  
  
  
  if((updated_alpha_1>0 && updated_alpha_1 <=4.5)&&(updated_lambda_1>0 && updated_lambda_1<=4.5)&&(updated_lambda_2>0 && updated_lambda_2<=4.5)&&(updated_alpha_2>0 && updated_alpha_2<=4.5)&& (updated_alpha_3>0 && updated_alpha_3 <=4.5)&&(updated_lambda_3>0 && updated_lambda_3<=4.5)&&(updated_pi1>0)&&(updated_pi2>0))
  {
    updated_alpha_1<-updated_alpha_1
    updated_lambda_1<-updated_lambda_1
    updated_lambda_2<-updated_lambda_2
    updated_alpha_2<-updated_alpha_2
    updated_lambda_3<-updated_lambda_3
    updated_alpha_3<-updated_alpha_3
    updated_pi1=updated_pi1
    updated_pi2=updated_pi2
  }  
  pi1_mle=c(old_pi1,updated_pi1)
  pi2_mle=c(old_pi2,updated_pi2)
  alpha_1_mle=c(old_alpha_1,updated_alpha_1)
  lambda_1_mle=c(old_lambda_1,updated_lambda_1)
  alpha_2_mle=c(old_alpha_2,updated_alpha_2)
  lambda_2_mle=c(old_lambda_2,updated_lambda_2)
  alpha_3_mle=c(old_alpha_3,updated_alpha_3)
  lambda_3_mle=c(old_lambda_3,updated_lambda_3)
  #print(c(count,updated_pi,updated_alpha_1,updated_lambda_1,updated_alpha_2,updated_lambda_2))
  if((abs(pi1_mle[2]-pi1_mle[1])<0.001)&&(abs(pi2_mle[2]-pi2_mle[1])<0.001)&&(abs(alpha_1_mle[2]-alpha_1_mle[1])<0.001)&&(abs(lambda_1_mle[2]-lambda_1_mle[1])<0.001)&&(abs(lambda_2_mle[2]-lambda_2_mle[1])<0.001)&&(abs(alpha_2_mle[2]-alpha_2_mle[1])<0.001)&&(abs(lambda_3_mle[2]-lambda_3_mle[1])<0.001)&&(abs(alpha_3_mle[2]-alpha_3_mle[1])<0.001))
  {
    indicatorr<-1
  }
}
MLL<-function(p1,p2,a,b,c,d,e,f)
{
  m1<-length(Data1)
  m2<-length(Data2)
  m3<-length(Data3)
  sum1<-0
  sum2<-0
  sum3<-0
  sum4<-0
  for(i in 1:length(Data1))
  {
    sum1=sum1+log(b)-b*Data1[i]+(a-1)*log(1-exp(-b*Data1[i]))
  }
  for(i in 1:length(Data2))
  {
    sum2=sum2+log(d)-d*Data2[i]+(c-1)*log(1-exp(-d*Data2[i]))
  }
  for(i in 1:length(Data3))
  {
    sum3=sum3+log(f)-f*Data3[i]+(e-1)*log(1-exp(-f*Data3[i]))
  }
  for(i in 1: length(Data))
  {
    sum4=sum4+(k*(R[i]+1)-1)*log(p1*(1-(1-exp(-b*Data[i]))^(a))+(p2)*(1-(1-exp(-d*Data[i]))^(c))+(1-p1-p2)*(1-(1-exp(-f*Data[i]))^(e)))
  }
  f<-m1*log(p1)+m2*log(p2)+ m2*log(1-p1-p2)+m1*log(a)+m2*log(c)+ m2*log(e)+sum1+sum2+sum3+sum4
}
aa<-c(updated_pi1, updated_pi2, updated_alpha_1, updated_lambda_1, updated_alpha_2, updated_lambda_2, updated_alpha_3, updated_lambda_3)
round(aa, digits = 4)
AIC<-10-2*(MLL(updated_pi1, updated_pi2,updated_alpha_1, updated_lambda_1,updated_alpha_2, updated_lambda_2,updated_alpha_3, updated_lambda_3))
AIC
proc.time() - ptm