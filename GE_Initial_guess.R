######################################### PFFCS Generation ###########################################
LL_Case_one=function(x)
{
  AA=length(Data1)
  sum_1=0
  sum_2=0
  for(i in 1:length(Data1))
  {
    sum_1=sum_1+(k*(R[counter_1[i]]+1)-1)*log(1- (1-exp(-x[2]*Data1[i]))^x[1])
    sum_2=sum_2+(x[1]-1)*log(1-exp(-x[2]*Data1[i]))-x[2]*Data1[i]
  }
  valueone=AA*log(x[1])+AA*log(x[2])+(sum_1+sum_2)
  return(valueone)
}

LL_Case_two=function(x)
{
  AAA=length(Data2)
  sum_11=0
  sum_22=0
  for(i in 1:length(Data2))
  {
    sum_11=sum_11+(k*(R[counter_2[i]]+1)-1)*log(1- (1-exp(-x[2]*Data2[i]))^x[1])
    sum_22=sum_22+(x[1]-1)*log(1-exp(-x[2]*Data2[i]))-x[2]*Data2[i]
  }
  valuetwo=AAA*log(x[1])+AAA*log(x[2])+(sum_11+sum_22)
  return(valuetwo)
}


n<-127
k<-2
m<-113
alpha=3
da <-read.csv("Cause_2.csv")

R=da[,3]
Data=da[,1]/3000
Indicator=da[,2]

Data1=Data[which(Indicator==1)]
Data2=Data[which(Indicator==3)]
counter_1=which(Indicator==1)
counter_2=which(Indicator==3)
count_1=length(which(Indicator==1))
count_2=length(which(Indicator==3))
####Estimation
M_pi=count_1/(count_1+count_2)
ayan_1=optim(c(1,1),LL_Case_one,control=list(fnscale=-1))$par
ayan_2=optim(c(1,1),LL_Case_two,control=list(fnscale=-1))$par
M_Alpha_1=ayan_1[1]
M_Lambda_1=ayan_1[2]
M_Alpha_2=ayan_2[1]
M_Lambda_2=ayan_2[2]
c(M_pi,M_Alpha_1,M_Lambda_1,M_Alpha_2,M_Lambda_2)