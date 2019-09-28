#----------SampleSize.Poisson.R-------------------------------------------------------------------

# Version of Dez/2016

SampleSize.Poisson <-
function(alpha=0.05,power=0.9,M=1,D=0,RR=2,precision=0.000001,Tailed="upper")
{  

if(sum(Tailed==c("upper","lower","two"))==0){stop(" 'Tailed' must be chosen among 'upper', 'lower' or 'two'.",call. =FALSE)}
if(is.numeric(RR)==FALSE){stop("'RR' must be a vector of numbers.",call. =FALSE)}
if(Tailed=="upper"&min(RR)<1){stop("For 'Tailed=upper' RR must be >=1",call. =FALSE)}
if(Tailed=="lower"&max(RR)>1){stop("For 'Tailed=lower' RR must be <=1",call. =FALSE)}

tai<- Tailed

teste1<- 0
MinCases<- M
Late<- D
####### Tests to verify the validity of the chosen parameters
if(teste1==0){if(alpha>0.5 | alpha<(10^(-7))){teste1<- 1; out<- c("alpha must be a number in the (1e-7,0.5] interval")}}



if(is.numeric(power)==FALSE){stop("'power' must be a vector of numbers greater than or equal to 'alpha' and smaller than 1.",call. =FALSE)}
if(sum(power<alpha)>0|sum(power>=1)>0){stop("'power' must be a vector of numbers greater than or equal to 'alpha' and smaller than 1.",call. =FALSE)}

if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer in the range [1,100]")}
if(Late<0 & teste1==0){teste1<- 1; out<- c("Negative values for D does not make sense. Use D>=0.") }
if(M<1 & teste1==0){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }
if(Late<0 & teste1==0){teste1<- 1; out<- c("Negative values for D does not make sense. Use 0<=D<=T.") }

if(length(RR)>10){stop("RR must be a vector of length smaller than 11.",call. =FALSE)}

RR<- as.numeric(names(table(RR)))
power<- as.numeric(names(table(power)))




##### Auxiliary function to be used in the search for sample size for each configuration of RR and power
find_N<- function(RRt,powert)
{

#---------CODE TO CALCULATE POWER, SIGNAL TIME AND SURVEILLANCE TIME FOR GIVEN T
#######--------------------------------------------------------------------------
power<- powert
RR<- RRt
faux<-
function(L=30,D=0,M=1,RR=1,alpha=0.05){

# ------------------- INPUT VARIABLES ----------------------------------------------------------
# L = maximum length of surveillance, defined in terms of expected counts under H)
# RR = relative risk, RR=1 corresponds to H0
# M = The minimum number of cases for which a signal is allowed to occur
# D = Time < T for first look at the data, defined in terms of the expected counts under H0
# alpha = significance level
# If Tailed="upper" (default), then the threshold is given as an upper boundary (H0: R<=1), Tailed="lower" for lower boundaries (H0: R>=1), and Tailed="two" for two-tailed testing (H0: R=1).



####### Tests to verify the validity of the chosen parameters
T<- L
teste1<- 0
MinCases<- M
Late<- D

if(T<=0){teste1<- 1; out<- c("T must be > 0")}
if(teste1==0){if(alpha>0.5 | alpha<(10^(-7))){teste1<- 1; out<- c("alpha must be a number in the (1e-7,0.5] interval")}}
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer.")}


if(M<1 & teste1==0){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }

# If the parameters are incorrect in any sense, the code is interrupted and an error message is informed according to the possibilities above
#------------------------------------------------------------------------------------------------------------------------------------------
if(teste1==1){stop(out,call.=FALSE)}

## calculates the critical value

cv<- CV.Poisson(SampleSize=T,D,M,alpha,Tailed=tai)

## calculates power and signal time for relative risk equal to 2

if(length(cv)==1){
PT<- Performance.Poisson(SampleSize=T,D,M,cv,RRt,Tailed=tai)
                }else{

cvc<- cv[1,1]
cvl<- cv[1,2]
resc<- Performance.Poisson(SampleSize=T,D,M,cvc,RRt,Tailed=tai)
resl<- Performance.Poisson(T,D,M,cvl,RRt,Tailed=tai)
powerc<- resc[1]
powerl<- resl[1]
signaltimec<- resc[2]
signaltimel<- resl[2]
SurveillanceTimec<- resc[3]
SurveillanceTimel<- resl[3]

PT<- matrix(c(powerc,powerl,signaltimec,signaltimel,SurveillanceTimec,SurveillanceTimel),ncol=2,byrow=T)
rownames(PT)<- c("Power","Signal Time","Surveillance Time")
colnames(PT)<- c("Conservative","Liberal")
                     }



# Output assigned as a vector
# ---------------------------


out=list(cv,PT)
names(out)<- c("CV","Power.SignalTime")
return(out)

}

##-------------------------------------------------------------------------------
######---------------------------------------------------------------------------

if(Tailed=="two"|Tailed=="upper"){
T1<- max(qexp(alpha,1),D)
T_min<- min(seq(0.001,M,0.01)[1-ppois(M-1,seq(0.001,M,0.01))>=alpha])
                                 }else{
                                       T1<- max(qexp(1-alpha,1),D)
                                       T_min<- max(seq(0.001,M,0.01)[ppois(M,seq(0.001,M,0.01))>=1-alpha])
                                      }

T1<- max(T1,T_min)
if(T1<30){T2<- 30}else{T2<- T1+30}
result<- faux(T2, D, M, RR,alpha)
pow<- result$Power.SignalTime[1]
while(pow<power){T2<- T2+30;result<- faux(T2, D, M, RR,alpha);pow<- result$Power.SignalTime[1]}
Told<- T2
Tm<- (T1+T2)/2
result<- faux(Tm, D, M, RR,alpha)
pow<- result$Power.SignalTime[1]
cont<- 0
poder<- matrix(0,31,1)
lim<- log((T2-T1)/precision)/log(2)+1
while((pow<power|power+precision<pow)&cont<lim){

                                       if(pow>power){T2<- Tm;Told<- Tm}else{T1<- Tm}
                                       Tm<- (T1+T2)/2; cont<- cont+1;result<- faux(Tm, D, M, RR,alpha); pow<- result$Power.SignalTime[1]
                                       poder[cont,1]<- pow
                                       
                                      }

#SignTime<- result$Power.SignalTime[2]
#SurvTime<- result$Power.SignalTime[3]
CV<- result$CV

res<- list(Tm,CV,alpha,pow)
#names(res)<- c("SampleSize","Critical value", "Type I error probability","Power")
return(res)

}############ Closes find_N, the function that obtains the solution for each configuration of RR and power
#############




# Finding the exact solution for each RR and for each power

RRs<- RR[order(RR)]
powers<- power[order(power)]

res<- matrix(,length(power)*length(RR),6)

colnames(res)<- c("Target RR","Target power","Sample Size","Critical value","Type I Error prob.","Actual power")
acon<- 1
for(i in 1:length(RR)){
 for(j in 1:length(power)){
                          rh<- find_N(RRs[i],powers[j])
                          res[acon,]<- c(RRs[i],powers[j],as.numeric(rh))
                          acon<- acon+1
                          }
                      }


# Ploting results

if(length(RR)*length(power)>1){  

if(length(power)>1){
for(i in 1:length(RR)){if(i==1){RRleg<- paste("RR=",RRs[i]);plot(res[res[,1]==RRs[i],6],res[res[,1]==RRs[i],3],type="l",xlab="Power",ylab="Sample Size",ylim=c(0,max(res[,3])),xlim=c(min(power),max(power)))}else{
lines(res[res[,1]==RRs[i],6],res[res[,1]==RRs[i],3],col=i)
RRleg<- cbind(RRleg,paste("RR=",RRs[i]))
                   } 
                                                                                                                                }
if(length(RR)<11){legend("topleft",RRleg[1:length(RR)],lty=1,col=seq(1,length(RR)),bty="n")} 
                   }

if(length(power)==1){plot(res[,1],res[,3],type="l",xlab="RR",ylab="Sample Size")}
                              }



SampleSize_by_RR_Power<- res
return(SampleSize_by_RR_Power)   

}

### Example
#SampleSize.Poisson(alpha=0.05,power=c(0.9,0.95),M=1,D=0,RR=c(0.7,1.5),precision=0.000001,Tailed="two")


