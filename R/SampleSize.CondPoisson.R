#----------SampleSize.Poisson.R-------------------------------------------------------------------

# Version of Dez/2016

SampleSize.CondPoisson<- function(cc,D=0,M=1,alpha=0.05,power=0.9,RR=2)
{                                

# ------------------- INPUT VARIABLE ----------------------------------------------------------
# cc = number of adverse events observed in the historical period
# D = Time < T for first look at the data, defined in terms of the expected counts under H0
# M = The minimum number of cases for which a signal is allowed to occur
# alpha= significance level
# power= a vector of target powers. This code looks for the minimum sample size in order to get power under the target relative risk RR
# RR= a vector of target relative risks 

## Tests
if(length(alpha)>1|length(cc)>1|length(D)>1|length(M)>1){stop("The inputs alpha, cc, D, and M cannot be vectors.",call. =FALSE)}
if(is.numeric(cc)!=TRUE){stop("'cc' must be a positive integer.",call. =FALSE)}
if(is.numeric(D)!=TRUE){stop("'D' must be a positive number.",call. =FALSE)}
if(is.numeric(M)!=TRUE){stop("'M' must be a positive integer.",call. =FALSE)}
if(is.numeric(alpha)!=TRUE){stop("'alpha' must be a number in the (0,0.05] range.",call. =FALSE)}

if(is.numeric(RR)==FALSE){stop("'RR' must be a vector of numbers each greater than 1.",call. =FALSE)}
if(sum(RR<=1)>0){stop("'RR' must be a vector of numbers each greater than 1.",call. =FALSE)}
if(is.numeric(power)==FALSE){stop("'power' must be a vector of numbers greater than or equal to 'alpha' and smaller than 1.",call. =FALSE)}
if(sum(power<alpha)>0|sum(power>=1)>0){stop("'power' must be a vector of numbers greater than or equal to 'alpha' and smaller than 1.",call. =FALSE)}

if(length(RR)>10){stop("RR must be a vector of length smaller than 11.",call. =FALSE)}

RR<- as.numeric(names(table(RR)))
power<- as.numeric(names(table(power)))

if(alpha<=0|alpha>=0.5){stop("'alpha' must be a number in the (0,0.05] range.",call. =FALSE)}
if(cc<=0|round(cc)!=cc){stop("'cc' must be a positive integer.",call. =FALSE)}
if(D<0){stop("'D' must be a positive number.",call. =FALSE)}
if(M<=0|round(M)!=M){stop("'M' must be a positive integer.",call. =FALSE)}


##### Auxiliary function to be used in the search for sample size for each configuration of RR and power
find_N<- function(RRt,powert)
{
RR<- RRt
power<- powert
mu2<- cc ; mu1<- 0 ; alphare<- 0 ; while(abs(alphare-alpha)>0.00001){mu0<- (mu2+mu1)/2 ; alphare<- 1-ppois(cc,mu0); if(alphare>alpha){mu2<- mu0}else{mu1<- mu0}}

powmax<- 1-ppois(cc,RR*mu0) ; powmax<- max(alpha,powmax-0.02)
if(powmax<power){
pow2<- 0
RRsol<- RR
while(pow2<power){RRsol<- RRsol+ 0.1; pow2<- 1-ppois(cc,RRsol*mu0); pow2<- max(alpha,pow2-0.02)} 
stop(c("For this cc it is not possible to reach the desired power under this target RR. For these parameters the maximum power is around ", round(powmax,2),". The desired power can be reached under RR= ",RRsol+0.5,"."),call. =FALSE)
}


## Redefining variables for cross-functions consistency

cc1<- cc
D1<- D
M1<- M
alpha1<- alpha
power1<- power
RR1<- RR
 


## Finding initial value for Tmed to be used in the bisection method

T2<- 0
pow<- 0

 while(pow<power){
                 T1<- T2
                 T2<- T1+0.5
                 cvaux<- CV.CondPoisson(Inference="liberal", StopType="Tal",T=T2,cc=cc1,D=D1,M=M1,alpha=alpha1)
                 res<- Performance.CondPoisson(Inference="liberal" , cv=as.numeric(cvaux[[2]]),StopType="Tal", T=T2,cc=cc1,D=D1,M=M1,RR=RR1) 
                 pow<- res[[1]]                               
                 }  

## Using initial value of Tmed for finding the sample size solution

epsilon<- 0.000001
pow<- 0
cont<- 1
scape<- ceiling( log((T2-T1)/epsilon)/log(2) )

while(abs(pow-power)>epsilon&cont<scape){
Tmed<- (T1+T2)/2
cvaux<- CV.CondPoisson(Inference="liberal",StopType="Tal",T=Tmed,cc=cc1,D=D1,M=M1,alpha=alpha1)
res<- Performance.CondPoisson(Inference="liberal",cv=as.numeric(cvaux[[2]]),StopType="Tal",T=Tmed,cc=cc1,D=D1,M=M1,RR=RR1)
pow<- res[[1]]
cont<- cont+1
if(pow>power){T2<- Tmed}else{T1<- Tmed}
                                        }

## Exact CV:

## cvaux<- CV.CondPoisson(Inference="exact",StopType="Tal",T=Tmed,cc=cc1,D=D1,M=M1,alpha=alpha1)



#### Finding K associated to Tal

####### Function to calculate cMaxSPRT
# ------------------------------------------------------------
cLLR<- function(k,cc,tal)
{
if(k/cc<=tal){return(0)}else{return(cc*log((cc*(1+tal)/(cc+k)))+k*log((k*(1+tal)/(tal*(cc+k)))))}
}

kt<- 0
cvt<- 0
while(cvt<as.numeric(cvaux[[2]])){kt<- kt+1; cvt<- cLLR(k=kt,cc,tal=Tmed)}  # finding the associated K for Tmed

#result<- list(kt,Tmed,as.numeric(cvaux[2]),pow)

result<- list(kt,as.numeric(cvaux[2]),alpha,pow)

names(result)<- c("K","CV","Type I error probability","Power")
return(result)


}############ Closes find_N, the function that obtains the solution for each configuration of RR and power
#############



# Finding the exact solution for each RR and for each power

RRs<- RR[order(RR)]
powers<- power[order(power)]

res<- matrix(,length(power)*length(RR),6)

colnames(res)<- c("Target RR","Target power","Sample Size (K)","Critical value","Type I Error prob.","Actual power")
acon<- 1
for(i in 1:length(RR)){
 for(j in 1:length(power)){
                          rh<- find_N(RRs[i],powers[j])
                          res[acon,]<- c(RRs[i],powers[j],as.numeric(rh))
                          acon<- acon+1
                          }
                      }


# Ploting the graph

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



}  # closes SampleSize.CondPoisson function
#############################################
#############################################



## Example for sample size calculation

#system.time(SS<- SampleSize.CondPoisson(cc=50,D=0,M=1,alpha=0.05,power=c(0.8,0.9),RR=c(2,3)))






