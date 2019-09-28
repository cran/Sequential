
#----------SampleSize.Binomial.R-------------------------------------------------------------------

# Version of Nov/2016

# -------------------------------------------------------------------------------------------
# Function produces critical value for the continuous Sequential Binomial MaxSPRT
# -------------------------------------------------------------------------------------------

SampleSize.Binomial<- function(RR,alpha=0.05,power=0.9,M=1,z="n",p="n",Tailed="upper"){

if(Tailed!="upper"){stop("For this version of the Sequential package, SampleSize.Binomial works only for 'Tailed=upper'.",call. =FALSE)}

MinCases<- M


if(p=="n"&z=="n"){stop("Please, at least z or p must be provided.",call. =FALSE)}

if( z!="n"){if(sum(is.numeric(z))!=1){stop("Symbols and texts are not applicable for 'z'. It must be a number greater than zero.",call. =FALSE)}}
if(z<=0){stop("'z' must be a number greater than zero.",call. =FALSE)}

if(p!="n"){
if(is.numeric(p)!=TRUE){stop("Symbols and texts are not applicable for 'p'. It must be a probability measure.",call. =FALSE)}
if(z!="n"&p!="n"){if(p!= 1/(1+z)){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed. .",call. =FALSE)}}
if(p<=0|p>=1){stop("p must be a number greater than zero and smaller than 1.",call. =FALSE)}
           }
if(p!="n"){z<- 1/p-1}

# alpha = desired alpha level
# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1
# z = matching ratio between exposed and unexposed cases  
# p = probability of having a case under the null hypothesis 
# RR is the relative risk
if(is.numeric(RR)==FALSE){stop("'RR' must be a vector of numbers each greater than 1.",call. =FALSE)}
if(sum(RR<=1)>0){stop("'RR' must be a vector of numbers each greater than 1.",call. =FALSE)}
if(MinCases<1){stop("'M' must be a positive integer.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("'M' must be a positive integer.",call. =FALSE)}
if(alpha<=0|alpha>0.5|is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(is.numeric(power)==FALSE){stop("'power' must be a vector of numbers greater than or equal to 'alpha' and smaller than 1.",call. =FALSE)}
if(sum(power<alpha)>0){stop("'power' must be a vector of numbers greater than or equal to 'alpha' and smaller than 1.",call. =FALSE)}
if(sum(power<alpha)>0){stop("'power' must be a vector of numbers greater than or equal to 'alpha' and smaller than 1.",call. =FALSE)}

if(is.numeric(MinCases)==FALSE|MinCases!=round(MinCases)|MinCases<1){stop("'M' must be a positive integer.",call. =FALSE)}

if(length(RR)>10){stop("RR must be a vector of length smaller than 11.",call. =FALSE)}

RR<- as.numeric(names(table(RR)))
power<- as.numeric(names(table(power)))

Nr<- 1
while(1-pbinom(Nr-1,Nr,1/(1+z))>alpha){Nr<- Nr+1}
N1<- max(Nr,M)

N2<- N1
pow<- 0
M<- as.integer(M)

# Finding a bound for all N's

while(pow< max(power) ){
   N2<- N2+100
   N2<- as.integer(N2)
   cv<- CV.Binomial(N=N2,alpha,M,z)[[1]] ; pow<- Performance.Binomial(N=N2,M,cv,z,p="n",min(RR) )[[1]]   
   
                 }

##### Auxiliary function
find_N<- function(RRt,powert)
{
N10<- N1
N20<- N2
tes<- 0
  while(N20-N10>1&tes==0){
   Nm<- round((N10+N20)/2)
   res1<- CV.Binomial(N=Nm,alpha,M,z)
   cv<- res1[[1]]
   res2<- Performance.Binomial(N=Nm,M,cv,z,p="n",RRt)   
   pow<- res2[[1]]
   if(pow==powert){tes<- 1;Ns<- Nm;pow1=pow}else{if(pow>powert){N20<- Nm;pow1<- pow;Ns<- Nm}else{N10<- Nm}}   
                       }
                                 
# Ns is the solution

error1<- res1[[2]]

out<- list(Ns,cv,error1,pow1)
names(out)<- c("Required_N","cv","Type_I_Error","Actual_power")
return(out)

} ###### Finish auxiliary function



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

} #end function SampleSize.Binomial


########################################################
## EXAMPLE.
## Note: remove the symbol "#" for running the lines below.
#system.time(result<- SampleSize.Binomial(RR=c(2,3,1.8,2.5),alpha=0.05,power=c(0.9,0.7,0.8),M=1,z="n",p=0.5))
#system.time(result<- SampleSize.Binomial(RR=c(2,5,1.8,2.5),alpha=0.05,power=c(0.9,0.7,0.8),M=1,z=1.2222222,p="n"))


