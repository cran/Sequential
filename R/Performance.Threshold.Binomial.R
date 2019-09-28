


#----------Performance.Threshold.Binomial.R-------------------------------------------------------------------

# Version of September/2019

# -------------------------------------------------------------------------
# Function produces alpha spending for user-defined thresholds 
# -------------------------------------------------------------------------

Performance.Threshold.Binomial<- function(N,CV.lower="n",CV.upper="n",z="n",p="n",GroupSizes="n",Tailed="upper",Statistic=c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis","Cases"),Delta="n",RR)
{


R<- RR
# N = maximum length of surveillance defined in terms of the total number of adverse events
# z = vector of matching ratios (between exposed and unexposed cases) for each group. If just a single number is given, then it will be used as a constant matching ratio for all groups. Otherwise, the dimension of z must coincide with the dimension of GroupSizes. 
# p = probability of having a case. If just a single number is given, then it will be used as a constant probability for all groups. Otherwise, the dimension of p must coincide with the dimension of GroupSizes.
# GroupSizes: Vector with the number of events (exposed+unexposed) between two looks at the data, i.e, irregular group sizes. Important: Must sums up N. For continuos sequential analysis, specify GroupSizes=1
# If Tailed="upper" (default), then the threshold is given as an upper tailed testing (H0: R<=RR0), Tailed="lower" for lower tailed (H0: R>=RR0), and Tailed="two" for two-tailed testing (H0: R=RR0). 
# RR= vector of relative risks for performance calculation

#### IMPORTANT OBSERVATION: Statistic can be different from "Cases" only for non-variable matching ratio. Otherwise, this function works only with "Statistic=Cases". 

Groups<- GroupSizes



####
# Some checks for input parameters
####

if(length(Statistic)>1){stop("'Statistic' must be selected among MaxSPRT, Pocock, OBrien-Fleming, Wang-Tsiatis and Cases.",call. =FALSE)}
if(sum(Statistic==c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis","Cases"))==0){stop("'Statistic' must be selected among MaxSPRT, Pocock, OBrien-Fleming, Wang-Tsiatis and Cases.",call. =FALSE)}

if(Statistic=="Cases"){
   cases.lower<- CV.lower; cases.upper<- CV.upper
   cvs.lower<- "n"; cvs.upper<- "n"
                      }else{
   cvs.lower<- CV.lower; cvs.upper<- CV.upper
   cases.lower<- "n"; cases.upper<- "n" 

                           }


if((is.numeric(N)==FALSE)){stop("'N' must be a positive integer.",call. =FALSE)}
if(N<=0|round(N)!=N){stop("'N' must be a positive integer.",call. =FALSE)}

if(length(Groups)==1){if(Groups=="n"){Groups<- rep(1,N)}}

if(length(p)==1&length(z)==1){if(p=="n"&z=="n"){stop("Please, at least one of the inputs, z or p, must be provided.",call. =FALSE)}}

if(length(z)>1){if(sum(is.numeric(z))!=1){stop("Symbols and texts are not applicable for 'z'. It must be a vector of positive numbers.",call. =FALSE)}}
if(length(z)==1){if(z!="n"&is.numeric(z)!=TRUE){stop("Symbols and texts are not applicable for 'z'. It must be a vector of positive numbers.",call. =FALSE)}}
if(sum(z<=0)>0){stop("'z' must be a vector of positive numbers.",call. =FALSE)}

if(length(p)==1){
if(p!="n"){
if(sum(is.numeric(p))!=1){stop("Symbols and texts are not applicable for 'p'. It must contain only probability measures.",call. =FALSE)}
if(sum(p<=0)>0|sum(p>=1)>0){stop("Each entry of p must be a number greater than zero and smaller than 1.",call. =FALSE)}

if(length(z)==1){
    if(z!="n"){
if(length(z)!=length(p)){stop("'z' and 'p' are vectors that must have the same dimension.",call. =FALSE)}
if(sum(p== 1/(1+z))!=length(p)&z!="n"){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed. .",call. =FALSE)}
              }
                }
if(length(z)>1){
    
if(length(z)!=length(p)){stop("'z' and 'p' are vectors that must have the same dimension.",call. =FALSE)}
if(sum(p== 1/(1+z))!=length(p)&z!="n"){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed. .",call. =FALSE)}
              
                }

          }
                }

if(length(p)>1){
if(sum(is.numeric(p))!=1){stop("Symbols and texts are not applicable for 'p'. It must contain only probability measures.",call. =FALSE)}
if(sum(p<=0)>0|sum(p>=1)>0){stop("Each entry of p must be a number greater than zero and smaller than 1.",call. =FALSE)}
if(length(z)==1){
    if(z!="n"){
if(length(z)!=length(p)){stop("'z' and 'p' are vectors that must have the same dimension.",call. =FALSE)}
if(sum(p== 1/(1+z))!=length(p)&z!="n"){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed. .",call. =FALSE)}
              }
                }
if(length(z)>1){
   
if(length(z)!=length(p)){stop("'z' and 'p' are vectors that must have the same dimension.",call. =FALSE)}
if(sum(p== 1/(1+z))!=length(p)){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed. .",call. =FALSE)}
              
                }
                }


if(length(z)==1){if(z=="n"){z<- 1/p-1}}

if(length(z)>1){if(length(z)!=length(GroupSizes)){stop("If the dimension of 'p', or equivalently of 'z', is greater than 1, then it must coincide with the dimension of 'GroupSizes'.",call. =FALSE)}}




if(length(Groups)==1){
if(is.numeric(Groups)==FALSE){stop("'Groups' must be an integer smaller than or equal to 'N'.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'Groups' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'Groups' must be a positive integer smaller than or equal to 'N'.",call. =FALSE)}

if(Groups==0){stop("'Groups' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(N/Groups!=round(N/Groups)){stop("The maximum length of surveillance, 'N', must be a multiple of 'Groups'.",call. =FALSE)}
if(Groups>N){stop("The maximum length of surveillance, 'N', must be a multiple of 'Groups'.",call.=FALSE)}
if(N/Groups==round(N/Groups)){Groups<- rep(Groups,N/Groups);GroupSizes<- Groups}
if(length(cases.upper)==1){if(is.numeric(cases.upper)==T){cases.upper<- rep(cases.upper,length(Groups))}}
if(length(cases.lower)==1){if(is.numeric(cases.lower)==T){cases.lower<- rep(cases.lower,length(Groups))}}
if(length(cvs.lower)==1){if(is.numeric(cvs.lower)==T){cvs.lower<- rep(cvs.lower,length(Groups))}}
if(length(cvs.upper)==1){if(is.numeric(cvs.upper)==T){cvs.upper<- rep(cvs.upper,length(Groups))}}
}

if(length(Groups)>1){
if(sum(is.numeric(Groups))==0){stop("'Groups' must be a vector of positive integers.",call. =FALSE)}else{
if(is.numeric(N)==FALSE){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}
if(N!=round(N)){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'Groups' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'Groups' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(sum(Groups)!=N){stop("'Groups' must sum up equal to 'N'.",call. =FALSE)}
}
}

if(max(z)==min(z)){
if(length(Statistic)>1){stop("'Statistic' must be selected among MaxSPRT, Pocock, OBrien-Fleming, Wang-Tsiatis and Cases.",call. =FALSE)}
if(sum(Statistic==c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis","Cases"))==0){stop("'Statistic' must be selected among MaxSPRT, Pocock, OBrien-Fleming, Wang-Tsiatis and Cases.",call. =FALSE)}
                  }else{Statistic<- "n"}

if(Statistic=="Wang-Tsiatis"&Delta=="n"){stop("Use a positive numeric value for 'Delta'",call. =FALSE)}
if(Statistic=="Wang-Tsiatis"&Delta<=0){stop("Use a positive numeric value for 'Delta'",call. =FALSE)}
if(Statistic=="Wang-Tsiatis"&Delta>0.5){stop("Use a positive numeric value smaller than or equal to 0.5 for 'Delta'",call. =FALSE)}

if(max(z)!=min(z)){           
if(length(cvs.lower)>1|length(cvs.upper)>1){    
 stop("For variable 'p', or equivalently variable 'z', the solution cannot be shown in a test statistic scale.",call. =FALSE)
                                           }else{
if(cvs.lower!="n"|cvs.upper!="n"){stop("For variable 'p', or equivalently variable 'z', the solution cannot be shown in a test statistic scale.",call. =FALSE)}
                                                }
                  }



if(Tailed=="upper"){
if(length(cvs.lower)==1&length(cases.lower)==1){if(cvs.lower!="n"|cases.lower!="n"){stop("For 'Tailed=upper' it is not possible to use 'cvs.lower' or 'cases.lower'.",call. =FALSE)}}
if(length(cvs.lower)>1|length(cases.lower)>1){stop("For 'Tailed=upper' it is not possible to use 'cvs.lower' or 'cases.lower'.",call. =FALSE)}
if(length(cases.upper)!=length(GroupSizes)&length(cvs.upper)!=length(GroupSizes)){stop("Signaling thresholds 'cvs.upper', or 'cases.upper' if the cases scale is used, must have the same dimension of 'GroupSizes'.",call. =FALSE)}
if(length(cases.upper)==1&length(cvs.upper)==1){if(cases.upper!="n"&cvs.upper!="n"){stop("Only one of the threshold scales, 'cvs.upper' or 'cases.upper', can be specified.",call. =FALSE)}}
if(length(cases.upper)>1&length(cvs.upper)>1){stop("Only one of the threshold scales, 'cvs.upper' or 'cases.upper', can be specified.",call. =FALSE)}
if(length(cases.upper)==1&length(cvs.upper)>1){if(cases.upper!="n"){stop("Only one of the threshold scales, 'cvs.upper' or 'cases.upper', can be specified.",call. =FALSE)}}
if(length(cases.upper)>1&length(cvs.upper)==1){if(cvs.upper!="n"){stop("Only one of the threshold scales, 'cvs.upper' or 'cases.upper', can be specified.",call. =FALSE)}}
                   }

if(Tailed=="lower"){
if(length(cvs.upper)==1&length(cases.upper)==1){if(cvs.upper!="n"|cases.upper!="n"){stop("For 'Tailed=lower' it is not possible to use 'cvs.upper' or 'cases.upper'.",call. =FALSE)}}
if(length(cvs.upper)>1|length(cases.upper)>1){stop("For 'Tailed=lower' it is not possible to use 'cvs.upper' or 'cases.upper'.",call. =FALSE)}
if(length(cases.lower)!=length(GroupSizes)&length(cvs.lower)!=length(GroupSizes)){stop("Signaling thresholds 'cvs.lower', or 'cases.lower' if the cases scale is used, must have the same dimension of 'GroupSizes'.",call. =FALSE)}
if(length(cases.lower)==1&length(cvs.lower)==1){if(cases.lower!="n"&cvs.lower!="n"){stop("Only one of the threshold scales, 'cvs.lower' or 'cases.lower', can be specified.",call. =FALSE)}}
if(length(cases.lower)>1&length(cvs.lower)>1){stop("Only one of the threshold scales, 'cvs.lower' or 'cases.lower', can be specified.",call. =FALSE)}
if(length(cases.lower)==1&length(cvs.lower)>1){if(cases.lower!="n"){stop("Only one of the threshold scales, 'cvs.lower' or 'cases.lower', can be specified.",call. =FALSE)}}
if(length(cases.lower)>1&length(cvs.lower)==1){if(cvs.lower!="n"){stop("Only one of the threshold scales, 'cvs.lower' or 'cases.lower', can be specified.",call. =FALSE)}}
                   }

### Checking the threshold scale choice 
aux_thres<- rep(0,4)
for(jj in 1:length(cvs.lower)){if(cvs.lower[jj]!="n"){aux_thres[1]<- 1}}
for(jj in 1:length(cvs.upper)){if(cvs.upper[jj]!="n"){aux_thres[2]<- 1}}
for(jj in 1:length(cases.lower)){if(cases.lower[jj]!="n"){aux_thres[3]<- 1}}
for(jj in 1:length(cases.upper)){if(cases.upper[jj]!="n"){aux_thres[4]<- 1}}

if(Tailed=="two"){
  if(length(GroupSizes)==1){
if(length(cvs.lower)>1|length(cases.lower)>1|length(cvs.upper)>1|length(cases.upper)>1){stop("Signaling thresholds must have the same dimension of 'GroupSizes'.",call. =FALSE)}
if(cvs.lower=="n"&cases.lower=="n"||cvs.upper=="n"&cases.upper=="n"){stop("Lower and upper signaling thresholds must be specified.",call. =FALSE)}
if(cvs.lower!="n"&cases.lower!="n"||cvs.upper!="n"&cases.upper!="n"){stop("Only one of the threshold scales, 'cvs.lower' or 'cases.lower', can be specified.",call. =FALSE)}
                           }

if(length(GroupSizes)>1){
if(aux_thres[1]+aux_thres[3]==2|aux_thres[1]+aux_thres[4]==2|aux_thres[2]+aux_thres[3]==2|aux_thres[2]+aux_thres[4]==2){stop("Only one of the threshold scales, 'cvs.lower' or 'cases.lower', can be specified.",call. =FALSE)}
if(aux_thres[1]+aux_thres[2]!=2&aux_thres[3]+aux_thres[4]!=2){stop("Lower and upper signaling thresholds must be specified.",call. =FALSE)}
if(aux_thres[1]+aux_thres[2]==2){if(length(cvs.lower)!=length(GroupSizes)|length(cvs.upper)!=length(GroupSizes)|length(cvs.upper)!=length(cvs.lower)){stop("Lower and upper signaling thresholds must have the same dimension of 'GroupSizes'.",call. =FALSE)}}
if(aux_thres[3]+aux_thres[4]==2){if(length(cases.lower)!=length(GroupSizes)|length(cases.upper)!=length(GroupSizes)|length(cases.upper)!=length(cases.lower)){stop("Lower and upper signaling thresholds must have the same dimension of 'GroupSizes'.",call. =FALSE)}}

                        }
                 }

if(aux_thres[1]==1|aux_thres[2]==1){stat_sca<- 1}else{stat_sca<- 0}


####
#### FUNCTIONS FOR THE TEST STATISTIC SCALE


########################################### 
#----- THE MAXSPRT STATISTIC (Wald-type)              

## LLR with variable z 

LLR<- function(zh,cch,nh,Tailed)
{

# zh: vector of macthing ratios
# cch: vector of cases
# nh: vector of events
# Tailed: possibilities are "upper", "lower" and "two" sided tests.

if(max(zh)==min(zh)){# 1

cc<- sum(cch)
n<- sum(nh)
z<- zh[1]

if(Tailed=="upper"){#2
if(cc>n){x<- 0}else{
       if(cc==n){x = n*log(1+z)}else{
         if(z*cc/(n-cc)<=1){x=0}else{
	       x = cc*log(cc/n)+(n-cc)*log((n-cc)/n)-cc*log(1/(z+1))-(n-cc)*log(z/(z+1))
                                    }
                                  } 
                            }	
      	if(x==0){res<- NA}else{res<- x}
                   }# 2 close

if(Tailed=="lower"){#3
if(cc<0){x<- 0}else{
       if(cc==0){x = n*log(1+z)}else{
         if(z*cc/(n-cc)>=1){x=0}else{
	       x = cc*log(cc/n)+(n-cc)*log((n-cc)/n)-cc*log(1/(z+1))-(n-cc)*log(z/(z+1))
                                    }
                                  } 
                            }	
      	if(x==0){res<- NA}else{res<- x}
                   }#3 close

if(Tailed=="two"){#4
if(cc<0|cc>n){x<- 0}else{
       if(cc==0|cc==n){x = n*log(1+z)}else{
         
	       x = cc*log(cc/n)+(n-cc)*log((n-cc)/n)-cc*log(1/(z+1))-(n-cc)*log(z/(z+1))
                                    
                                    } 
                       }	
     if(is.numeric(x)==FALSE){res<- NA}else{if(x==0){res<- NA}else{res<- x}}
                   }#4 close


                    }else{ # from 1
ord<- length(zh)
lenZs<- length(zh)
counta<- 1
rrest<- 0
for(yy in 1:ord){rrest<- max(rrest,max(zh[1:yy]*cch[1:yy]/(nh-cch)[1:yy]))}

LR<- function(Rcand){return(prod( dbinom(cch,nh,1/(1+zh/Rcand) ) ) )}
if(rrest<Inf){Rsup<- rrest}else{Rsup<- 200} 
if(sum(cch==nh)==length(cch)){rrest<- Inf}else{
tesaux<- 0 ; Rinf<- 0
while(tesaux==0){Rs<- matrix(seq(Rinf,Rsup,min(0.01,Rsup)),,1); LRs<- apply(Rs,1,LR); rrest<- Rs[LRs==max(LRs)];if(rrest< counta*200|counta>20){tesaux<- 1}else{counta<- counta+1 ; Rinf<- Rsup; Rsup<- Rsup+200}}
                                                                         }
if((Tailed=="upper"&rrest<1) || (Tailed=="lower"&rrest>1) ){rrest<-1} 
res<- log(LR(rrest))-log(LR(1))
                          }# 1 close


return(res)
}# CLOSE LLR Wald function
###########################
###########################

#----- THE Pocock (1977) statistic
LLR2 <- function(zh,cch,nh,Tailed){
cc<- sum(cch)
n<- sum(nh)
z<- zh[1]
p0<- 1/(1+z)
if(Tailed=="upper"){
   if(cc/n>p0){   
   x<- 1/sqrt(n)*abs( (cc-n*p0)/sqrt(p0*(1-p0)) )
              }else{x<- 0}
                   }
if(Tailed=="lower"){
   if(cc/n<p0){   
   x<- 1/sqrt(n)*abs( (cc-n*p0)/sqrt(p0*(1-p0)) )
              }else{x<- 0}
                   }
if(Tailed=="two"){      
   x<- 1/sqrt(n)*abs( (cc-n*p0)/sqrt(p0*(1-p0)) )            
                 }
      	x
                       }


#----- THE OBrien and Fleming (1972) statistic
LLR3 <- function(zh,cch,nh,Tailed){
cc<- sum(cch)
n<- sum(nh)
z<- zh[1]
p0<- 1/(1+z)
if(Tailed=="upper"){
      if(cc/n>p0){
      x<-  1/sqrt(N) * abs( (cc-n*p0)/sqrt(p0*(1-p0)) )
                }else{x<- 0}
                   }

if(Tailed=="lower"){
      if(cc/n<p0){
      x<-  1/sqrt(N) * abs( (cc-n*p0)/sqrt(p0*(1-p0)) )
                }else{x<- 0}
                   }
if(Tailed=="two"){      
      x<-  1/sqrt(N) * abs( (cc-n*p0)/sqrt(p0*(1-p0)) )                
                 }

          	x
                       }

#---- Wang e Tsiatis(1987), which is also a Pocock type, but with threshold multiplied by ((i/N)^(Delta-0.5)). Values of Delta can be 0.1, 0.25, and 0.4, page 40 by Jennison and Turniball(2000)

LLR4 <- function(zh,cch,nh,Tailed){
cc<- sum(cch)
n<- sum(nh)
z<- zh[1]
 	p0<- 1/(1+z)

if(Tailed=="upper"){
     if(cc/n>p0){
      x<-  (n/N)^(0.5-Delta)*(1/sqrt(n))* abs( (cc-n*p0)/sqrt(p0*(1-p0)) )
                }else{x<- 0}
                   }
if(Tailed=="lower"){
     if(cc/n<p0){
      x<-  (n/N)^(0.5-Delta)*(1/sqrt(n))* abs( (cc-n*p0)/sqrt(p0*(1-p0)) )
                 }else{x<- 0}
                   }
if(Tailed=="two"){     
      x<-  (n/N)^(0.5-Delta)*(1/sqrt(n))* abs( (cc-n*p0)/sqrt(p0*(1-p0)) )
                 }

      	x
                       }


if(Statistic=="Pocock"){LLR<- LLR2}; if(Statistic=="OBrien-Fleming"){LLR<- LLR3}; if(Statistic=="Wang-Tsiatis"){LLR<- LLR4} 


#######################
####################### AUXILIAR FUNCTION FOR PERFORMANCE CALCULATION UNDER ONE-TAILED TESTING

perf<- function(RR,res=0)
{

if(Tailed=="upper"|Tailed=="two"){pp<- 1/(1+z/RR)}; if(Tailed=="lower"){pp<- 1/(1+1/(z*RR))}

# Auxiliar functions to run the binomial Markov Chain in a fast way:
func_aux2<- function(j,i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*dbinom(j-k,1,pp[i])))}
func_aux3<- function(i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*(1-pbinom(absorb[i]-k,1,pp[i]))))}
func_aux1<- function(i){ j<- matrix(seq(1,absorb[i]),ncol=1) ; return(apply(j,1,func_aux2,i))}

p<- matrix(0,N,N+2)    	

for(s in 1:absorb[1]){ p[1,s]=dbinom(s-1,1,pp[1])}		
p[1,absorb[1]+1]=1-pbinom(absorb[1]-1,1,pp[1])			
for(s in 1:absorb[1]){ p[1,s]=dbinom(s-1,1,pp[1])}		
p[1,absorb[1]+1]=1-pbinom(absorb[1]-1,1,pp[1])			

if(N>1){
i<- 1

while(i<N){
i<- i+1
       p[i,1:absorb[i]]<- func_aux1(i)
       p[i,absorb[i]+1]<- func_aux3(i) 	
           }

       }

if(res==0){
power<- 0 ; AUX<- rep(0,N)
for(i in 1:N){power<- power+p[i,absorb[i]+1]; AUX[i]<- p[i,absorb[i]+1]}
ETS<- sum(seq(1,N)*AUX)/power; ELS<- sum(seq(1,N)*AUX)+N*(1-power)

measures<- c(power,ETS,ELS)
return(measures)
          }else{
AlphaSpend<- rep(0,N)
AlphaSpend[1]<- p[1,absorb[1]+1]
for(i in 2:N){AlphaSpend[i]<- AlphaSpend[i-1]+p[i,absorb[i]+1]}
return(AlphaSpend[an])
               }
}


#######################
####################### AUXILIAR FUNCTION FOR PERFORMANCE CALCULATION UNDER TWO-TAILED TESTING

perf2<- function(RR,res=0)
{

# res=0 for returning performance measures, and res=1 for returning alpha spending.

if(Tailed=="upper"|Tailed=="two"){pp<- 1/(1+z/RR)}; if(Tailed=="lower"){pp<- 1/(1+RR/z)}

func_aux2<- function(j,i){ k<- seq(uc1[i-1]+1,uc2[i-1]+1); return(sum(p[i-1,k]*dbinom(j-k,1,pp[i])))}
func_aux3<- function(i){ k<- seq(uc1[i-1]+1,uc2[i-1]+1); return(sum(p[i-1,k]*(1-pbinom(absorb2[i]-k,1,pp[i]))))}
func_aux4<- function(i){ k<- seq(uc1[i-1]+1,uc2[i-1]+1); return(sum(p[i-1,k]*(pbinom(absorb1[i]-k+1,1,pp[i]))))}
func_aux1<- function(i){ j<- matrix(seq(absorb1[i]+2,absorb2[i]),ncol=1) ; return(apply(j,1,func_aux2,i))}

p<- matrix(0,N,N+2)    	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's


for(s in absorb1[1]:absorb2[1]){ p[1,s]=dbinom(s-1,1,pp[1])}
pabs2<- 1-pbinom(absorb2[1]-1,1,pp[1])
pabs1<- pbinom(absorb1[1],1,pp[1])

pabs2 =1-pbinom(absorb2[1]-1,1,pp[1])			# probability of rejecting H0 by the upper side
                             
pabs1 =pbinom(absorb1[1],1,pp[1])			# probability of rejecting H0 by the lower side
                             

for(s in (absorb1[1]+2):absorb2[1]){ p[1,s]=dbinom(s-1,1,pp[1])}
p[1,absorb2[1]+1]=pabs2+pabs1

if(N>1){
i<- 1

jj<- 1
while(i<N){

i<- i+1
      	
       pabs2<- func_aux3(i) # Calculates the diagonal absorbing states where H0 is rejected            

       pabs1<- func_aux4(i) # Calculates the diagonal absorbing states where H0 is rejected

p[i,(absorb1[i]+2):absorb2[i]]<- func_aux1(i) # Calculates the standard p[][] cell values
p[i,absorb2[i]+1]<- pabs2+pabs1
               
          } # end for i
       }

if(res==0){
power<- 0 ; AUX<- rep(0,N)
for(i in 1:N){power<- power+p[i,absorb2[i]+1]; AUX[i]<- p[i,absorb2[i]+1]}
ETS<- sum(seq(1,N)*AUX)/power; ELS<- sum(seq(1,N)*AUX)+N*(1-power)

measures<- c(power,ETS,ELS)
return(measures)
          }else{
AlphaSpend<- rep(0,N)
AlphaSpend[1]<- p[1,absorb2[1]+1]
for(i in 2:N){AlphaSpend[i]<- AlphaSpend[i-1]+p[i,absorb2[i]+1]}
return(AlphaSpend[an])
               }
}


##########################################################################
### ANOTHER AUXILIAR FUNCTION

inv_upper<- function(cvs)
{
cases_sc<- ii
cvsaux<- LLR(z[1:ii],rep(1,ii),rep(1,ii),Tailed)
if(cvsaux<cvs){cases_sc_old<- ii+1}else{
cases_sc_old<- cases_sc
while(cvsaux>=cvs){cases_sc<- cases_sc-1;  cvsaux<- LLR(z[1:ii],c(rep(1,cases_sc),rep(0,ii-cases_sc)),rep(1,ii),Tailed);if(cvsaux>=cvs){cases_sc_old<- cases_sc}}
                                       }
return(cases_sc_old)
}

inv_lower<- function(cvs)
{
cases_sc<- 0
cvsaux<- LLR(z[1:ii],rep(0,ii),rep(1,ii),Tailed)
if(cvsaux<cvs){cases_sc_old<- -1}else{
cases_sc_old<- cases_sc
while(cvsaux>=cvs){cases_sc<- cases_sc+1; cvsaux<- LLR(z[1:ii],c(rep(1,cases_sc),rep(0,ii-cases_sc)),rep(1,ii),Tailed);if(cvsaux>=cvs){cases_sc_old<- cases_sc}}
                                     }
return(cases_sc_old)
}


##########################################################################
############################################################################


#### FINDING THE ALPHA SPENDING

if(min(z)==max(z)){z<- rep(z,N)}else{z<- rep(z,GroupSizes)}

cvs<- rep(0,N)

if(Tailed=="upper"|Tailed=="two"){ps<- 1/(1+z)}; if(Tailed=="lower"){ps<- 1/(1+1/z)}

G<- length(GroupSizes) ; an<- GroupSizes%*%(upper.tri(matrix(0,G,G),diag=T)*1);an<- an[1,1:G]


######################### ONE-SIDED TEST
if(Tailed!="two"){  ## HERE IS THE POINT WHERE THE CODE FOR TWO-TAILED AND ONE-TAILED TESTS DIFFER

### Establishing the absorb values

absorb = rep(0,N)
uc= rep(0,N)
for(ii in 1:N){
if(Tailed=="upper"&stat_sca==0){
if(sum(ii==an)>0){if(is.na(cases.upper[an==ii])==F){absorb[ii]<- cases.upper[an==ii]}else{absorb[ii]<- ii+1}}else{absorb[ii]<- ii+1}
                               }
if(Tailed=="lower"&stat_sca==0){
if(sum(ii==an)>0){if(is.na(cases.lower[an==ii])==F){absorb[ii]<- cases.lower[an==ii]}else{absorb[ii]<- -1}}else{absorb[ii]<- -1}
absorb[ii]<- ii-absorb[ii]
                               }

if(Tailed=="upper"&stat_sca==1){
if(sum(ii==an)>0){if(is.na(cvs.upper[an==ii])==F){absorb[ii]<- inv_upper(cvs.upper[an==ii])}else{absorb[ii]<- ii+1}}else{absorb[ii]<- ii+1}
                               }
if(Tailed=="lower"&stat_sca==1){
if(sum(ii==an)>0){if(is.na(cases.lower[an==ii])==F){absorb[ii]<- inv_lower(cvs.lower[an==ii])}else{absorb[ii]<- -1}}else{absorb[ii]<- -1}
absorb[ii]<- ii-absorb[ii]
                               }

uc[ii]<- absorb[ii]-1
              }

AlphaSpend<- perf(1,1)

#### SAVING PERFORMANCE MEASURES 
R<- matrix(R[order(R)],length(R),1)
Measures<- apply(R,1,perf)
Measures<- cbind(R,matrix(as.numeric(Measures),ncol=3,byrow=T)) 
colnames(Measures)<- c("R","Power","ESignalTime","ELS")

              }  ## CLOSE if(Tailed!="two")



##### TWO-SIDED CASE

if(Tailed=="two"){

### Establishing the absorb values

absorb1 = rep(0,N) ; absorb2 = rep(0,N)
uc1= rep(0,N); uc2= rep(0,N)
for(ii in 1:N){

if(sum(ii==an)>0){
 if(stat_sca==0){
  if(is.na(cases.upper[an==ii])==F){absorb2[ii]<- cases.upper[an==ii]}else{absorb2[ii]<- ii+1}
  if(is.na(cases.lower[an==ii])==F){absorb1[ii]<- cases.lower[an==ii]}else{absorb1[ii]<- -1}
                }else{
  if(is.na(cvs.upper[an==ii])==F){absorb2[ii]<- inv_upper(cvs.upper[an==ii])}else{absorb2[ii]<- ii+1}
  if(is.na(cvs.lower[an==ii])==F){absorb1[ii]<- inv_lower(cvs.lower[an==ii])}else{absorb1[ii]<- -1}
                     }
                 }else{absorb2[ii]<- ii+1; absorb1[ii]<- -1}
         
uc2[ii]<- absorb2[ii]-1 ; uc1[ii]<- absorb1[ii]+1
              }

AlphaSpend<- perf2(1,1)

#### SAVING PERFORMANCE MEASURES 
R<- matrix(R[order(R)],length(R),1)
Measures<- apply(R,1,perf2)
Measures<- cbind(R,matrix(as.numeric(Measures),ncol=3,byrow=T)) 
colnames(Measures)<- c("RR","Power","ESignalTime","ESampleSize")

              }  ## CLOSE if(Tailed=="two")


res<- list(AlphaSpend,Measures); names(res)<- c("AlphaSpend","Performance")

return(res)

} #end function AlphaSpend.Binomial


#### EXAMPLE

#res<- Performance.Threshold.Binomial(N=50,CV.upper=c(12,25,35,45),z=c(1,1.5,2,1.3),GroupSizes=c(15,15,10,10),Tailed="upper",Statistic="Cases",RR=c(1.2,1.5,2))


