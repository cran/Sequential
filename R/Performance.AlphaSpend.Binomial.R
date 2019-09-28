


#----------Performance.AlphaSpend.Binomial.R-------------------------------------------------------------------

# Version of September/2019

# -------------------------------------------------------------------------
# Function produces performance and threshold for a user-defined alpha spending
# -------------------------------------------------------------------------

Performance.AlphaSpend.Binomial<- function(N,alpha,AlphaSpend="n",z="n",p="n",GroupSizes="n",Tailed="upper",rho="n",gamma="n",RR=2,Statistic=c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis"),Delta="n")
{

R<- RR
# N = maximum length of surveillance defined in terms of the total number of adverse events
# alpha = desired alpha level
# AlphaSpend: cumulative type I error probability spending. It can also be choosen among four classical alpha spending shapes, which are indicated with numbers 1 to 4.
# z = vector of matching ratios (between exposed and unexposed cases) for each group. If just a single number is given, then it will be used as a constant matching ratio for all groups. Otherwise, the dimension of z must coincide with the dimension of GroupSizes. 
# p = probability of having a case. If just a single number is given, then it will be used as a constant probability for all groups. Otherwise, the dimension of p must coincide with the dimension of GroupSizes.
# GroupSizes: Vector with the number of events (exposed+unexposed) between two looks at the data, i.e, irregular group sizes. Important: Must sums up N. For continuos sequential analysis, specify GroupSizes=1
# If Tailed="upper" (default), then the threshold is given as an upper boundary (H0: R<=1), Tailed="lower" for lower boundaries (H0: R>=1), and Tailed="two" for two-tailed testing ((H0: R=1). 
# Statistic= choose the test statistic scale for the threshold
# RR= vector of relative risks for calculation of power, expected time to signal, and expected length of surveillance.


Groups<- GroupSizes
AlphaSpend2<- AlphaSpend


############################
### Four classical alpha spending shapes

# AlphaSpend=1: power-type t^rho, Kim and DeMets (1987a), and Jennison & Turnball(1989,1990), noted that this function produces Pocock and O'Brien & Fleming tests. Values of 1,2 and 3 for mimicking Wang & Tsiatis(1987) test.

alpha_spendT1<- function(N,alpha,rho){
x<- seq(1/N,by=1/N,1)
sum_sa<- alpha*(x^rho)
return(sum_sa[an])
}


# AlphaSpend= 2: Gaussian-type, Lan & DeMets(1983) suggested this one for mimicking OBrien and Fleming test.  

alpha_spendT2<- function(N,alpha){
x<- seq(1/N,by=1/N,1)
za<- qnorm(1-alpha/2)
sum_sa<- 2-2*pnorm(za/sqrt(x))
return(sum_sa[an])
}


# AlphaSpend= 3: LogExp-type, Lan & DeMets(1983) indicated this one for mimicking Pocock's test.

alpha_spendT3<- function(N,alpha){
x<- seq(1/N,by=1/N,1)
sum_sa<- alpha*log(1+(exp(1)-1)*x)
return(sum_sa[an])
}


# AlphaSpend= 4: Gamma-type, Hwang, Shih & DeCani(1990). 

alpha_spendT4<- function(N,alpha,gamma){
x<- seq(1/N,by=1/N,1)
if(gamma==0){sum_sa<- alpha*x}else{
sum_sa<- alpha*(1-exp(-gamma*x))/(1-exp(-gamma))
                                  }
return(sum_sa[an])
}

####
# Some checks for input parameters
####

if(length(z)>1){if(length(z)!=length(GroupSizes)&length(GroupSizes)>1&min(z)!=max(z)){stop("For variable 'p', or equivalently for variable 'z', if the dimension of 'GroupSizes' is greater than 1, then these vectors must have the same dimension.",call. =FALSE)}}


alpha1<- alpha
if(length(Groups)==1){
if(is.numeric(Groups)==FALSE){stop("'GroupSizes' must be an integer smaller than or equal to 'N'.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}

if(Groups==0){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}
if(N/Groups!=round(N/Groups)){stop("The maximum length of surveillance, 'N', must be a multiple of 'GroupSizes'.",call. =FALSE)}
if(Groups>N){stop("The maximum length of surveillance, 'N', must be a multiple of 'GroupSizes'.",call.=FALSE)}
}

if(length(Groups)>1){
if(sum(is.numeric(Groups))==0){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}else{
if(is.numeric(N)==FALSE){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}
if(N!=round(N)){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}
if(sum(Groups)!=N){stop("'GroupSizes' must sum up equal to 'N'.",call. =FALSE)}
}
}



### IMPORTANT VARIABLE
if(length(Groups)==1){if(Groups=="n"|Groups==1){Groups<- rep(1,N)}else{Groups<- rep(Groups,N/Groups)}}
G<- length(Groups) ; an<- Groups%*%(upper.tri(matrix(0,G,G),diag=T)*1);an<- an[1,1:G]


if((is.numeric(alpha)==FALSE)){stop("'alpha' must be a number in the '(0,0.5]' interval.",call. =FALSE)}
if((alpha<=0|alpha>0.5)){stop("'alpha' must be a number in the '(0,0.5]' interval.",call. =FALSE)}

if((is.numeric(N)==FALSE)){stop("'N' must be a positive integer.",call. =FALSE)}
if(N<=0|round(N)!=N){stop("'N' must be a positive integer.",call. =FALSE)}


if(length(AlphaSpend2)==1){ 
                          if( AlphaSpend2=="n" ){stop(" 'AlphaSpend' must be a vector containing a non-decreasing and positive sequence of numbers with maximum at alpha, with same length as 'GroupSizes', or a single integer among 1 to 4.",call. =FALSE)}
                          if(sum(AlphaSpend2)!=alpha&sum(AlphaSpend2==c(1,2,3,4))==0){stop(" 'AlphaSpend' must be a vector containing a non-decreasing and positive sequence of numbers with maximum at alpha, with same length as 'GroupSizes', or a single integer among 1 to 4.",call. =FALSE)}
                          if(AlphaSpend2==alpha){if(length(GroupSizes)>1){stop(" 'AlphaSpend' and 'GroupSizes' must have the same length.",call. =FALSE)};if(GroupSizes=="n"){GroupSizes<- N}}
                          if( AlphaSpend2==1 ){
                                              if(length(rho)>1){stop(" 'rho' must be a single number greater than zero and smaller than 5.",call. =FALSE)}               
                                              if(rho=="n"){stop(" Please, choose 'rho' between 0 and 5.",call. =FALSE)} 
                                              if(sum(is.numeric(rho))!=1){stop("Symbols and texts are not applicable for 'rho'. It must be a number greater than zero and smaller than 5.",call. =FALSE)}
                                              if(rho<=0|rho>5){stop(" 'rho' must be a single number greater than zero and smaller than 5.",call. =FALSE)}
                                              AlphaSpend<- alpha_spendT1(N,alpha,rho)
                                             }
                          if( AlphaSpend2==2 ){AlphaSpend<- alpha_spendT2(N,alpha)}
                          if( AlphaSpend2==3 ){AlphaSpend<- alpha_spendT3(N,alpha)}
                          if( AlphaSpend2==4 ){
                                              if(length(gamma)>1){stop(" 'gamma' must be a single number greater than -10 and smaller than 10.",call. =FALSE)}               
                                              if(gamma=="n"){stop(" Please, choose 'gamma' between -10 and 10.",call. =FALSE)} 
                                              if(sum(is.numeric(gamma))!=1){stop("Symbols and texts are not applicable for 'gamma'. It must be a number greater than -10 and smaller than 10.",call. =FALSE)}
                                              if(gamma<=0|gamma>5){stop(" 'gamma' must be a single number greater than -10 and smaller than 10.",call. =FALSE)}
                                              AlphaSpend<- alpha_spendT4(N,alpha,gamma)
                                             }
                         }else{
                               if(is.numeric(AlphaSpend2)==FALSE){stop(" 'AlphaSpend' must be a vector containing a non-decreasing and positive sequence of numbers with maximum at alpha or a single integer among 1 to 4.",call. =FALSE)}
                               if(min(AlphaSpend2)<0|max(AlphaSpend)!=alpha){stop(" 'AlphaSpend' must be a vector containing a non-decreasing and positive sequence of numbers with maximum at alpha or a single integer among 1 to 4.",call. =FALSE)}
                              }






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



if((is.numeric(alpha)==FALSE)){stop("'alpha' must be a number in the '(0,0.5]' interval.",call. =FALSE)}
if((alpha<=0|alpha>0.5)){stop("'alpha' must be a number in the '(0,0.5]' interval.",call. =FALSE)}

if(min(z)==max(z)){
if(length(Statistic)>1){stop("'Statistic' must be selected among MaxSPRT, Pocock, OBrien-Fleming, Wang-Tsiatis.",call. =FALSE)}
if(sum(Statistic==c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis"))==0){stop("'Statistic' must be selected among MaxSPRT, Pocock, OBrien-Fleming, Wang-Tsiatis.",call. =FALSE)}
                  }else{Statistic<- "n"}

if(Statistic=="Wang-Tsiatis"&Delta=="n"){"Use a positive numeric value for 'Delta'"}
if(Statistic=="Wang-Tsiatis"&Delta<=0){"Use a positive numeric value for 'Delta'"}
if(Statistic=="Wang-Tsiatis"&Delta>0.5){"Use a positive numeric value smaller than or equal to 0.5 for 'Delta'"}

if(Tailed=="upper"|Tailed=="two"){pst<- 1/(1+max(z))}
if(Tailed=="lower"){pst<- 1/(1+max(1/z))}
if(1-pbinom(N-1,N,pst)>alpha){
Nr<- N
while(1-pbinom(Nr-1,Nr,pst)>alpha){Nr<- Nr+1}
stop(c("For this 'N' there is no solution with prob of Type I error smaller than"," ",alpha,". Use 'N' of at least"," ",Nr,"."),call. =FALSE)
                             }
                  



####
#### FUNCTIONS FOR TEST STATISTIC SCALE


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

perf<- function(RR)
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

power<- 0 ; AUX<- rep(0,N)
for(i in 1:N){power<- power+p[i,absorb[i]+1]; AUX[i]<- p[i,absorb[i]+1]}
ETS<- sum(seq(1,N)*AUX)/power; ELS<- sum(seq(1,N)*AUX)+N*(1-power)

measures<- c(power,ETS,ELS)
return(measures)
}


#######################
####################### AUXILIAR FUNCTION FOR PERFORMANCE CALCULATION UNDER TWO-TAILED TESTING

perf2<- function(RR)
{

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

power<- 0 ; AUX<- rep(0,N)
for(i in 1:N){power<- power+p[i,absorb2[i]+1]; AUX[i]<- p[i,absorb2[i]+1]}
ETS<- sum(seq(1,N)*AUX)/power; ELS<- sum(seq(1,N)*AUX)+N*(1-power)

measures<- c(power,ETS,ELS)
return(measures)
}


##########################################################################
############################################################################

#### FINDING THE THRESHOLD

if(min(z)==max(z)){z<- rep(z,N)}else{z<- rep(z,GroupSizes)}

cvs<- rep(0,N)

if(Tailed=="upper"|Tailed=="two"){ps<- 1/(1+z)}; if(Tailed=="lower"){ps<- 1/(1+1/z)}


#########################
if(Tailed!="two"){  ## HERE IS THE POINT WHERE THE CODE FOR TWO-TAILED AND ONE-TAILED TESTES DIFFER

# Auxiliar functions to run the binomial Markov Chain in a fast way:
func_aux2<- function(j,i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*dbinom(j-k,1,ps[i])))}
func_aux3<- function(i){ k<- seq(1,uc[i-1]+1); return(sum(p[i-1,k]*(1-pbinom(absorb[i]-k,1,ps[i]))))}
func_aux1<- function(i){ j<- matrix(seq(1,absorb[i]),ncol=1) ; return(apply(j,1,func_aux2,i))}

p<- matrix(0,N,N+2)    	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's

absorb = rep(0,N)		# Contains the number of events needed at time mu[i] in order to reject H0
uc<- rep(0,N)
actual<- rep(0,N)	

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------


if(AlphaSpend[1]==0){
absorb[1]<- 2
for(s in 1:absorb[1]){ p[1,s]=dbinom(s-1,1,ps[1])}		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-pbinom(absorb[1]-1,1,ps[1])			# probability of rejecting H0 at time mu[1]
                                   }else{
absorb[1]<- 0
for(s in 1:absorb[1]){ p[1,s]=dbinom(s-1,1,ps[1])}		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-pbinom(absorb[1]-1,1,ps[1])			# probability of rejecting H0 at time mu[1]

while(p[1,absorb[1]+1]> AlphaSpend[1]){
p[1,]<- 0
absorb[1]<- absorb[1]+1 ; uc[1]<- absorb[1]-1

for(s in 1:absorb[1]){ p[1,s]=dbinom(s-1,1,ps[1])}		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-pbinom(absorb[1]-1,1,ps[1])			# probability of rejecting H0 at time mu[1]
                                       }

                                         }
uc[1]<- absorb[1]-1


if(N>1){
i<- 1
cumalpha<- p[1,absorb[1]+1]
actual[1]<- cumalpha
jj<- 1

while(i<N){
i<- i+1
alphai<- 1
teste<- 0
while(alphai> AlphaSpend[jj]&teste==0){
if(sum(i==an)==0){teste<- 1; absorb[i]<- i}
alphai<- cumalpha
p[i,]<- 0
absorb[i]<- absorb[i]+1 ; uc[i]<- absorb[i]-1
       p[i,1:absorb[i]]<- func_aux1(i) # Calculates the standard p[][] cell values
       p[i,absorb[i]+1]<- func_aux3(i) # Calculates the diagonal absorbing states where H0 is rejected	

                             alphai<- alphai + p[i,absorb[i]+1] 
                                      }
cumalpha<- alphai
actual[i]<- cumalpha
if(sum(i==an)>0){jj<- jj+1}                
          } # end for i
       }

#### SAVING PERFORMANCE MEASURES 
R<- matrix(R[order(R)],length(R),1)
Measures<- apply(R,1,perf)
Measures<- cbind(R,matrix(as.numeric(Measures),ncol=3,byrow=T)) 
colnames(Measures)<- c("R","Power","ESignalTime","ESampleSize")

absorb2<- rep(0,N)
marks<- seq(1,N)[absorb<=seq(1,N)]
for(i in 1:N){
if(sum(i==marks)==0){cvs[i]<- NA;absorb2[i]<- NA}else{
absorba<- absorb
for(kk in (i-1):1){absorba[kk]<- max(absorba[kk+1]-1,0)}

if(Tailed=="lower"){
if(max(z)==min(z)){cvs[i]<- LLR(z[1:i],1-diff(c(0,absorba[1:i])),rep(1,i),Tailed)}
absorb2[i]<- i-absorba[i]
                   }
if(Tailed=="upper"){
if(max(z)==min(z)){cvs[i]<- LLR(z[1:i],diff(c(0,absorba[1:i])),rep(1,i),Tailed)}
absorb2[i]<- absorba[i]
                   }
                  
                                                     }
             }

if(max(z)!=min(z)){
res<- list(absorb2[an],actual[an],Measures) ; names(res)<- c("cvs.cases","ActualSpend","Performance")
                  }else{
res<- list(cvs[an],absorb2[an],actual[an],Measures) ; names(res)<- c("cvs","cvs.cases","ActualSpend","Performance")
                        }

                       }else{ #### HERE IS THE PART FOR TWO-TAILED TESTS



# Auxiliar functions to run the binomial Markov Chain in a fast way:
func_aux2<- function(j,i){ k<- seq(uc1[i-1]+1,uc2[i-1]+1); return(sum(p[i-1,k]*dbinom(j-k,1,ps[i])))}
func_aux3<- function(i){ k<- seq(uc1[i-1]+1,uc2[i-1]+1); return(sum(p[i-1,k]*(1-pbinom(absorb2[i]-k,1,ps[i]))))}
func_aux4<- function(i){ k<- seq(uc1[i-1]+1,uc2[i-1]+1); return(sum(p[i-1,k]*(pbinom(absorb1[i]-k+1,1,ps[i]))))}
func_aux1<- function(i){ j<- matrix(seq(absorb1[i]+2,absorb2[i]),ncol=1) ; return(apply(j,1,func_aux2,i))}

p<- matrix(0,N,N+2)    	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's



absorb2 = rep(0,N)		# Contains the number of events needed at time mu[i] in order to reject H0
absorb1 = seq(1,N)
uc1<- rep(0,N)
uc2<- rep(0,N)

actual<- rep(0,N)	

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------


if(AlphaSpend[1]==0){
absorb2[1]<- 2; absorb1[1]<- -1
for(s in absorb1[1]:absorb2[1]){ p[1,s]=dbinom(s-1,1,ps[1])}		# Probability of having s-1 cases at time mu[1]
pabs2<- 1-pbinom(absorb2[1]-1,1,ps[1])			# probability of rejecting H0 at time mu[1]
pabs1<- pbinom(absorb1[1],1,ps[1])
p[1,absorb2[1]+1]=pabs2+pabs1

                                                 }else{# OPEN ELSE
absorb2[1]<- 0
absorb1[1]<- 1

for(s in absorb1[1]:absorb2[1]){ p[1,s]=dbinom(s-1,1,ps[1])}		# Probability of having s-1 cases at time mu[1]
pabs2<- 1-pbinom(absorb2[1]-1,1,ps[1])			# probability of rejecting H0 at time mu[1]
pabs1<- pbinom(absorb1[1],1,ps[1])

while(pabs2> AlphaSpend[1]/2){
pabs2<- 0
absorb2[1]<- absorb2[1]+1
pabs2 =1-pbinom(absorb2[1]-1,1,ps[1])			# probability of rejecting H0 by the upper side
                             }
while(pabs1> AlphaSpend[1]/2){
pabs1<- 0
absorb1[1]<- absorb1[1]-1
pabs1 =pbinom(absorb1[1],1,ps[1])			# probability of rejecting H0 by the lower side
                             }

for(s in (absorb1[1]+2):absorb2[1]){ p[1,s]=dbinom(s-1,1,ps[1])}		# Probability of having s-1 cases at time mu[1]
p[1,absorb2[1]+1]=pabs2+pabs1
                                                       }# CLOSE ELSE
uc2[1]<- absorb2[1]-1
uc1[1]<- absorb1[1]+1

#########


if(N>1){
i<- 1
cumalpha2<- pabs2
cumalpha1<- pabs1 
actual[1]<- cumalpha2+cumalpha1


jj<- 1
while(i<N){

i<- i+1

alphai<- 1
teste<- 0
while(alphai> AlphaSpend[jj]/2&teste==0){
if(sum(i==an)==0){teste<- 1; absorb2[i]<- i}
alphai<- cumalpha2
       absorb2[i]<- absorb2[i]+1 ; uc2[i]<- absorb2[i]-1       	
       pabs2<- func_aux3(i) # Calculates the diagonal absorbing states where H0 is rejected
                             alphai<- alphai + pabs2 
                                      }
cumalpha2<- alphai

alphai<- 1
teste<- 0
while(alphai> AlphaSpend[jj]/2&teste==0){
if(sum(i==an)==0){teste<- 1; absorb1[i]<- 0}
alphai<- cumalpha1
       absorb1[i]<- absorb1[i]-1 ; uc1[i]<- absorb1[i]+1       	
       pabs1<- func_aux4(i) # Calculates the diagonal absorbing states where H0 is rejected
                             alphai<- alphai + pabs1 
                                      }
cumalpha1<- alphai


p[i,(absorb1[i]+2):absorb2[i]]<- func_aux1(i) # Calculates the standard p[][] cell values
p[i,absorb2[i]+1]<- pabs2+pabs1


actual[i]<- cumalpha2+cumalpha1
if(sum(i==an)>0){jj<- jj+1}                
          } # end for i
       }

#### SAVING PERFORMANCE MEASURES
R<- matrix(R[order(R)],length(R),1)
Measures<- apply(R,1,perf2)
Measures<- cbind(R,matrix(as.numeric(Measures),ncol=3,byrow=T)) 
colnames(Measures)<- c("RR","Power","ESignalTime","ESampleSize")

absorbl<- rep(0,N)
absorbu<- rep(0,N)

if(max(z)==min(z)){
marks<- seq(1,N)[absorb1>=0]
cvsl<- rep(0,N)
for(i in 1:N){
if(sum(i==marks)==0){cvsl[i]<- NA;absorb1[i]<- NA}else{
absorbal<- absorb1
for(kk in (i-1):1){if(kk>=absorb1[i]){absorbal[kk]<- absorb1[i]}else{absorbal[kk]<- kk}}
cvsl[i]<- LLR(z[1:i],diff(c(0,absorbal[1:i])),rep(1,i),Tailed)
                                      }
             }
                  }



if(max(z)==min(z)){
marks<- seq(1,N)[absorb2<=seq(1,N)]
cvsu<- rep(0,N)
for(i in 1:N){
if(sum(i==marks)==0){cvsu[i]<- NA;absorb2[i]<- NA}else{
absorbau<- absorb2
for(kk in (i-1):1){absorbau[kk]<- max(absorbau[kk+1]-1,0)}
cvsu[i]<- LLR(z[1:i],diff(c(0,absorbau[1:i])),rep(1,i),Tailed)
                                      }
             }
                 }

if(max(z)!=min(z)){
res<- list(absorb1[an],absorb2[an],actual[an],Measures) ; names(res)<- c("cases.lower","cases.upper","ActualSpend","Performance")
                  }else{
res<- list(cvsl[an],cvsu[an],absorb1[an],absorb2[an],actual[an],Measures) ; names(res)<- c("cvs.lower","cvs.upper","cases.lower","cases.upper","ActualSpend","Performance")
                         }
                            }#### HERE CLOSES THE POINT WHERE THE CODE FOR TWO-TAILED AND ONE-TAILED TESTS DIFFER 



return(res)

} #end function Performance.AlphaSpend.Binomial


#### EXAMPLES

#res<- Performance.AlphaSpend.Binomial(N=50,alpha=0.05,AlphaSpend=1,z=c(1,1.5,2,1.3),p="n",GroupSizes=c(15,15,10,10),Tailed="upper",rho=0.5,gamma="n",RR=2,Statistic= "MaxSPRT",Delta="n")

#res<- Performance.AlphaSpend.Binomial(N=50,alpha=0.05,AlphaSpend=1,z=c(1,1.5,2,1.3),p="n",GroupSizes=c(15,15,10,10),Tailed="upper",RR=2,rho=0.5)

#res<- Performance.AlphaSpend.Binomial(N=50,alpha=0.05,AlphaSpend=c(0.03,0.04,0.05),z=c(1,2,1.5),p="n",GroupSizes=c(20,20,10),Tailed="two",RR=2,Statistic="MaxSPRT")

#res<- Performance.AlphaSpend.Binomial(N=50,alpha=0.05,AlphaSpend=c(0.03,0.04,0.05),z=c(1,2,1.5),p="n",GroupSizes=c(20,20,10),Tailed="two",RR=2,Statistic="Pocock")

#res<- Performance.AlphaSpend.Binomial(N=50,alpha=0.05,AlphaSpend=c(0.03,0.04,0.05),z=c(1,1,1),p="n",GroupSizes=c(20,20,10),Tailed="two",RR=2,Statistic="Pocock")

#res<- Performance.AlphaSpend.Binomial(N=55,alpha=0.05,AlphaSpend=1,z=c(1,1.5,2),p="n",GroupSizes=c(20,20,15),Tailed="lower",RR=c(1.5,1.2),Statistic="MaxSPRT",rho=2)


