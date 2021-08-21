Performance.AlphaSpend.Poisson<- function(SampleSize, alpha=0.05,D=0,M=1,RR,alphaSpend=1,rho=0.5,gamma="n",Statistic=c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis"),Delta="n",Tailed="upper"){


# ------------------- INPUT VARIABLE ----------------------------------------------------------
# SampleSize = maximum length of surveillance, defined in terms of expected counts under H0. It is used only when cvs.lower and cvs.upper are non-empty. Otherwise, SampleSize is ignored since cases.lower and/or cases.upper define the maximum sample size in the scale of the events. 
# alpha = overall size of the test
# alphaSpend: is a number among 1, 2, 3, 4. Each representing one of the shapes shown in the paper.
# RR= relative risk for performance calculation
# rho: is a parameter when alphaSpend=1 
# gamma: is a parameter when alphaSpend=4 
# Statistic = test statistic scale that the user wants for the threshold output 
# Delta = number in (0, 0.5] for the test statistic of Wang-Tsiatis
# If Tailed="upper" (default), then the threshold is given as an upper tailed testing (H0: R<=RR0), Tailed="lower" for lower tailed (H0: R>=RR0), and Tailed="two" for two-tailed testing (H0: R=RR0).

pho<- rho

# CHECKING THE VALIDITY OF INPUT PARAMETERS

if(length(Tailed)!=1){stop("For this version of the Sequential package, Threshold.Poisson works only for 'Tailed=upper'.",call. =FALSE)}
if(Tailed!="upper"){stop("For this version of the Sequential package, Threshold.Poisson works only for 'Tailed=upper'.",call. =FALSE)}
if(sum(Tailed==c("upper","lower","two"))==0){stop(" 'Tailed' must be chosen among 'upper', 'lower' or 'two'.",call. =FALSE)}

if(length(SampleSize)>1){stop("The maximum length of surveillance, 'SampleSize', must be a single positive number.",call. =FALSE)}
if(is.numeric(SampleSize)==FALSE){stop("The maximum length of surveillance, 'SampleSize', must be a positive number.",call. =FALSE)}
if(SampleSize<=0){stop("The maximum length of surveillance, 'SampleSize', must be a positive number.",call. =FALSE)}

if(length(alphaSpend)>1){stop("The 'alphaSpend' parameter must contain a single number among 1, 2, 3 and 4. Read the user manual for more details on the usage of the Threshold.Poisson function.",call. =FALSE)}
if(length(alphaSpend)>1){stop("The 'alphaSpend' parameter must contain a single number among 1, 2, 3 and 4. Read the user manual for more details on the usage of the Threshold.Poisson function.",call. =FALSE)}


if(length(alpha)>1){stop("The 'alpha' level must be a single positive number in the (0,1) interval.",call. =FALSE)}
if(is.numeric(alpha)==FALSE){stop("The 'alpha' level must be a single positive number in the (0,1) interval.",call. =FALSE)}
if(alpha<=0|alpha>=1){stop("The 'alpha' level must be a single positive number in the (0,1) interval.",call. =FALSE)}

if(length(rho)>1){stop("The parameter 'rho' for 'alphaSpend=1' must be a single positive number.",call. =FALSE)}
if(is.numeric(rho)==FALSE){stop("The parameter 'rho' for 'alphaSpend=1' must be a single positive number.",call. =FALSE)}
if(rho<=0){stop("The parameter 'rho' for 'alphaSpend=1' must be a single positive number.",call. =FALSE)}

if(alphaSpend==4){
if(length(gamma)>1){stop("The parameter 'gamma' for 'alphaSpend=4' must be a single positive number.",call. =FALSE)}
if(is.numeric(gamma)==FALSE){stop("The parameter 'gamma' for 'alphaSpend=4' must be a single positive number.",call. =FALSE)}
if(gamma<=0){stop("The parameter 'gamma' for 'alphaSpend=4' must be a single positive number.",call. =FALSE)}
                 }


if(Statistic=="Wang-Tsiatis"){
if(length(Delta)>1){stop("The parameter 'Delta' for 'Statistic=Wang-Tsiatis' must be a single number in the (0, 0.5] interval.",call. =FALSE)}
if(is.numeric(Delta)==FALSE){stop("The parameter 'Delta' for 'Statistic=Wang-Tsiatis' must be a single number in the (0, 0.5] interval.",call. =FALSE)}
if(Delta<=0|Delta>0.5){stop("The parameter 'Delta' for 'Statistic=Wang-Tsiatis' must be a single number in the (0, 0.5] interval.",call. =FALSE)}
                             }



T<- SampleSize

#### Here are the functions for each shape

#1 Type t^pho, Kim and DeMets (1987a) and Jennison & Turnball(1989,1990) noted that this function produces Pocock and O'Brien & Fleming. Values 1,2 e 3 for similarity with Wang & Tsiatis(1987) test

alpha_spendT1<- function(tt)
{
x<- tt/T
return(alpha*(x^pho))
}


#2 Type Normal distribution. Lan & DeMets(1983) suggested this for similarity with O'Brien & Fleming  

alpha_spendT2<- function(tt)
{
x<- tt/T
za<- qnorm(1-alpha/2)
return(2-2*pnorm(za/sqrt(x)))
}


#3 Type LOG-EXP, Lan & DeMets(1983) suggested this to mimic Pocock's test

alpha_spendT3<- function(tt)
{
x<- tt/T
return(alpha*log(1+(exp(1)-1)*x))
}


#4 Type gamma, Hwang, Shih & DeCani(1990) 

alpha_spendT4<- function(tt)
{
x<- tt/T
if(gamma==0){return(alpha*x)}else{
return(alpha*(1-exp(-gamma*x))/(1-exp(-gamma)))
                                  }
}

## Defining the type of alpha spending

if(length(alphaSpend)>1){alphashape<- alphaSpend}else{
if(alphaSpend==1){alphashape<- alpha_spendT1}
if(alphaSpend==2){alphashape<- alpha_spendT2}
if(alphaSpend==3){alphashape<- alpha_spendT3}
if(alphaSpend==4){alphashape<- alpha_spendT4}
                                                     }
#----- MAXSPRT STATISTIC

LLR <- function(cc,uu,Tai=Tailed) {

if(Tai=="upper"){
	if(cc<=uu) x=0
	if(cc>uu) x = (uu-cc) + cc*log(cc/uu)
                   }
if(Tai=="lower"){
	if(cc>=uu) x=0
	if(cc<uu) x = (uu-cc) + cc*log(cc/uu)
                   }
if(Tai=="two"){	
	x = (uu-cc) + cc*log(cc/uu)
                   }
	return(x)
	}
             
#--------------------------

#----- THE Pocock (1977) statistic
LLR2 <- function(cc,uu,Tai=Tailed)
{

if(Tai=="upper"){
	if(cc<=uu) x=0
     if(cc>uu) x<- (cc-uu)/sqrt(uu)
                   }

if(Tai=="lower"){
	if(cc>=uu) x=0
      if(cc<uu) x<- (cc-uu)/sqrt(uu)
                   }   

if(Tai=="two"){
    x<- (cc-uu)/sqrt(uu)
                 } 
return(abs(x))

}

#----- THE OBrien and Fleming (1972) statistic
LLR3 <- function(cc,uu,Tai=Tailed)
{

if(Tai=="upper"){
	if(cc<=uu) x=0
     if(cc>uu) x<- (cc-uu)/sqrt(uu)
                   }

if(Tai=="lower"){
	if(cc>=uu) x=0
      if(cc<uu) x<- (cc-uu)/sqrt(uu)
                   }   

if(Tai=="two"){
    x<- (cc-uu)/sqrt(uu)
                 } 
return(sqrt(uu/SampleSize)*abs(x))

}


#---- Wang e Tsiatis(1987), which is also a Pocock type, but with threshold multiplied by ((i/N)^(Delta-0.5)). Values of Delta can be 0.1, 0.25, and 0.4, page 40 by Jennison and Turniball(2000)

LLR4 <- function(cc,uu,Tai=Tailed)
{

if(Tai=="upper"){
	if(cc<=uu) x=0
     if(cc>uu) x<- (cc-uu)/sqrt(uu)
                   }

if(Tai=="lower"){
	if(cc>=uu) x=0
      if(cc<uu) x<- (cc-uu)/sqrt(uu)
                   }   

if(Tai=="two"){
    x<- (cc-uu)/sqrt(uu)
                 } 
return(((uu/SampleSize)^(0.5-Delta))*abs(x))

}



##### Defining the test statistic to work in the rest of the code

if(Statistic=="Pocock"){LLR<- LLR2}; if(Statistic=="OBrien-Fleming"){LLR<- LLR3}; if(Statistic=="Wang-Tsiatis"){LLR<- LLR4}






mu<- rep(0,round(T*4))
auxD<- 0
while(auxD==0){
mu1<- 0 ; mu2<- T
mut<- (mu1+mu2)/2
alphas<- alphashape(mut)
while(abs(1-ppois(M-1,mut)-alphas)>0.0000001&abs(mut-T)>0.000001){
if(1-ppois(M-1,mut)>alphas){mu2<- mut}else{mu1<- mut}
mut<- (mu1+mu2)/2
alphas<- alphashape(mut)
                                                                 }
mu[M]<- mut; if(mu[M]<D){M<- M+1}else{auxD<- 1}
              }


# Defining the p[][] matrix
# -------------------------

p = matrix(0,round(T*4),round(T*4)+1)				# p[i,j] is the probability of having j-1 cases at time mu[i]
p1<- matrix(0,round(T*4),round(T*4)+1)

# Calculating the M-th row in the p[][] matrix for which there is a chance to reject H0

for(s in 1:M){ 
p[M,s]=dpois(s-1,mu[M])		  # Probability of having s-1 cases at time mu[M] under H0
p1[M,s]=dpois(s-1,RR*mu[M])	  # Probability of having s-1 cases at time mu[M] under the alternative	
             }
p[M,M+1]=1-ppois(M-1,mu[M])     # Probability of rejecting H0 at time mu[M] under H0
p1[M,M+1]=1-ppois(M-1,RR*mu[M]) # Probability of rejecting H0 at time mu[M] under the alternative

# Defining the starting target alpha spend
target<- rep(0,ncol(p))
target[M]<- alphas

# Looping to find the solutions for each i>=M+1
i<- M+1
TypeI<- p[M,M+1]
while(abs(mu[i-1]-T)>0.000001){


mu1<- mu[i-1] ; mu2<- T
mut<- (mu1+mu2)/2
alphas<- alphashape(mut)
teste<- 0
while(teste==0&abs(mu1-mu2)>0.00000000001){ # open while 2

for(j in 1:(i-1)){							# This loop calculates the p[][] matrix, one column at a time, from left to right
		for(k in 1:j){ 
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,mut-mu[i-1])	# Calculates the standard p[][] cell values
                  p1[i,j]=p1[i,j]+p1[i-1,k]*dpois(j-k,(mut-mu[i-1])*RR)
                         } 
                 }
	for(k in 1:(i-1)){
		p[i,i]=p[i,i]+p[i-1,k]*dpois(i-k,mut-mu[i-1])		# Calculates the diagonal under the absorbing states, which requires a unique formula
            p1[i,i]=p1[i,i]+p1[i-1,k]*dpois(i-k,(mut-mu[i-1])*RR)
                 }
	for(k in 1:(i-1)){ 
		p[i,i+1]=p[i,i+1]+p[i-1,k]*(1-ppois(i-k,mut-mu[i-1]))# Calculates the diagonal absorbing states where H0 is rejected
            p1[i,i+1]=p1[i,i+1]+p1[i-1,k]*(1-ppois(i-k,(mut-mu[i-1])*RR))
                       }

if(abs(TypeI+p[i,i+1]-alphas)>0.0000001&abs(mu1-mu2)>0.00000000001){
if(TypeI+p[i,i+1]>alphas){mu2<- mut}else{mu1<- mut}
mut<- (mu1+mu2)/2
alphas<- alphashape(mut)
p[i,]<- 0 ; p1[i,]<- 0
                                                                   }else{teste<- 1}

                                          } # close while 2

target[i]<- alphas
TypeI<- TypeI+p[i,i+1]
mu[i]<- mut
i<- i+1

               }# CLOSE while(abs(mu[i-1]-T)>0.000001))

# Calculating expected sample size, expected number of cases to signal, power, and threshold
power<- p1[1,1+1]
ETS<- mu[1]*p1[1,1+1] 

CVS<- rep(0,i)
alphaspend<- rep(0,ncol(p)); alphaspend[1]<- p[1,2]
for(j in 2:i){power<- power+p1[j,j+1]; ETS<- ETS+mu[j]*p1[j,j+1];alphaspend[j]<- p[j,j+1]+alphaspend[j-1];CVS[j]<- LLR(j,mu[j])}
ETS<- ETS/power

# Calculating expected sample size
ESS<- ETS+(1-power)*SampleSize

perf<- c(RR,power,ETS,ESS)
names(perf)<- c("RR","power","ESignalTime","ESampleSize")
results<- list(CVS[-length(CVS)],perf)
names(results)<- c("cvs","Performance")

return(results)

} # CLOSES FUNCTION
##################################


# EXAMPLE

#res<- Performance.AlphaSpend.Poisson(SampleSize=30, alpha=0.05,alphaSpend=1,RR=1.5,rho=0.5,gamma="n",Statistic="MaxSPRT")









