


## This version was implemented on April 2018

Performance.Threshold.Poisson <- function(SampleSize,CV.lower="n",CV.upper="n",GroupSizes="n",Tailed="upper",Statistic=c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis"),Delta="n",RR)
{

PRECISION=0.00000001
cvs.lower<- CV.lower
cvs.upper<- CV.upper

# ------------------- INPUT VARIABLE ----------------------------------------------------------
# SampleSize = maximum length of surveillance, defined in terms of expected counts under H0. It is used only when CV.lower and CV.upper are non-empty. Otherwise, SampleSize is ignored since CV.lower and/or CV.upper define the maximum sample size in the scale of the events. 
# cv = critical value
# M = The minimum number of cases for which a signal is allowed to occur
# D = Time < T for first look at the data, defined in terms of the expected counts under H0
# If Tailed="upper" (default), then the threshold is given as an upper tailed testing (H0: R<=RR0), Tailed="lower" for lower tailed (H0: R>=RR0), and Tailed="two" for two-tailed testing (H0: R=RR0).
# GroupSizes: Vector with the number of events (exposed+unexposed) between two looks at the data, i.e, irregular group sizes. Important: Must sums up N. For continuos sequential analysis, specify GroupSizes=1
# RR= vector of relative risks for performance calculation
# Statistic = scale that the user is using for the cvs inputs
# Delta = number in (0, 0.5] for the test statistic of Wang-Tsiatis

if(sum(Tailed==c("upper","lower","two"))==0){stop(" 'Tailed' must be chosen among 'upper', 'lower' or 'two'.",call. =FALSE)}

if(length(GroupSizes)==1){if(GroupSizes=="n"){method<- "continuous"}else{method<- "group"}}
if(length(GroupSizes)>1){method<- "group"}

M<- 1
D<- 0
Groups<- GroupSizes
T<- SampleSize
L<- T
MinCases<- M
Late<- D

if(length(SampleSize)>1){stop("The maximum length of surveillance, 'SampleSize', must be a single positive number.",call. =FALSE)}

if(length(Groups)==1){# 1

   if(Groups!="n"){
if(is.numeric(Groups)==FALSE){stop("'GroupSizes' must be a single number or a vector summing up 'SampleSize'.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'GroupSizes' must be a single number or a vector summing up 'SampleSize'.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'Groups' must be a single number or a vector summing up 'SampleSize'.",call. =FALSE)}

if(Groups==0){stop("'Groups' must be a single number or a vector summing up 'SampleSize'.",call. =FALSE)}
if(SampleSize/Groups!=round(SampleSize/Groups)){stop("If 'GroupSizes' is a single number, then the maximum length of surveillance, 'SampleSize', must be a multiple of 'GroupSizes'.",call. =FALSE)}
if(Groups>SampleSize){stop("If 'GroupSizes' is a single number, then the maximum length of surveillance, 'SampleSize', must be a multiple of 'GroupSizes'.",call.=FALSE)}
if(SampleSize/Groups==round(SampleSize/Groups)){Groups<- rep(Groups,SampleSize/Groups);GroupSizes<- Groups}
if(length(cvs.lower)==1){if(is.numeric(cvs.lower)==TRUE){cvs.lower<- rep(cvs.lower,length(Groups))}}
if(length(cvs.upper)==1){if(is.numeric(cvs.upper)==TRUE){cvs.upper<- rep(cvs.upper,length(Groups))}}
                   }
                    }# close 1



if(length(Groups)>1){
if(sum(is.numeric(Groups))==0){stop("'GroupSizes' must be a single number or a vector summing up 'SampleSize'.",call. =FALSE)}else{
if(is.numeric(SampleSize)==FALSE){stop("The maximum length of surveillance, 'SampleSize', must be a positive number.",call. =FALSE)}
if(SampleSize<=0){stop("The maximum length of surveillance, 'SampleSize', must be a positive number.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'GroupSizes' must be a single number or a vector summing up 'SampleSize'.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'GroupSizes' must be a single number or a vector summing up 'SampleSize'.",call. =FALSE)}
if(sum(Groups)!=SampleSize){stop("'GroupSizes' must be a single number or a vector summing up 'SampleSize'.",call. =FALSE)}
if(length(cvs.lower)==1){if(is.numeric(cvs.lower)==TRUE){cvs.lower<- rep(cvs.lower,length(Groups))}}
if(length(cvs.upper)==1){if(is.numeric(cvs.upper)==TRUE){cvs.upper<- rep(cvs.upper,length(Groups))}}
}
}



if(Tailed=="upper"){
if(length(cvs.lower)==1){if(cvs.lower!="n"){stop("For 'Tailed=upper' it is not possible to use 'cvs.lower'.",call. =FALSE)}}
if(length(cvs.lower)>1){stop("For 'Tailed=upper' it is not possible to use 'cvs.lower'.",call. =FALSE)}
                   }

if(Tailed=="lower"){
if(length(cvs.upper)==1){if(cvs.upper!="n"){stop("For 'Tailed=lower' it is not possible to use 'cvs.upper'.",call. =FALSE)}}
if(length(cvs.upper)>1){stop("For 'Tailed=lower' it is not possible to use 'cvs.upper'.",call. =FALSE)}
                   }

### Checking the threshold scale choice 
aux_thres<- rep(0,2)
for(jj in 1:length(cvs.lower)){if(cvs.lower[jj]!="n"){aux_thres[1]<- 1}}
for(jj in 1:length(cvs.upper)){if(cvs.upper[jj]!="n"){aux_thres[2]<- 1}}

if(Tailed=="two"){
  if(length(GroupSizes)==1){
if(length(cvs.lower)>1|length(cvs.upper)>1&method=="group"){stop("Signaling thresholds must have the same dimension of 'GroupSizes'.",call. =FALSE)}
if(cvs.lower=="n"|cvs.upper=="n"){stop("Lower and upper signaling thresholds must be specified.",call. =FALSE)}
                           }

if(length(GroupSizes)>1){
if(aux_thres[1]+aux_thres[2]!=2&aux_thres[3]+aux_thres[4]!=2){stop("Lower and upper signaling thresholds must be specified.",call. =FALSE)}
if(aux_thres[1]+aux_thres[2]==2){if(length(cvs.lower)!=length(GroupSizes)|length(cvs.upper)!=length(GroupSizes)|length(cvs.upper)!=length(cvs.lower)&method=="group"){stop("Lower and upper signaling thresholds must have the same dimension of 'GroupSizes'.",call. =FALSE)}}
                        }
                 }

if(length(GroupSizes)==1){
if(GroupSizes=="n"){if(Tailed=="two"){if(length(cvs.lower)==1){cvs.lower<- rep(cvs.lower,2*T)};if(length(cvs.upper)==1){cvs.upper<- rep(cvs.upper,2*T)}}}
if(GroupSizes=="n"){if(Tailed=="lower"){if(length(cvs.lower)==1){cvs.lower<- rep(cvs.lower,2*T)}}}
if(GroupSizes=="n"){if(Tailed=="upper"){if(length(cvs.upper)==1){cvs.upper<- rep(cvs.upper,2*T)}}}

if(GroupSizes=="n"&Tailed=="two"){if(length(cvs.lower)<2*T|length(cvs.upper)<2*T){stop("For variable lower and upper signaling thresholds, the dimension of cvs.lower and cvs.upper must be at least twice the value of SampleSize.",call. =FALSE)}}
if(GroupSizes=="n"&Tailed=="lower"){if(length(cvs.lower)<2*T){stop("For variable lower signaling threshold, the dimension of cvs.lower must be at least twice the value of SampleSize.",call. =FALSE)}}
if(GroupSizes=="n"&Tailed=="upper"){if(length(cvs.upper)<2*T){stop("For variable lower signaling threshold, the dimension of cvs.upper must be at least twice the value of SampleSize.",call. =FALSE)}}
                         }


if(aux_thres[1]==1|aux_thres[2]==1){stat_sca<- 1}else{stat_sca<- 0}



####
#### FUNCTIONS FOR THE TEST STATISTIC SCALE


########################################### 
#----- THE MAXSPRT STATISTIC (Wald-type) 

LLR <- function(cc,uu,Tailed) {

if(Tailed=="upper"){
	if(cc<=uu) x=0
	if(cc>uu) x = (uu-cc) + cc*log(cc/uu)
                   }
if(Tailed=="lower"){
	if(cc>=uu) x=0
	if(cc<uu){ if(cc>0){x = (uu-cc) + cc*log(cc/uu)}else{x= uu} }
                   }
if(Tailed=="two"){	
	if(cc>0){x = (uu-cc) + cc*log(cc/uu)}else{x= uu}
                 }
	return(x)
	                        }
             
#--------------------------


#----- THE Pocock (1977) statistic
LLR2 <- function(cc,uu,Tailed)
{

if(Tailed=="upper"){
	if(cc<=uu) x=0
     if(cc>uu) x<- (cc-uu)/sqrt(uu)
                   }

if(Tailed=="lower"){
	if(cc>=uu) x=0
      if(cc<uu) x<- (cc-uu)/sqrt(uu)
                   }   

if(Tailed=="two"){
    x<- (cc-uu)/sqrt(uu)
                 } 
return(abs(x))

}

#----- THE OBrien and Fleming (1972) statistic
LLR3 <- function(cc,uu,Tailed)
{

if(Tailed=="upper"){
	if(cc<=uu) x=0
     if(cc>uu) x<- (cc-uu)/sqrt(uu)
                   }

if(Tailed=="lower"){
	if(cc>=uu) x=0
      if(cc<uu) x<- (cc-uu)/sqrt(uu)
                   }   

if(Tailed=="two"){
    x<- (cc-uu)/sqrt(uu)
                 } 
return(sqrt(uu/SampleSize)*abs(x))

}


#---- Wang e Tsiatis(1987), which is also a Pocock type, but with threshold multiplied by ((i/N)^(Delta-0.5)). Values of Delta can be 0.1, 0.25, and 0.4, page 40 by Jennison and Turniball(2000)

LLR4 <- function(cc,uu,Tailed)
{

if(Tailed=="upper"){
	if(cc<=uu) x=0
     if(cc>uu) x<- (cc-uu)/sqrt(uu)
                   }

if(Tailed=="lower"){
	if(cc>=uu) x=0
      if(cc<uu) x<- (cc-uu)/sqrt(uu)
                   }   

if(Tailed=="two"){
    x<- (cc-uu)/sqrt(uu)
                 } 
return(((uu/SampleSize)^(0.5-Delta))*abs(x))

}



##### Defining the test statistic to work in the rest of the code

if(Statistic=="Pocock"){LLR<- LLR2}; if(Statistic=="OBrien-Fleming"){LLR<- LLR3}; if(Statistic=="Wang-Tsiatis"){LLR<- LLR4}



#######################
####################### AUXILIAR FUNCTION FOR PERFORMANCE CALCULATION FOR "method=group"

perf<- function(RR,res=0)
{

ress<- res

#####################################################
####### HERE STARTS THE CODE FOR UPPER-TAILED TESTING
#####################################################

if(Tailed=="upper"){

CV<- cvs.upper
					
imax<- length(GroupSizes) # The number of tests performed, including final time T.
mmu<- GroupSizes
mu0<- mmu%*%upper.tri(matrix(0,length(mmu),length(mmu)),diag=T)*1
mu0<- as.vector(mu0)                  # An array of the expected counts at each of the tests
		
								# mu[i] is the commulative expected count at the i'th test
mu0[imax]=T							# Sets the expected count at the last test to equal T

absorb = seq(length=imax,from=0,by=0)		# Contains the number of events needed at time mu[i] in order to reject H0
for(i in 1:imax)
	while( LLR(absorb[i],mu0[i],Tailed="upper")<CV[i] ) 
		absorb[i]=absorb[i]+1;

mu=mu0*RR

AlphaSpend<- rep(0,imax)

absorb[absorb<MinCases]<- MinCases

p = seq(length=imax*(absorb[imax]+1), from=0, by=0)	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's
dim(p) = c(imax,absorb[imax]+1)				# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

for(s in 1:absorb[1]) p[1,s]=dpois(s-1,mu[1])		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-ppois(absorb[1]-1,mu[1])			# probability of rejecting H0 at time mu[1]
AlphaSpend[1]<- p[1,absorb[1]+1]

# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------
if(imax>1){
for(i in 2:imax) {
	for(j in 1:absorb[i])					# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:min(j,absorb[i-1]))
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,mu[i]-mu[i-1])	# Calculates the standard p[][] cell values
	for(k in 1:absorb[i-1]) p[i,absorb[i]+1]=p[i,absorb[i]+1]+p[i-1,k]*(1-ppois(absorb[i]-k,mu[i]-mu[i-1]))

AlphaSpend[i]<- p[i,absorb[i]+1]
} # end for i	
          }


# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

power=0
for(i in 1:imax){ power=power+p[i,absorb[i]+1]	}				
time=0
for(i in 1:imax){time=time+mu0[i]*p[i,absorb[i]+1]}
signaltime=time/power					
surveillancetime<- time+(1-power)*L

         }#  close if related to if(Tailed=="upper")

#####################################################
####### HERE CLOSES THE CODE FOR UPPER-TAILED TESTING
#####################################################




#####################################################
####### HERE STARTS THE CODE FOR LOWER-TAILED AND TWO-TAILED TESTING
#####################################################

if(Tailed=="two"|Tailed=="lower"){# open



imin<- M

mmu<- GroupSizes       # An array of the expected counts at each test

mu<- mmu%*%upper.tri(matrix(0,length(mmu),length(mmu)),diag=T)*1
mu<- as.vector(mu)
N<- length(mu)
								# mu[i] is the commulative expected count at the i'th test
mu[N]=T							# Sets the expected count at the last test equal to T
mtemp = c(0,mu)
mmu = diff(mtemp) 		                  #The marginal difference of the mu[] vector

RRmu<- mu*RR



####----------------------------
#### FINDING absorb
####----------------------------

absorb<- rep(-1,N); absorb2<- rep(-1,N)


for(l in 1:length(mu)){# 1

llraux<- 0
cc<- mu[l]; while( llraux<cvs.lower[l]& cc>=0 ){cc<- cc-1; llraux<- LLR(cc,mu[l],Tailed="lower")}; absorb[l]<- cc  

if(Tailed=="two"){
llraux<- 0
cc<- mu[l]; while(llraux<cvs.upper[l]){cc<- cc+1; llraux<- LLR(cc,mu[l],Tailed="two")}; absorb2[l]<- cc 
                 }

                     }# close 1

imax<- max(absorb,absorb2)



p = matrix(0,N,imax+1)	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's


# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------
AlphaSpend_lower<- rep(0,N)
AlphaSpend_upper<- rep(0,N)

if(Tailed=="upper"){
for(s in 1:absorb2[1]) p[1,s]=dpois(s-1,RRmu[1])		# Probability of having s-1 cases at time mu[1]
#p[1,absorb2[1]+1]=1-ppois(absorb2[1]-1,RRmu[1])			# probability of rejecting H0 at time mu[imin]
AlphaSpend_upper[1]<- 1-ppois(absorb2[1]-1,RRmu[1])
                   }

if(Tailed=="lower"){
for(s in absorb[1]:(absorb[2])) p[1,s]=dpois(s-1,RRmu[1]) # Probability of having s-1 cases at time mu[1] 
AlphaSpend_lower[1]<- ppois(absorb[1]-1,RRmu[1]) # probability of rejecting H0 at time mu[imin] 
                   }

if(Tailed=="two"){
 for(s in 1:max(absorb[2],absorb2[1])) p[1,s]=dpois(s-1,RRmu[1])
    AlphaSpend_upper[1]<- 1-ppois(absorb2[1]-1,RRmu[1])
    AlphaSpend_lower[1]<- ppois(absorb[1]-1,RRmu[1])
                 }


# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------

if(N>1){# 1

j1_old<- 1; if(Tailed=="two"|Tailed=="upper"){j2_old<- absorb2[1]}else{j2_old<- absorb[2]+1}

for(i in 2:N){
 
 
     j1<- max(1,absorb[i-1]+1,j1_old) 
  
      if(absorb2[i]> -1){j2<- absorb2[i]}else{ii<- i; while(absorb2[ii]== -1&ii<N){ii<- ii+1}; j2<- absorb2[ii]}; if(j2== -1){j2<- j1+1}

	for(j in j1:j2)					# This loop calculates the p[][] matrix, one column at a time, from left to right
		for(k in j1:min(j,j2_old))
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,RRmu[i]-RRmu[i-1])	# Calculates the standard p[][] cell values
           

	 if(absorb2[i]> -1){for(k in j1:j2_old){ AlphaSpend_upper[i]=AlphaSpend_upper[i]+p[i-1,k]*(1-ppois(j2-k,RRmu[i]-RRmu[i-1]))}}
      # if(absorb[i]> -1){for(k in j1_old:j1){ AlphaSpend_lower[i]=AlphaSpend_lower[i]+p[i-1,k]*ppois(j1_old-1-(k-1),RRmu[i]-RRmu[i-1])}}
      if(absorb2[i]== -1){AlphaSpend_lower[i]= p[i,i]}

j2_old<- j2 ; j1_old<- j1

} # end for i	
        }# close 1


# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

AlphaSpend<- AlphaSpend_upper+AlphaSpend_lower

power=0
for(i in 1:N){ power=power+ AlphaSpend[i]	}				
time=0
for(i in 1:N){time=time+mu[i]*( AlphaSpend[i])}
signaltime=time/power					
surveillancetime<- time+(1-power)*L




                                 }# close if(Tailed=="two"|Tailed=="lower")

#####################################################
####### HERE CLOSES THE CODE FOR LOWER-TAILED AND TWO-TAILED TESTING
#####################################################


if(ress==0){
result<- c(power,signaltime,surveillancetime)
          }else{
if(Tailed=="two"){
events_lower<- absorb; if(sum(events_lower==-1)>0){events_lower[events_lower==-1]<- NA}
events_upper<- absorb2; if(sum(events_upper==-1)>0){events_upper[events_upper==-1]<- NA}
result<- list(AlphaSpend_lower,AlphaSpend_upper,events_lower,events_upper); names(result)<- c("AlphaSpend_lower","AlphaSpend_upper","events_lower","events_upper")
                 }
if(Tailed=="lower"|Tailed=="upper"){result<- AlphaSpend}
               }
return(result)

}


##############################################################################################
#######################
####################### AUXILIAR FUNCTIONS FOR PERFORMANCE CALCULATION FOR "method=continuous"
#######################
##############################################################################################

##### Function for finding threshold in the scale of the Poisson time index given a fixed CV
####------------------------------------------------------------

tals<- function(cc,CV,Tai)
{

if(Tai=="upper"){

tau1<- 0; tau2<- 2*T
   CVaux<- 0
while(abs(CVaux-CV)>PRECISION){
    tauM<- (tau1+tau2)/2; if(tauM<0.0000001){stop(c("Use cvs.upper smaller than ",CVaux),call. =FALSE)}    
    CVaux<- LLR(cc,tauM,Tai)
    if(CVaux> CV){tau1<- tauM}else{tau2<- tauM}
                               }
                 }

if(Tai=="lower"){   
   tau1<- 0; tau2<- T
   CVaux<- 0
   scape<- ceiling(log(T/PRECISION)/log(2))
   count<- 1
   auxscape<- 0
while(auxscape==0){
    tauM<- (tau1+tau2)/2
    CVaux<- LLR(cc,tauM,Tai)
        if(CVaux> CV){tau2<- tauM}else{tau1<- tauM}
         if(count==scape&abs(CVaux-CV)>PRECISION){
          T<- T+1.2*T; tau1<- 0; tau2<- T; CVaux<- 0; scape<- log(T/PRECISION)/log(2); count<- 1  
                                                 }else{if(abs(CVaux-CV)<=PRECISION){auxscape<- 1}else{count<- count+1} }
    
                  }
               }
                   
return(tauM)
}


##### Function for finding threshold in the scale of the Poisson couting for a fixed CV
####------------------------------------------------------------

absorb_find<- function(CV)
{

if(Tailed=="upper")
{


tauaux<- 0; cc<- 1 
taulsthes<- tals(cc,cvs.upper[cc],"upper")

while(tauaux<T){
cc<- cc+1
taulsthes<- c(taulsthes,tals(cc,cvs.upper[cc],"upper"))
tauaux<- max(taulsthes)
               }
imax<- cc-1;  mumax<- max(taulsthes) ; taulsthes<- taulsthes[-length(taulsthes)]

return(list(taulsthes,seq(1,imax),imax,mumax))
}


if(Tailed=="lower")
{

tauaux<- 0; cc<- 1 
taulsthes<- tals(cc,cvs.lower[cc],"lower")

while(tauaux<T){
cc<- cc+1
taulsthes<- c(taulsthes,tals(cc,cvs.lower[cc],"lower"))
tauaux<- max(taulsthes)
               }
imax<- cc-1;  mumax<- max(taulsthes) ; taulsthes<- taulsthes[-length(taulsthes)]

return(list(taulsthes,seq(1,imax),imax,mumax))

}


if(Tailed=="two")
{

tauaux<- 0; cc<- 1 
taulsthes1<- tals(cc,cvs.upper[cc],"upper")  
taulsthes2<- tals(cc,cvs.lower[cc],"lower")

while(tauaux<T){
cc<- cc+1
taulsthes1<- c(taulsthes1,tals(cc,cvs.upper[cc],"upper"))
taulsthes2<- c(taulsthes2,tals(cc,cvs.lower[cc],"lower"))
tauaux<- max(taulsthes2)
               }
imax<- cc-1;  taulsthes1<- taulsthes1[-length(taulsthes1)]; taulsthes2<- taulsthes2[-length(taulsthes2)]
mumax<- tauaux

absrob2<- c(rep(-1,length(taulsthes2)),seq(1,imax))
absrob1<- c(seq(1,imax),rep(-1,length(taulsthes2)))
mus<- c(taulsthes1,taulsthes2)
absorb_lower<- absrob2[order(mus)] # -1 means no test by the lower side at that time 
absorb_upper<- absrob1[order(mus)] # -1 means no test by the upper side at that time

taulsthes<- mus[order(mus)]


return(list(taulsthes,absorb_lower,absorb_upper,imax,mumax))
}

}### Close function for finding the absorb states


perf2<- function(RR,res=0)
{

ress<- res

##################################################################
######### HERE STARTS THE CODE FOR UPPER-TAILED TESTING
##################################################################

if(Tailed=="upper"){

#---------------------------------------------------------------------
# Function that calculates the product log through a recursive formula
#---------------------------------------------------------------------
ProdLog <- function(z) {
	x = z-z^2+1.5*z^3-(8/3)*z^4+(125/24)*z^5-(54/5)*z^6+(16807/720)*z^7
	for(i in 1:10) x = x-(x*exp(x)-z)/(exp(x)+x*exp(x))
	x
	}

if(Statistic=="MaxSPRT"){# 3
c = 1:max(round(2*T),1200)
CV<- cvs.upper
z<- rep(0,length(c))

for(ii in 1:max(c)){if(ii>length(CV)){z[ii] = -exp(-1-CV[length(CV)]/ii)}else{ z[ii] = -exp(-1-CV[ii]/ii)}}
mu = -c * ProdLog(z) 		#The expected counts under H0 that is needed to reject the null with i number of adverse events
                        }# close 3



if(Statistic!="MaxSPRT"){# 4
res<- absorb_find(CV)
mu<- res[[1]]
                        }# close 4
mtemp = c(0,mu)
mmu = diff(mtemp) 		#The marginal difference of the mu[] vector
RRmmu=mmu*RR			#The marginal difference of the mu[] vector
RRmu=mu*RR				#The expected counts under HA that is needed to reject the null with ii number of adverse events

imin=MinCases
while (mu[imin] < Late) imin=imin+1
if(imin>MinCases) { 
	mu[imin-1]=Late
	mmu[imin]=mu[imin]-Late
	} # end if 

imax=1						
while (mu[imax] < T) imax=imax+1 	# imax is the maximum number of cases that will generate a signal. 
	
if(imin<imax){
				              
# Defining the p[][] matrix
# -------------------------

p = seq(length=(imax-1)*imax, from=0, by=0)
dim(p) = c(imax-1,imax)


# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# When MinCases=1, there is no skipping, and it is the first row in the matrix (p[1][]).
# --------------------------------------------------------------------------------------

if(imin==MinCases) {
	for(s in 1:imin) p[imin,s]=dpois(s-1,RRmu[imin])		# Probability of having s-1 cases at time mu[MinCases]
	p[imin,imin+1]=1-ppois(imin-1,RRmu[imin])
	} # end if 

if(imin>MinCases) {
	for(s in 1:imin) p[imin-1,s]=dpois(s-1,RRmu[imin-1])		# Probability of having s-1 cases at time mu[imin-1], not rejecting H0
	p[imin-1,imin+1] = 1-ppois(imin-1,RRmu[imin-1])			# Probability of having s+ cases at time mu[imin-1], rejecting H0
	for(s in 1:imin) 								# Probability of having s-1 cases at time mu[imin], not rejectinh H0
		for(k in 1:s) 
			p[imin,s]=p[imin,s]+p[imin-1,k]*dpois(s-k,RRmmu[imin])	
	for(k in 1:imin) 
		p[imin,imin+1] = p[imin,imin+1] + p[imin-1,k]*(1-ppois(imin-k,RRmmu[imin]))
} # end if 


funcaux1<- function(ii){j<- matrix(seq(1,(ii-1)),,1); ptes<- apply(j,1,funcaux2,ii); return(ptes)}
funcaux2<- function(jj,ii){k<- seq(1,jj); return(sum(p[ii-1,k]*dpois(jj-k,RRmmu[ii])) ) }
funcaux3<- function(ii){k<- seq(1,ii-1); return(sum(p[ii-1,k]*dpois(ii-k,RRmmu[ii])) ) }
funcaux4<- function(ii){k<- seq(1,ii-1); return(sum(p[ii-1,k]*(1-ppois(ii-k,RRmmu[ii])) ) ) }

# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------

if(MinCases+1<=imax-1)
for(i in (imin+1):(imax-1)) {

p[i,1:((i-1))]<- funcaux1(i)  # This loop calculates the p[][] matix, one column at a time, from left to right
p[i,i]<- funcaux3(i)
p[i,i+1]<- funcaux4(i) # Calculates the diagonal absorbing states where H0 is rejected

                            } # end for i	
pp=0
for(k in 1:(imax-1)) pp=pp+p[imax-1,k]*(1-ppois(imax-k,(T-mu[imax-1])*RR)) #Calculates the last probability to signal before time T



# Calculating the power
# ---------------------

AlphaSpend<- rep(0,imax); AlphaSpend[imax]<- pp
power=0
if(imin>MinCases){ power=p[imin-1,imin+1] }
for(i in imin:(imax-1)){ power=power+p[i,i+1]; AlphaSpend[i]<- p[i,i+1] }		# Sums up the probabilities when a signal occurs, to get total power
power=power+pp

            


# Calculates the time until a signal occurs
#------------------------------------------

etime=0
if(imin==MinCases)
	for(n in imin:1000) etime=etime+dpois(n,RRmu[imin])*mu[imin]*imin/(n+1) 
if(imin>MinCases) {
	etime=etime+(1-ppois(imin-1,RRmu[imin-1]))*Late  						# (Late=mu[imin-1]) 
	etime=etime+mu[imin]*p[imin,imin+1]
	for(k in 1:imin) 
		for(n in 0:100)
			etime=etime+p[imin-1,k]*dpois(imin-k+1+n,RRmmu[imin])*mmu[imin]*(imin-k+1)/(imin-k+1+n+1)	# Adding expected times to signal, using a beta distribution
	} # end if


if(imin+1<=imax-1)
for(i in (imin+1):(imax-1)) {			
	etime=etime+mu[i-1]*p[i,i+1]
	for(k in 1:(i-1)) 
		for(n in 0:100)
			etime=etime+p[i-1,k]*dpois(i-k+1+n,RRmmu[i])*mmu[i]*(i-k+1)/(i-k+1+n+1)	# Adding expected times to signal, using a beta distribution
	} # end for i

etime=etime+mu[imax-1]*pp
margin=(T-mu[imax-1])

for(k in 1:(imax-1)) 
	for(n in 0:1000)
		etime=etime+p[imax-1,k]*dpois(imax-k+1+n,margin*RR)*margin*(imax-k+1)/(imax-k+1+n+1)	# Adding expected times to signal, using a beta distribution


# The expected time until signal, given that a signal occurs
# ----------------------------------------------------------

signaltime = etime/power


# The expected length of surveillance
# -----------------------------------

surveillancetime = etime+(1-power)*T


               }else{
                     power<- 1-ppois(imax-1,mu[imax]*RR)
                     
                     etime<- sum(dpois(seq(0,imax-1,1),mu[imax-1]*RR)*(mu[imax-1])*(1-ppois(imax-1-seq(0,imax-1,1),RR*(mu[imax]-mu[imax-1]))) )
                     etime<- etime + sum(dpois(seq(0,imax-1,1),mu[imax-1]*RR)*(mu[imax]-mu[imax-1])*(1-ppois(imax-seq(0,imax-1,1),RR*(mu[imax]-mu[imax-1])))/(RR*(mu[imax]-mu[imax-1])) )
                     signaltime<- etime + (1-ppois(imax-1,mu[imax-1]*RR))*mu[imax-1]
                    
                     surveillancetime = signaltime+(1-power)*mu[imax] 

                    } # end if(imin<imax)




      } ### CLOSES IF() RELATED TO THE FUNCTION FOR UPPER-TAILED TESTING

##################################################################
######### HERE CLOSES THE CODE FOR UPPER-TAILED TESTING
##################################################################




##################################################################
######### HERE STARTS THE CODE FOR TWO-TAILED OR LOWER-TAILED TESTING
################################################################## 

if(Tailed=="two"|Tailed=="lower"){

PRECISION=0.00000001

imin<- M


res<- absorb_find(CV)

mu<- res[[1]]

if(Tailed=="two"){absorb<- res[[2]]; absorb2<- res[[3]]; imax<- res[[4]] ; mumax<- res[[5]]}   # Contains the number of events needed at time mu[i] in order to reject H0
if(Tailed=="lower"){absorb<- res[[2]]; absorb2<- rep(-1,length(absorb)); imax<- res[[3]]; mumax<- res[[4]]}
if(Tailed=="upper"){absorb2<- res[[2]]; absorb<- rep(-1,length(absorb2)); imax<- res[[3]]; mumax<- res[[4]]}

## absorb: lower threshold for rection of H0 in favor of RR<1
## absorb2: upper threshold for rection of H0 in favor of RR>1

mtemp = c(0,mu)
mmu = diff(mtemp) 		#The marginal difference of the mu[] vector
RRmmu=mmu*RR			#The marginal difference of the mu[] vector
RRmu=mu*RR				#The expected counts under HA that is needed to reject the null with ii number of adverse events

Late<- 0
#imin=MinCases
#while (mu[imin] < Late) imin=imin+1
#if(imin>MinCases) { 
#	mu[imin-1]=Late
#	mmu[imin]=mu[imin]-Late
#	} # end if 



N<- length(mu)

if(imin<imax){# OPEN 1

p = matrix(0,N,imax+1)	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's


# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

if(Tailed=="two"|Tailed=="upper"){
for(s in 1:absorb2[1]) p[1,s]=dpois(s-1,RRmu[1])		# Probability of having s-1 cases at time mu[1]
p[1,absorb2[1]+1]=1-ppois(absorb2[1]-1,RRmu[1])			# probability of rejecting H0 at time mu[imin]
AlphaSpend<- rep(0,N); AlphaSpend[1]<- p[1,absorb2[1]+1]
                                 }else{
for(s in absorb[1]:(absorb[2])) p[1,s]=dpois(s-1,RRmu[1]) # Probability of having s-1 cases at time mu[1] 
AlphaSpend<- rep(0,N) 
AlphaSpend[1]<- ppois(absorb[1]-1,RRmu[1]) # probability of rejecting H0 at time mu[imin]
                                      }

# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------

if(N>1){# 1

j1_old<- 1; if(Tailed=="two"|Tailed=="upper"){j2_old<- absorb2[1]}else{j2_old<- absorb[2]+1} 

for(i in 2:N){
 
 
     j1<- max(1,absorb[i-1]+1,j1_old) 
  
      if(absorb2[i]> -1){j2<- absorb2[i]}else{ii<- i; while(absorb2[ii]== -1&ii<N){ii<- ii+1}; j2<- absorb2[ii]}; if(j2== -1){j2<- j1+1} 

	for(j in j1:j2)					# This loop calculates the p[][] matrix, one column at a time, from left to right
		for(k in j1:min(j,j2_old))
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,RRmu[i]-RRmu[i-1])	# Calculates the standard p[][] cell values
           

	 if(absorb[i]== -1){for(k in j1:j2_old){ AlphaSpend[i]=AlphaSpend[i]+p[i-1,k]*(1-ppois(j2-k,RRmu[i]-RRmu[i-1]))}}
       #if(absorb2[i]== -1){for(k in j1_old:j1){ AlphaSpend[i]=AlphaSpend[i]+p[i-1,k]*ppois(j1_old-1-(k-1),RRmu[i]-RRmu[i-1])}}
       if(absorb2[i]== -1){AlphaSpend[i]= p[i,absorb[i]]}

j2_old<- j2 ; j1_old<- j1

} # end for i	
        }# close 1

#### Calculating the last probability to signal before time T

pp=0   ;  j1<- max(1,absorb[N]+1,j1_old)
if(imax>imin&max(mu)<T){
if(absorb[i]== -1){for(k in j1:j2_old){ pp<- pp+p[N,k]*(1-ppois(j2-k,T*RR-RRmu[N]))}}
if(absorb2[i]== -1){for(k in j1_old:j1){ pp<- pp+p[N,k]*ppois(imax-1-(k-1),T*RR-RRmu[N])}}
                          } 


# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

alpha=pp
for(i in 1:N) alpha=alpha+ AlphaSpend[i]

### Separating alpha spending from lower and upper sides
AlphaSpend_lower<- AlphaSpend[absorb> -1]; mus_lower<- mu[absorb> -1]
AlphaSpend_upper<- AlphaSpend[absorb2> -1]; mus_upper<- mu[absorb2> -1]
  

# Calculating the power
# ---------------------

power=0
#if(imin>MinCases) power=p[imin-1,imin+1]
for(i in 1:N) power=power+AlphaSpend[i]			# Sums up the probabilities when a signal occurs, to get total power
power=power+pp

            


# Calculates the time until a signal occurs
#------------------------------------------
etime=0

if(Tailed=="lower"){
etime=AlphaSpend[1]*mu[1]/2
for(i in 2:length(mu)){etime<- etime + AlphaSpend[i]*(mu[i]+mu[i-1])/2}
                   }else{
                       etime=AlphaSpend_lower[1]*mus_lower[1]/2 + AlphaSpend_upper[1]*mus_upper[1]/2
                       for(i in 2:length(mus_lower)){etime<- etime + AlphaSpend_lower[i]*(mus_lower[i]+mus_lower[i-1])/2}
                       for(i in 2:length(mus_upper)){etime<- etime + AlphaSpend_upper[i]*(mus_upper[i]+mus_upper[i-1])/2}
                        }

etime<- etime+(T-mu[i])*pp

signaltime = etime/power


# The expected length of surveillance
# -----------------------------------

surveillancetime = etime+(1-power)*T



               }else{ # RELATED TO 1
                     power<- 1-ppois(imax-1,mu[imax]*RR)
                     
                     etime<- sum(dpois(seq(0,imax-1,1),mu[imax-1]*RR)*(mu[imax-1])*(1-ppois(imax-1-seq(0,imax-1,1),RR*(mu[imax]-mu[imax-1]))) )
                     etime<- etime + sum(dpois(seq(0,imax-1,1),mu[imax-1]*RR)*(mu[imax]-mu[imax-1])*(1-ppois(imax-seq(0,imax-1,1),RR*(mu[imax]-mu[imax-1])))/(RR*(mu[imax]-mu[imax-1])) )
                     signaltime<- etime + (1-ppois(imax-1,mu[imax-1]*RR))*mu[imax-1]
                    
                     surveillancetime = signaltime+(1-power)*mu[imax] 



                    } # end if(imin<imax) 


# Output assigned as a vector
# ---------------------------



#--------------------------------------------------------------####

  }################## TWO-SIDED AND LOWER-SIDED TESTS CLOSE HERE
############################



if(ress==0){
result<- c(power,signaltime,surveillancetime)
          }else{
if(Tailed=="two"){
events_lower<- absorb; if(sum(events_lower==-1)>0){events_lower[events_lower==-1]<- NA}
events_upper<- absorb2; if(sum(events_upper==-1)>0){events_upper[events_upper==-1]<- NA}
result<- list(AlphaSpend_lower,AlphaSpend_upper,events_lower,events_upper); names(result)<- c("AlphaSpend_lower","AlphaSpend_upper","events_lower","events_upper")
                 } 
if(Tailed=="lower"|Tailed=="upper"){result<- AlphaSpend}
               }


return(result)
}



####### Tests to verify the validity of the other input parameters

teste1<- 0

if(T<=0){teste1<- 1; out<- c("SampleSize must be > 0")}
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer in the range [1,100]")}
#if(teste1==0 & T>1000){teste1<- 1; out<- c("Use SampleSize<=1000")}


if(Late>T & teste1==0){teste1<- 1; out<- c("D must be <= SampleSize") }
if(Late<0 & teste1==0){teste1<- 1; out<- c("Negative values for D does not make sense. Use 0<=D<=SampleSize.") }
if(M<1 & teste1==0){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }
if(round(M)!=M){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }

# If the parameters are incorrect in any sense, the code is interrupted and an error message is informed according to the possibilies above
#------------------------------------------------------------------------------------------------------------------------------------------
if(teste1==1){stop(out,call.=FALSE)}


#####################
##################### COLOCAR A FUNCAO QUE CALCULA O ALPHA SPENDING ABAIXO

if(method=="continuous"){perf<- perf2}

if(Tailed=="two"){
Res<- perf(1,1)
AlphaSpend_lower<- as.numeric(Res$AlphaSpend_lower%*%upper.tri(matrix(0,length(Res$AlphaSpend_lower),length(Res$AlphaSpend_lower)),diag=TRUE))
AlphaSpend_upper<- as.numeric(Res$AlphaSpend_upper%*%upper.tri(matrix(0,length(Res$AlphaSpend_upper),length(Res$AlphaSpend_upper)),diag=TRUE))
events_lower<- Res$events_lower; events_upper<- Res$events_upper 
                 }else{
AlphaSpend<- perf(1,1)
AlphaSpend<- as.numeric(AlphaSpend%*%upper.tri(matrix(0,length(AlphaSpend),length(AlphaSpend)),diag=TRUE))
                      }


#### SAVING PERFORMANCE MEASURES 
R<- matrix(RR[order(RR)],length(RR),1)
Measures<- apply(R,1,perf,res=0)
Measures<- cbind(R,matrix(as.numeric(Measures),ncol=3,byrow=T)) 
colnames(Measures)<- c("RR","Power","ESignalTime","ESampleSize")

if(Tailed=="two"){res<- list(AlphaSpend_lower,AlphaSpend_upper,events_lower,events_upper,Measures); names(res)<- c("AlphaSpend_lower","AlphaSpend_upper","events_lower","events_upper","Performance") }else{res<- list(AlphaSpend,Measures); names(res)<- c("AlphaSpend","Performance")}

return(res)

}## CLOSE DE WHOLE FUNCTION

# Example

#res<- Performance.Threshold.Poisson(SampleSize=90,CV.lower=c(2.5,2.6,2.7,2.8),CV.upper=c(3,3.1,3.2,3.3),GroupSizes=c(25,20,20,25),Tailed="two",Statistic="MaxSPRT",Delta="n",RR=c(2,1.2,1,3))




