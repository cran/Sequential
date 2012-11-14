CV.Poisson <-
function(SampleSize=30,D=0,M=1,alpha=0.05){
T<- SampleSize
L<- T
# ------------------- INPUT VARIABLE ----------------------------------------------------------
# T = maximum length of surveillance, defined in terms of expected counts under H0
# alpha = desired alpha level
# CVstart = CV start value, a good guess reduces computing time
# MinCases = The minimum number of cases for which a signal is allowed to occur
# Late = Time < T for first look at the data, defined in terms of the expected counts under H0

Late<- D
MinCases<- M

data(TableCV.PoissonD)
data(TableCV.PoissonM)

colnames(TableCV.PoissonD)<- c(0,1,2,3,4,6,8,10)
colnames(TableCV.PoissonM)<- c(1,2,3,4,6,8,10)

TV<- c(1,1.5,2,2.5,3,4,5,6,8,10,12,15,20,25,30,40,50,60,80,100,120,150,200,250,300,400,500,600,800,1000)
MV<- c(1,2,3,4,6,8,10)
DV<- c(0,1,2,3,4,6,8,10)

tDescont<- matrix(c(5,10,10,10,15,20,60,60,80,800,1000,1000,1,2,4,8,10,3,4,6,8,3,1,8),nrow=2,byrow=T)
# ------------------------------------------------------------

####### Tests to verify the validity of the chosen parameters

teste1<- 0

if(T<=0){teste1<- 1; out<- c("SampleSize must be > 0")}
if(teste1==0){if(alpha>0.5 | alpha<(10^(-7))){teste1<- 1; out<- c("alpha must be a number in the (1e-7,0.5] interval")}}
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer in the range [1,100]")}
#if(teste1==0 & T>1000){teste1<- 1; out<- c("Use T<=1000")}


if(teste1==0){if(1-ppois(MinCases-1,T)<alpha){
                                teste1<- 1
                                T_min<- min(seq(T,M,0.01)[1-ppois(M-1,seq(T,M,0.01))>=alpha])
                                out<- list("Does not have solution. For this M and alpha, SampleSize must be >=",T_min,".")
                                             }
             }
if(Late>T & teste1==0){teste1<- 1; out<- c("D must be <= SampleSize") }
if(Late<0 & teste1==0){teste1<- 1; out<- c("Negative values for D does not make sense. Use 0<=D<=SampleSize.") }
if(M<1 & teste1==0){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }

if(teste1==1){stop(out,call.=FALSE)}

####### end parameters validity tests
# ------------------------------------------------------------

teste2<- teste1

####### Choosing CVstart by using the known CV table to save time.
# -------------------------------------------------------------------------

x1<- 1
x2<- 0
while(x1<=ncol(tDescont)&x2==0){if(T==tDescont[1,x1]&D==tDescont[2,x1]){x2<- 1};x1<- x1+1}

if(x2==0){
if(alpha==0.05 & T>=2 & teste1==0){
               il<- sum(TV<=T)
               ir<- sum(MV<=MinCases)
               ir2<- sum(DV<=Late)
             
 if(sum(T==TV)>0 & sum(M==MV)>0){if(Late==0){teste1<- 1;out<- TableCV.PoissonM[il,ir]}else{if(M==1 & sum(Late==DV)>0 & TableCV.PoissonD[il,ir2]>0){teste1<- 1;out<- TableCV.PoissonD[il,ir2]}else{CVstart<- 3}}}else{
           il2<- 1+sum(TV<T) 
           if(Late==0){CVstart<- (TableCV.PoissonM[il2,ir]+TableCV.PoissonM[il,ir])/2}else{if((TableCV.PoissonD[il2,ir]+TableCV.PoissonD[il,ir2])/2>0){CVstart<- (TableCV.PoissonD[il2,ir]+TableCV.PoissonD[il,ir2])/2}else{CVstart<-3}}
                                                                                                               }
                                 }else{CVstart<- 3}
         }else{CVstart<- 3}

####### end the choice of CVstart by using the known CV table to save time.  
# -------------------------------------------------------------------------


if(teste1==0){ # calculating CV only if the parameters tests have not found a warning message and the known CV table does not have the corresponding CV solution. 

#---------------------------------------------------------------------
# Function that calculates the product log through a recursive formula
#---------------------------------------------------------------------
ProdLog <- function(z){
	x = z-z^2+1.5*z^3-(8/3)*z^4+(125/24)*z^5-(54/5)*z^6+(16807/720)*z^7
	for(i in 1:10) x = x-(x*exp(x)-z)/(exp(x)+x*exp(x))
	x
	                } # end ProdLog function 

#----------------------------------------------------------------------------------------------
# Function that calculates the probability of type I error for a given set of IMPUT parameters
#----------------------------------------------------------------------------------------------
Perror_I<- function(cv){

z = -exp(-1-cv/c)
mu = -c * ProdLog(z) 		#The expected counts under H0 that is needed to reject the null with i number of adverse events
mtemp = c(0,mu)
mmu = diff(mtemp) 		#The marginal difference of the mu[] vector


imin=MinCases
while (mu[imin] < Late) imin=imin+1
if(imin>MinCases) { 
	mu[imin-1]=Late
	mmu[imin]=mu[imin]-Late
	} #END if imin>MinCases

imax=1
while (mu[imax] < T) imax=imax+1    		# imax is the maximum number of cases that will generate a signal.            

# NOTE: If imax=1, this code will not work


if(imin<imax){

# Defining the p[][] matrix
# -------------------------

p = seq(length=(imax-1)*imax, from=0, by=0)				# p[i,j] is the probability of having j-1 cases at time mu[i]
dim(p) = c(imax-1,imax)								# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row p[imin][] in the matrix for which there is a chance to reject H0
# When MinCases=1, there is no skipping, and it is the first row in the matrix (p[1][]).
# ------------------------------------------------------------------------------------------

if(imin==MinCases) {
	for(s in 1:imin) p[imin,s] = dpois(s-1,mu[imin])			# Probability of having s-1 cases at time mu[imin], not rejectinh H0
	p[imin,imin+1] = 1-ppois(imin-1,mu[imin])				# Probability of having s+ cases at time mu[imin], rejectinh H0
	} # end if

if(imin>MinCases) {
	for(s in 1:imin) p[imin-1,s]=dpois(s-1,mu[imin-1])		# Probability of having s-1 cases at time mu[imin-1], not rejecting H0
	p[imin-1,imin+1] = 1-ppois(imin-1,mu[imin-1])				# Probability of having s+ cases at time mu[imin-1], rejecting H0
	for(s in 1:imin) 								# Probability of having s-1 cases at time mu[imin], not rejectinh H0
		for(k in 1:s) 
			p[imin,s]=p[imin,s]+p[imin-1,k]*dpois(s-k,mmu[imin])	
	for(k in 1:imin) 
		p[imin,imin+1] = p[imin,imin+1] + p[imin-1,k]*(1-ppois(imin-k,mmu[imin]))
} # end if 



# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------

if(MinCases+1<=imax-1&((imin+1)<=(imax-1)))
for(i in (imin+1):(imax-1)) {
	for(j in 1:(i-1))							# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:j) 
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,mmu[i])	# Calculates the standard p[][] cell values
	for(k in 1:(i-1))
		p[i,i]=p[i,i]+p[i-1,k]*dpois(i-k,mmu[i])		# Calculates the diagonal under the absorbing states, which requires a unique formula
	for(k in 1:(i-1)) 
		p[i,i+1]=p[i,i+1]+p[i-1,k]*(1-ppois(i-k,mmu[i]))# Calculates the diagonal absorbing states where H0 is rejected
} # end for i	


pp=0
if(imax>imin)
for(k in 1:(imax-1)) pp=pp+p[imax-1,k]*(1-ppois(imax-k,T-mu[imax-1])) #Calculates the last probability to signal before time T


# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

alpha=0
if(imin>MinCases) alpha=p[imin-1,imin+1]
for(i in imin:(imax-1)) alpha=alpha+p[i,i+1]					
alpha=alpha+pp


}else{alpha<- 1-ppois(imax-1,mu[imax])} # end if(imin<imax)



                      } # end Perror_I
#############################################################################################




###### Here starts the numerical procedure to find CV as a solution for a type I error equal to alpha 
#----------------------------------------------------------------------------------------------------

ALPHAFIX=alpha
PRECISION=0.00000001

LLRold=0				#Smart start values has little effect on computing time.
alphaold=1
LLR=CVstart				#Smarter start vaules reduces computing time.
alpha=0

c = 1:1200
CV1<- LLRold
CV2<- LLR
teste<- 0

loop=0
while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION) { #1
loop=loop+1	

if(loop>100 | alpha==alphaold | LLR<0 ){teste<- 1}								

if(teste==0){alpha<- Perror_I(LLR)} 

# Searching for the critical value that gives the exact alpha level. It tries first a linear interpolation, but, if it does not works, it is used the bisection method
# -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

if(teste==0){ # this if() serves to verifies if the linear interpolation has failed 

LLRnew = LLR - (alpha-ALPHAFIX)*(LLR-LLRold)/(alpha-alphaold)
alphaold=alpha
LLRold=LLR
LLR=LLRnew
            }else{ # starting the bisection method in case of the linear interpolation has failed

CV1<- 0
CV2<- 20


                  while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION){#2

                                                              LLR<- (CV1+CV2)/2
                                                              alpha<- Perror_I(LLR)                                                               
                                                              if(alpha<ALPHAFIX){CV2<- LLR}else{CV1<- LLR}  
                                                               

                                                                               }# end while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION) #2                                                                             


                 }# end else if(teste==0)       
                                                             } #end while (abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION) #1


## The outcome of this program will be the CV, or, in case of not being possible an exact solution, it will return a table with conservative and liberal critical values.

if(abs(alpha-ALPHAFIX)>PRECISION){
teste2<- 1
CV2<- ceiling((10^6)*CV2)/(10^6)
CV1<- floor((10^6)*CV1)/(10^6)
alpha1<- Perror_I(CV1)
alpha2<- Perror_I(CV2)

out<- matrix(c(CV2,alpha2,CV1,alpha1),ncol=2,byrow=F)
out<- as.data.frame(out)
colnames(out)<- c("Conservative","Liberal")
rownames(out)<- c("CV","alpha")
                                 }else{out<- LLR}

} # end if(if teste1==0)

if(teste2==0 & teste1==0 | teste2==0 & teste1==1){return(out)}
if(teste2==1 & teste1==0){
message("For these parameters there is no exact alpha level.")
message("Follows a conservative critical value:")
message("--------------------------------")
return(out[1,1])
                         }                         

}
