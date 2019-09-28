

## This version has the two-tail testing. September 2019


CV.Poisson<- function(SampleSize,D=0,M=1,alpha=0.05,GroupSizes="n",Tailed="upper"){


if(length(GroupSizes)>1){auxgroup<- 1}else{

if(sum(is.numeric(GroupSizes))==0& GroupSizes!="n"){stop("'GroupSizes' must be a vector of positive integers or used with the default.",call. =FALSE)}
if(GroupSizes=="n"){auxgroup<- 2}else{auxgroup<- 1}
                                          }


#################################################################################################################################
######################### HERE IS THE PART FOR CONTINUOUS SEQUENTIAL ANALYSIS (GroupSizes=1) ####################################

if(auxgroup==2){ # OPENS THE PART FOR THE CONTINUOUS CASE




alphin<- alpha
T<- SampleSize
L<- T
# ------------------- INPUT VARIABLE ----------------------------------------------------------
# T = maximum length of surveillance, defined in terms of expected counts under H0
# alpha = desired alpha level
# CVstart = CV start value, a good guess reduces computing time
# MinCases = The minimum number of cases for which a signal is allowed to occur
# Late = Time < T for first look at the data, defined in terms of the expected counts under H0
# If Tailed="upper" (default), then the threshold is given as an upper tailed testing (H0: R<=RR0), Tailed="lower" for lower tailed (H0: R>=RR0), and Tailed="two" for two-tailed testing (H0: R=RR0).
# Important: for any choice of Tailed, the interpretation is: reject H0 if the observed log-likelihood test statistic is greater than the retuned solution.
# RR0 is the relative risk under the null hypothesis. Default is RR0=1.


if(sum(Tailed==c("upper","lower","two"))==0){stop(" 'Tailed' must be chosen among 'upper', 'lower' or 'two'.",call. =FALSE)}



Late<- D
MinCases<- M

RR0=1  # This will be implemented in the next version

##################################################################
######### HERE STARTS THE CODE FOR UPPER-TAILED TESTING
##################################################################

if(Tailed=="upper"){

##### SOME PRE-CALCULATED CV's TO SAVE EXECUTION TIME IN THE UPPER-TAILED SITUATION, GIVEN M, T, ALPHA=0.05, AND D=0.

TableCV.PoissonM<-data.frame(matrix(c(2.85393700326128,2.36663822744395,1.7742178463146,0,0,0,0
,2.96497144429574,2.57638960908667,2.15070686609869,1.68320857573985,0,0,0
,3.0469774628286,2.68935364101625,2.34967932492501,2.00015758807598,0,0,0
,3.11041922987987,2.77748334336629,2.47487306554246,2.18732770846946,0,0,0
,3.16210581650405,2.8493269885449,2.56531970873171,2.31713928988974,1.76648463065627,0,0
,3.24500386808042,2.93740985148835,2.69918199672857,2.49889219887528,2.08947296098257,1.56463589358817,0
,3.29718279728549,3.01290868171239,2.80395477091307,2.62366847225542,2.26759506082247,1.93644717718729,0
,3.34272874946839,3.08209924178235,2.87390384946275,2.69934978606631,2.40680950082531,2.09383466767681,1.74055087831132
,3.41378224992083,3.17006191246864,2.98556046816797,2.82925856210701,2.57262664682324,2.33777083613016,2.08603188356668
,3.46795158033368,3.23800870652229,3.06424796686384,2.92156132978567,2.69058633468102,2.4848340056375,2.2814406139523
,3.51174870474631,3.29055057983813,3.12525250048283,2.99310550216978,2.7814346207046,2.58938849090246,2.41540178483103
,3.56259057516111,3.35326528420839,3.19995342291952,3.07561262580597,2.87793940574268,2.71199552491449,2.55663428340283
,3.62812301587201,3.4301407625587,3.28821566346814,3.17637032177521,2.99779165442139,2.84685765103372,2.71713745193818
,3.67631995875227,3.48796061486674,3.3566769512246,3.24963354903075,3.0810513452118,2.94727030411798,2.82771061367823
,3.71576438223103,3.53415028882636,3.40671487050728,3.3071350801791,3.14780120177428,3.01963940268373,2.91122208653942
,3.77466349620378,3.60505631658695,3.48595979980634,3.39197395348475,3.24661941308076,3.13049529955633,3.03073524477903
,3.81990269391777,3.65714243579444,3.54482617261637,3.4555214311217,3.31795545593698,3.21042821525469,3.11755331874543
,3.85575466209975,3.69888480478571,3.59056669774142,3.50521969535492,3.37419401277773,3.27148598919451,3.18419630429743
,3.91085346466642,3.7624739866046,3.65993869307721,3.58089967455607,3.45808683624598,3.36288832681964,3.28403027659113
,3.95232090599911,3.81014086931301,3.71199319899898,3.63650814981748,3.52008147061398,3.43006496056302,3.35579446678864
,3.98557719058308,3.84774767167223,3.75332945022298,3.68058354716942,3.56867928724568,3.48296604520948,3.41123480524391
,4.02533791597334,3.89271466068751,3.80241164747275,3.73238627490875,3.62614993442272,3.54430757139021,3.4766545706701
,4.07482766406364,3.94892970910266,3.86276222895003,3.79683534705039,3.69651138482233,3.61982463789605,3.55679903763937
,4.11223437028032,3.990901126178,3.90806497513453,3.84484673513135,3.74875734849489,3.67570325845673,3.61551308623033
,4.14213356396756,4.02415316719417,3.94413519839946,3.88271021939305,3.79014311617218,3.71945215752241,3.66182960350214
,4.18803098306215,4.07529664413041,3.99894986501106,3.94056304984504,3.85265830465641,3.78592974865935,3.73152355787215
,4.22263170622893,4.11369239777545,4.04002057942574,3.98377804171546,3.899238807548,3.83526454689286,3.78312575686608
,4.25030952971488,4.1443167419865,4.07263799414012,4.01808951759544,3.93617480815424,3.87418262220993,3.82390844776879
,4.29282949701916,4.19116679996068,4.12255869710559,4.07046611590226,3.9922717187982,3.93336365702725,3.88559998964833
,4.32491732176988,4.22641197548099,4.16002223662034,4.10966536848407,4.03421016321423,3.97745325368572,3.93152895644531),ncol=7,byrow=T))

###### END TABLE

##### SOME PRE-CALCULATED CV's TO SAVE EXECUTION TIME, GIVEN D, T, ALPHA=0.05, AND M=1.

TableCV.PoissonD<-data.frame(matrix(c(0,0,0,0,0,0,0,0
,2.964971,1.683209,0,0,0,0,0,0
,3.046977,2.000158,0,0,0,0,0,0
,3.110419,2.187328,1.600544,0,0,0,0,0
,3.162106,2.317139,1.766485,0,0,0,0,0
,3.245004,2.498892,2.089473,1.84232,0,0,0,0
,3.297183,2.545178,2.267595,1.936447,1.611553,0,0,0
,3.342729,2.546307,2.40681,2.093835,1.921859,0,0,0
,3.413782,2.694074,2.572627,2.337771,2.211199,1.829011,0,0
,3.467952,2.799333,2.591675,2.484834,2.298373,2.087405,1.834621,0
,3.511749,2.88072,2.683713,2.589388,2.415402,2.254018,1.96566,1.755455
,3.562591,2.970411,2.794546,2.711996,2.556634,2.347591,2.203782,2.020681
,3.628123,3.082511,2.918988,2.846635,2.717137,2.542045,2.425671,2.260811
,3.67632,3.15949,3.011001,2.886783,2.827711,2.668487,2.527763,2.432668
,3.715764,3.223172,3.080629,2.963485,2.911222,2.765594,2.634068,2.553373
,3.774663,3.313966,3.186878,3.078748,3.030735,2.903286,2.789967,2.68473
,3.819903,3.381606,3.261665,3.162197,3.117553,2.99958,2.897811,2.802863
,3.855755,3.434748,3.320749,3.226113,3.162908,3.05147,2.978063,2.890933
,3.910853,3.515052,3.407923,3.321868,3.247872,3.15182,3.090356,3.019184
,3.952321,3.574091,3.47261,3.391376,3.321971,3.232345,3.155596,3.109251
,3.985577,3.620223,3.523446,3.445695,3.379278,3.294843,3.222053,3.177847
,4.025338,3.675035,3.583195,3.509028,3.446674,3.367227,3.298671,3.238461
,4.074828,3.742844,3.655984,3.587079,3.528662,3.454679,3.391821,3.336012
,4.112234,3.792978,3.710128,3.644349,3.588871,3.518954,3.459256,3.406929
,4.142134,3.832686,3.752749,3.689355,3.636272,3.568952,3.512138,3.462111
,4.222632,3.938105,3.835265,3.808087,3.760123,3.700033,3.649189,3.605012
,4.25031,3.97371,3.874183,3.847892,3.801678,3.743656,3.694832,3.652326
,4.292829,4.028089,3.933364,3.887512,3.864597,3.809685,3.763627,3.723608
,4.324917,4.047191,3.977453,3.931529,3.911308,3.858669,3.814122,3.776275),ncol=8,byrow=T))

###### END TABLE


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


if(teste1==0){if(1-ppois(MinCases-1,RR0*T)<alpha){
                                teste1<- 1
                                T_min<- min(seq(T,M,0.01)[1-ppois(M-1,seq(RR0*T,M,0.01))>=alpha])
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


if(RR0==1){ # open 1

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

if(T>1000){CVstart<- 4}

         }else{CVstart<- 40} # close 1

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
RRmmu=mmu*RR0			#The marginal difference of the mu[] vector unde the null hypothesis
RRmu=mu*RR0				#The expected counts under H0 that is needed to reject the null with ii number of adverse events

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
	for(s in 1:imin) 								# Probability of having s-1 cases at time mu[imin], not rejecting H0
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
for(k in 1:(imax-1)) pp=pp+p[imax-1,k]*(1-ppois(imax-k,(T-mu[imax-1])*RR0)) #Calculates the last probability to signal before time T


# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

alpha=0
if(imin>MinCases) alpha=p[imin-1,imin+1]
for(i in imin:(imax-1)) alpha=alpha+p[i,i+1]					
alpha=alpha+pp


}else{alpha<- 1-ppois(imax-1,RRmu[imax])} # end if(imin<imax)



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

cc<- 0
mum<- 0
while(mum<T){
cc<- cc+1
if(RR0==1){zm = -exp(-1-10/cc)}else{zm = -exp(-1-40/cc)}
mum = -cc * ProdLog(zm)
mum<- mum
            }

c = 1:cc


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
CV2<- 40


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

      } ### CLOSES IF() RELATED TO THE FUNCTION FOR LOWER AND TWO-TAILED TESTING

##################################################################
######### HERE CLOSES THE CODE FOR UPPER-TAILED TESTING
##################################################################






##################################################################
######### HERE STARTS THE CODE FOR TWO-TAILED OR LOWER-TAILED TESTING
################################################################## 

if(Tailed=="two"|Tailed=="lower"){

####### Tests to verify the validity of the chosen parameters

if(T<=0){stop("SampleSize must be > 0",call. =FALSE)}
if(alpha>0.5 | alpha<(10^(-7))){stop("alpha must be a number in the (1e-7,0.5] interval",call. =FALSE)}
if(M>100){stop("M must be a positive integer in the range [1,100]",call. =FALSE)}



if(Late>T){stop("D must be <= SampleSize",call. =FALSE) }
if(Late<0){stop("Negative values for D does not make sense. Use 0<=D<=SampleSize.",call. =FALSE) }
if(M<1){stop("M must be a positive integer in the range[1,100].",call. =FALSE) }


PRECISION=0.00000001

imin<- M


#----- MAXSPRT STATISTIC

LLR <- function(cc,uu,Tai) {

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



##### Function for finding threshold in the scale of the Poisson time index given a fixed CV
####------------------------------------------------------------

tals<- function(cc,CV,Tai)
{

if(Tai=="upper"){

tau1<- 0; tau2<- 2*RR0*T
   CVaux<- 0
while(abs(CVaux-CV)>PRECISION){
    tauM<- (tau1+tau2)/2
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


tauaux<- 0; cc<- imin 
taulsthes<- tals(cc,CV,"upper")

while(tauaux<T){
cc<- cc+1
taulsthes<- c(taulsthes,tals(cc,CV,"upper"))
tauaux<- max(taulsthes)
               }
imax<- cc-1;  mumax<- max(taulsthes) ; taulsthes<- taulsthes[-length(taulsthes)]

return(list(taulsthes,seq(1,imax),imax,mumax))
}


if(Tailed=="lower")
{

tauaux<- 0; cc<- imin 
taulsthes<- tals(cc,CV,"lower")

while(tauaux<T){
cc<- cc+1
taulsthes<- c(taulsthes,tals(cc,CV,"lower"))
tauaux<- max(taulsthes)
               }
imax<- cc-1;  mumax<- max(taulsthes) ; taulsthes<- taulsthes[-length(taulsthes)]

return(list(taulsthes,seq(1,imax),imax,mumax))

}


if(Tailed=="two")
{

tauaux<- 0; cc<- imin 
taulsthes1<- tals(cc,CV,"upper")  ## NOTE THAT DIFFERENT CV's CAN BE SPECIFIED FOR LOWER AND UPPER CASES
taulsthes2<- tals(cc,CV,"lower")

while(tauaux<T){
cc<- cc+1
taulsthes1<- c(taulsthes1,tals(cc,CV,"upper"))
taulsthes2<- c(taulsthes2,tals(cc,CV,"lower"))
tauaux<- max(taulsthes2)
               }
imax<- cc-1;  taulsthes1<- taulsthes1[-length(taulsthes1)]; taulsthes2<- taulsthes2[-length(taulsthes2)]
mumax<- tauaux

absrob2<- c(rep(-1,length(taulsthes2)),seq(imin,imax))
absrob1<- c(seq(imin,imax),rep(-1,length(taulsthes2)))
mus<- c(taulsthes1,taulsthes2)
absorb_lower<- absrob2[order(mus)] # -1 means no test by the lower side at that time 
absorb_upper<- absrob1[order(mus)] # -1 means no test by the upper side at that time

taulsthes<- mus[order(mus)]


return(list(taulsthes,absorb_lower,absorb_upper,imax,mumax))
}

}### Close function for finding the absorb states


#--------------------------------------------------------------####

Perro_I<- function(CV){

res<- absorb_find(CV)

mu<- res[[1]]
mu<- RR0*mu

if(Tailed=="two"){absorb<- res[[2]]; absorb2<- res[[3]]; imax<- res[[4]] ; mumax<- RR0*res[[5]]}   # Contains the number of events needed at time mu[i] in order to reject H0
if(Tailed=="lower"){absorb<- res[[2]]; absorb2<- rep(-1,length(absorb)); imax<- res[[3]]; mumax<- RR0*res[[4]]}
if(Tailed=="upper"){absorb2<- res[[2]]; absorb<- rep(-1,length(absorb2)); imax<- res[[3]]; mumax<- RR0*res[[4]]}

## absorb: lower threshold for rection of H0 in favor of RR<1
## absorb2: upper threshold for rection of H0 in favor of RR>1

N<- length(mu)

p = matrix(0,N,imax+1)	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's


if(N==1){# 1
if(Tailed=="two"|Tailed=="upper"){
p[1,absorb2[1]+1]=1-ppois(absorb2[1]-1,mu[1])			# probability of rejecting H0 at time mu[imin]
AlphaSpend<- rep(0,N); AlphaSpend[1]<- p[1,absorb2[1]+1]
                                 }else{
AlphaSpend<- rep(0,N) 
AlphaSpend[1]<- ppois(absorb[1]-1,mu[1]) # probability of rejecting H0 at time mu[imin]
                                      }
        }# close 1


if(N>1){# 2

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

if(Tailed=="two"|Tailed=="upper"){
for(s in 1:absorb2[1]) p[1,s]=dpois(s-1,mu[1])		# Probability of having s-1 cases at time mu[1]
p[1,absorb2[1]+1]=1-ppois(absorb2[1]-1,mu[1])			# probability of rejecting H0 at time mu[imin]
AlphaSpend<- rep(0,N); AlphaSpend[1]<- p[1,absorb2[1]+1]
                                 }else{
for(s in absorb[1]:(absorb[2])) p[1,s]=dpois(s-1,mu[1]) # Probability of having s-1 cases at time mu[1]
AlphaSpend<- rep(0,N) 
AlphaSpend[1]<- ppois(absorb[1]-1,mu[1]) # probability of rejecting H0 at time mu[imin]

                                      }

# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------

j1_old<- 1; if(Tailed=="two"|Tailed=="upper"){j2_old<- absorb2[1]}else{j2_old<- absorb[2]+1}

for(i in 2:N){
 
 
     j1<- max(1,absorb[i-1]+1,j1_old) 
  
      if(absorb2[i]> -1){j2<- absorb2[i]}else{ii<- i; while(absorb2[ii]== -1&ii<N){ii<- ii+1}; j2<- absorb2[ii]}; if(j2== -1){j2<- j1+1}

	for(j in j1:j2)					# This loop calculates the p[][] matrix, one column at a time, from left to right
		for(k in j1:min(j,j2_old))
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,mu[i]-mu[i-1])	# Calculates the standard p[][] cell values
           

	 if(absorb[i]== -1){for(k in j1:j2_old){ AlphaSpend[i]=AlphaSpend[i]+p[i-1,k]*(1-ppois(j2-k,mu[i]-mu[i-1]))}}
       #if(absorb2[i]== -1){for(k in j1_old:j1){ AlphaSpend[i]=AlphaSpend[i]+p[i-1,k]*ppois(j1_old-1-(k-1),mu[i]-mu[i-1])}}
       if(absorb2[i]== -1){AlphaSpend[i]= p[i,absorb[i]]}

j2_old<- j2 ; j1_old<- j1

} # end for i	
        }# close 2


if(N>1){ # 3
#### Calculating the last probability to signal before time T


pp=0   ;  j1<- max(1,absorb[N]+1,j1_old)
if(imax>imin&max(mu/RR0)<T){
if(absorb[i]== -1){for(k in j1:j2_old){ pp<- pp+p[N,k]*(1-ppois(j2-k,RR0*T-mu[N]))}}
if(absorb2[i]== -1){for(k in j1_old:j1){ pp<- pp+p[N,k]*ppois(imax-1-(k-1),RR0*T-mu[N])}}
                          } 


# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------


alpha=pp
for(i in 1:N) alpha=alpha+ AlphaSpend[i]

### Separating alpha spending from lower and upper sides
AlphaSpend_lower<- AlphaSpend[absorb> -1]
AlphaSpend_upper<- AlphaSpend[absorb2> -1]
       }else{alpha<- sum(AlphaSpend); AlphaSpend_lower<- AlphaSpend; AlphaSpend_upper<- 0} # close 3

#return(list(alpha,AlphaSpend_lower,AlphaSpend_upper))
return(alpha)

}# end function P_erro_I
################################################################


###################################
#### Now we apply bisection method for finding the critical value



CV1<- 0
CV2<- 10
alphaux<- 0
while(abs(CV1-CV2)>PRECISION & abs(alpha-alphaux)>PRECISION){

                                                              CV<- (CV1+CV2)/2
                                                              alphaux<- Perro_I(CV)                                                              
                                                              if(alphaux<alpha){CV2<- CV}else{CV1<- CV}                                                  

                                                             }

return(CV)

 }### CLOSES IF() RELATED TO THE FUNCTION FOR LOWER AND TWO-TAILED TESTING

##################################################################
######### HERE CLOSES THE CODE FOR TWO-TAILED OR LOWER-TAILED TESTING
##################################################################


             } # CLOSES THE PART FOR THE CONTINUOUS CASE
##########################################################################################
##########################################################################################



#################################################################################################################################
######################### HERE IS THE PART FOR GROUP SEQUENTIAL ANALYSIS (GroupSizes!=1)####################################


if(auxgroup==1){ # OPENS THE PART FOR THE GROUP CASE


alphin<- alpha 
#----------ExactPoisson.R-------------------------------------------------------------------
# Function that calculates the LLR for a given observed (c) and expected (u) number of cases
# This version enables two-tailed and lower-tailed testing since April 4, 2018
#-------------------------------------------------------------------------------------------


# T = maximum length of surveillance, defined in terms of expected counts under H0
# alpha = desired alpha level
# CVstart = LLR/CV start value, a good guess reduces computing time
# MinCases = The minimum number of cases for which a signal is allowed to occur, default=1
# GroupSizes = Expected number of events under H0 between two looks at the data
# Late = Time at first look at the data, defined in terms of the expected counts under H0, default=0
# NOTE: Only one of MinCases and Late can be different from the default value
# If Tailed="upper" (default), then the threshold is given as an upper tailed testing (H0: R<=RR0), Tailed="lower" for lower tailed (H0: R>=RR0), and Tailed="two" for two-tailed testing (H0: R=RR0).
# Important: for any choice of Tailed, the interpretation is: reject H0 if the observed log-likelihood test statistic is greater than the retuned solution.
# RR0 is the relative risk under the null hypothesis. Default is RR0=1.


CVstart=40
MinCases<- M
T<- SampleSize

RR0=1 # This will be implemented in the next version. <<<<<<<<<<

if(sum(Tailed==c("upper","lower","two"))==0){stop(" 'Tailed' must be chosen among 'upper', 'lower' or 'two'.",call. =FALSE)}


if(alpha<=0|alpha>0.5|is.numeric(alpha)==FALSE){stop("alpha must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(is.numeric(MinCases)==FALSE|MinCases!=round(MinCases)|MinCases<1){stop("M must be a positive integer.",call. =FALSE)}
if(T<=0|is.numeric(T)==FALSE){stop("SampleSize must be a positive number.",call. =FALSE)}
if(is.numeric(GroupSizes)==FALSE){stop("Symbols and texts are not applicable for GroupSizes. It must contain only positive numbers.",call. =FALSE)}
if(sum(GroupSizes<=0)>1){stop("GroupSizes must contain only positive numbers.",call. =FALSE)}
if(sum(GroupSizes)!=SampleSize){stop("The entries of GroupSizes must sum up SampleSize.",call. =FALSE)}


#####################################################
####### HERE STARTS THE CODE FOR UPPER-TAILED TESTING
#####################################################

if(Tailed=="upper"){

LLR <- function(cc,uu,Tai) {
	if(cc<=uu) x=0
	if(cc>uu) x = (uu-cc) + cc*log(cc/uu)
	x
	}






                         if(1-ppois(0,RR0*T)<alpha|qpois(1-alpha,RR0*T)<MinCases){
                            out<- LLR(MinCases,T,Tai="upper")                         
                            return(out)
                                      
                                               }else{ ## number 1
#--------------------------------------------------------------####


Perro_I<- function(CV){


absorb = seq(length=imax,from=0,by=0)		# Contains the number of events needed at time mu[i] in order to reject H0
for(i in 1:imax)
	while(LLR(absorb[i],RRmu[i],Tai="upper")<CV) 
		absorb[i]=absorb[i]+1;



imin<- 1
count<-0
while(count==0){if(absorb[imin]<MinCases){absorb[imin]<- MinCases;imin<- imin+1};if(absorb[imin]>=MinCases | imin>imax){count<- 1}}



p = seq(length=imax*(absorb[imax]+1), from=0, by=0)	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's
dim(p) = c(imax,absorb[imax]+1)				# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------

for(s in 1:absorb[1]) p[1,s]=dpois(s-1,RRmu[1])		# Probability of having s-1 cases at time mu[1]
p[1,absorb[1]+1]=1-ppois(absorb[1]-1,RRmu[1])			# probability of rejecting H0 at time mu[imin]


# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------
if(imax>1){
for(i in 2:imax) {
	for(j in 1:absorb[i])					# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:min(j,absorb[i-1]))
			p[i,j]=p[i,j]+p[i-1,k]*dpois(j-k,RRmu[i]-RRmu[i-1])	# Calculates the standard p[][] cell values
	for(k in 1:absorb[i-1]) p[i,absorb[i]+1]=p[i,absorb[i]+1]+p[i-1,k]*(1-ppois(absorb[i]-k,RRmu[i]-RRmu[i-1]))
} # end for i	
          }

# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

alpha=0
for(i in 1:imax) alpha=alpha+p[i,absorb[i]+1]					
time=0
for(i in 1:imax) time=time+sum(GroupSizes[1:i])*p[i,absorb[i]+1]
#for(i in 1:imax) time=time+i*Group*p[i,absorb[i]+1]
signaltime=time/alpha					
return(alpha)
}# end function P_erro_I
################################################################



#T=10
#alpha=0.05
#CVstart=3
#Group=1
#RR=2

ALPHAFIX=alpha
PRECISION=0.00000001

alphacont<- 0
cont<- 0
aux<- 2
CVold=0
CVnew<- CVstart			#Smart start values has little effect on computing time.
alphaold=1
CV=CVstart				#Smarter start vaules reduces computing time.
alpha=0
teste<- 0
hist<- c(CVold,alpha)

CV1<- CVold
CV2<- CV


#imax = ceiling(T/Group)
imax<- length(GroupSizes)                       # The number of tests performed, including final time T.					
#mu = seq(length=imax,from=Group,by=Group)	
mmu<- 	GroupSizes                                # An array of the expected counts at each of the tests

mu<- mmu%*%upper.tri(matrix(0,length(mmu),length(mmu)),diag=T)*1
mu<- as.vector(mu)


								# mu[i] is the commulative expected count at the i'th test
mu[imax]=T							# Sets the expected count at the last test to equal T
mtemp = c(0,mu)
mmu = diff(mtemp) 		                  #The marginal difference of the mu[] vector
RRmmu=mmu*RR0			#The marginal difference of the mu[] vector unde the null hypothesis
RRmu=mu*RR0				#The expected counts under H0 that is needed to reject the null with ii number of adverse events

loops=0
while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION) {
loops=loops+1
						# Just a counter, not really needed
if(CV>0){alpha<- Perro_I(CV)} 

if(alpha==alphaold){teste<-1}

# Estimates a new critical value in the search of the one that gives alpha=0.05
# ---------------------------------------------------------------------------------


if(teste==0){

CVnew = CV - (alpha-ALPHAFIX)*(CV-CVold)/(alpha-alphaold)
alphaold=alpha
CVold=CV
CV=CVnew

hist<- rbind(hist,c(CVold,alpha))

            }else{if(CV>0){
if(alpha>ALPHAFIX & alphaold>ALPHAFIX){ # if 1
                                       CV1<- max(CV, CVold)
                                       CV<- CV1+abs(CV-CVold)+(loops-1)*0.1
                                       CVold<- CV1
                                       alpha<- Perro_I(CV)
                                       hist<- rbind(hist,c(CV,alpha))
                                      } # end if 1
if(alpha<ALPHAFIX & alphaold<ALPHAFIX){ # if 2
                                       CV2<- min(CV, CVold)
                                       CV<- CV2-abs(CV-CVold)-(loops-1)*0.1
                                       CVold<- CV2
                                       alpha<- Perro_I(CV)
                                       hist<- rbind(hist,c(CV,alpha))
                                      } # end if2
                          } # END IF(CV>0)


if((alpha<ALPHAFIX & alphaold>ALPHAFIX) | (alpha>ALPHAFIX & alphaold<ALPHAFIX) | CV<=0){ # if 3

CV1<- max(0,min(CV,CVold))
CV2<- max(CV, CVold)

while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION){

                                                              CV<- (CV1+CV2)/2
                                                              alpha<- Perro_I(CV) 
                                                              hist<- rbind(hist,c(CV,alpha))
                                                              if(alpha<ALPHAFIX){CV2<- CV}else{CV1<- CV} 

                                                             }# end while(abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION)
                                                                               } # end if 3


                 }# end else if(teste==0)

                

} #end while (abs(CV1-CV2)>PRECISION & abs(alpha-ALPHAFIX)>PRECISION)



out<- matrix(0,1,4)
if(teste==1){
hist<- hist[-1,]
CV2<- round(min(hist[hist[,2]<ALPHAFIX,1]),10)+0.000001
CV1<- round(max(hist[hist[,2]>ALPHAFIX,1]),10)
alpha1<- min(hist[hist[,2]>ALPHAFIX,2])
alpha2<- max(hist[hist[,2]<ALPHAFIX,2])


out[1,1:2]<- c(CV1,CV2)
out[1,3:4]<- c(alpha1,alpha2)

            }else{out[1,1:2]<- c(CV,alpha)}
if(length(out)>1){out<- out[2]}
return(out)
#return(hist[-1,])
                                            } ## close number 1
#--------------------------------------------------------------####

         }#  close if related to if(Tailed=="upper")

#####################################################
####### HERE CLOSES THE CODE FOR UPPER-TAILED TESTING
#####################################################




#####################################################
####### HERE STARTS THE CODE FOR LOWER-TAILED AND TWO-TAILED TESTING
#####################################################

if(Tailed=="two"|Tailed=="lower"){# open


#----- MAXSPRT STATISTIC

LLR <- function(cc,uu,Tai) {

if(Tai=="upper"){
	if(cc<=uu) x=0
	if(cc>uu) x = (uu-cc) + cc*log(cc/uu)
                   }
if(Tai=="lower"){
	if(cc>=uu) x=0
	if(cc<uu){ if(cc>0){x = (uu-cc) + cc*log(cc/uu)}else{x= uu} }
                   }
if(Tai=="two"){	
	if(cc>0){x = (uu-cc) + cc*log(cc/uu)}else{x= uu}
                   }
	return(x)
	}
             
#--------------------------



imin<- M
PRECISION=0.00000001
mmu<- GroupSizes       # An array of the expected counts at each test

mu<- mmu%*%upper.tri(matrix(0,length(mmu),length(mmu)),diag=T)*1
mu<- as.vector(mu)
N<- length(mu)
								# mu[i] is the commulative expected count at the i'th test
mu[N]=T							# Sets the expected count at the last test equal to T
mtemp = c(0,mu)
mmu = diff(mtemp) 		                  #The marginal difference of the mu[] vector
RRmmu=mmu*RR0			#The marginal difference of the mu[] vector unde the null hypothesis
RRmu=mu*RR0				#The expected counts under H0 that is needed to reject the null with ii number of adverse events

###################################################################
Perro_I<- function(CV){

####----------------------------
#### FINDING absorb
####----------------------------

absorb<- rep(-1,N); absorb2<- rep(-1,N)

for(l in 1:length(mu)){# 1

llraux<- 0
cc<- mu[l]; while( llraux<CV& cc>=0 ){cc<- cc-1; llraux<- LLR(cc,mu[l],Tai="lower")}; absorb[l]<- cc  

if(Tailed=="two"){
llraux<- 0
cc<- mu[l]; while(llraux<CV){cc<- cc+1; llraux<- LLR(cc,mu[l],Tai="two")}; absorb2[l]<- cc 
                 }

                     }# close 1

imax<- max(absorb,absorb2)



p = matrix(0,N,max(N,imax+1))	# p[i,j] is the probability of having j-1 cases at time mu[i]
									# starting probabilities are all set to zero's


# Calculating the first row in the p[][] matrix for which there is a chance to reject H0
# --------------------------------------------------------------------------------------
AlphaSpend_lower<- rep(0,N)
AlphaSpend_upper<- rep(0,N)

if(Tailed=="two"|Tailed=="upper"){
for(s in 1:absorb2[1]) p[1,s]=dpois(s-1,RRmu[1])		# Probability of having s-1 cases at time mu[1]
#p[1,absorb2[1]+1]=1-ppois(absorb2[1]-1,RRmu[1])			# probability of rejecting H0 at time mu[imin]
AlphaSpend_upper[1]<- 1-ppois(absorb2[1]-1,RRmu[1])
                                 }else{
for(s in absorb[1]:(absorb[2])) p[1,s]=dpois(s-1,RRmu[1]) # Probability of having s-1 cases at time mu[1] 
AlphaSpend_lower[1]<- ppois(absorb[1]-1,RRmu[1]) # probability of rejecting H0 at time mu[imin] 
                                      }


# Calculating the remaining rows in the p[][] matrix
# -------------------------------------------------

if(N>1){# 1

j1_old<- 1; if(Tailed=="two"|Tailed=="upper"){j2_old<- absorb2[1]}else{j2_old<- max(1,absorb[2]+1)}

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

#### Calculating the last probability to signal before time T

pp_lower=0; pp_upper=0   ;  j1<- max(1,absorb[N]+1,j1_old)
if(imax>imin&max(mu)<T){
if(absorb2[i]> -1){for(k in j1:j2_old){ pp_upper<- pp_upper+p[N,k]*(1-ppois(j2-k,T*RR0-RRmu[N]))}}
if(absorb[i]> -1){for(k in j1_old:j1){ pp_lower<- pp_lower+p[N,k]*ppois(imax-1-(k-1),T*RR0-RRmu[N])}}
                          } 


# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

alpha=pp_lower+pp_upper
for(i in 1:N) alpha=alpha+ AlphaSpend_upper[i]+AlphaSpend_lower[i]

#return(list(alpha,AlphaSpend_lower,AlphaSpend_upper))
return(alpha)


}# end function P_erro_I
################################################################

if(1-ppois(0,RR0*T)<alpha|qpois(1-alpha,RR0*T)<MinCases){
                            out<- LLR(MinCases,T,Tai="two")                         
                            return(out)
                                                 }else{
                                        
CV1<- 0
CV2<- 40
alphaux<- 0
while(abs(CV1-CV2)>PRECISION & abs(alpha-alphaux)>PRECISION){

                                                              CV<- (CV1+CV2)/2
                                                              alphaux<- Perro_I(CV)                                                              
                                                              if(alphaux<alpha){CV2<- CV}else{CV1<- CV}                                                  

                                                             }

return(CV)


                                                      }

                                 }# close if(Tailed=="two"|Tailed=="lower")

#####################################################
####### HERE CLOSES THE CODE FOR LOWER-TAILED AND TWO-TAILED TESTING
#####################################################


             } # CLOSES THE PART FOR THE GROUP CASE
##########################################################################################
##########################################################################################
                           

} ### CLOSES THE FUNCTION CV.Poisson
##################################################################

# Example
#system.time(res<- CV.Poisson(SampleSize=501,D=0,M=1,alpha=0.05,Tailed="upper",RR0=1))




