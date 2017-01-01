
# -------------------------------------------------------------------------
# Function to perform the unpredictable conditional Poisson MaxSPRT surveillance - Version edited at Dez-12-2016
# -------------------------------------------------------------------------

AnalyzeSetUp.CondPoisson<- function(name,SampleSizeType="Events",T="n",K="n",cc,alpha=0.05,M=1,AlphaSpendType="Wald",rho="n",title="n",address="n")
{


# SampleSizeType = "PersonTimeRatio" or "Events". With SampleSizeType="PersonTimeRatio", the upper limit on the time of surveillance
# is given in the scale of the ratio Pk/V, and the parameter is "T". With SampleSizeType="Events", the upper limit is given in the scale 
# of the number of events, and the parameter is "K".  
# T = maximum length of surveillance, defined in terms of the relative person-time observed during the surveillance period
# cc = number of adverse events observed in the historical period
# K = maximum length of surveillance, defined in terms of the number of events in the surveillance period


# cc: the number of events observed in the historical sample

# SampleSize: must be in the scale of the number of Events observed in the surveillance data. It was denoted by K in the paper of Li and Kulldorff(2009).

# Example of address: "C:/Users/Ivair/Documents"

if(address=="n"){stop(c("Please, provide a valid directory address to save the setup information."),call. =FALSE)}
                                    

if(SampleSizeType=="PersonTimeRatio"){Tal<- "PersonTimeRatio"}
pho<- rho

phoref<- rho


if(sum(SampleSizeType==c("PersonTimeRatio","Events"))==0){stop("SampleSizeType must contain a label equal to 'PersonTimeRatio' or 'Events'.",call. =FALSE)}
if(SampleSizeType=="PersonTimeRatio"){
                    if(is.numeric(T)==FALSE){stop("T must be a positive number.",call. =FALSE)}else{if(T<=0){stop("T must be a positive number.",call. =FALSE)}}
                    L<- T
                   }

if(SampleSizeType=="Events"){if(M>K){stop("M must be a positive integer in the '[0,K]' interval.",call. =FALSE)}
                    if(is.numeric(K)==FALSE){stop("K must be a positive integer.",call. =FALSE)}else{if(K<=0){stop("K must be a positive integer.",call. =FALSE)}}
                   }

if(SampleSizeType=="PersonTimeRatio"){Tal<- "PersonTimeRatio"}

if(AlphaSpendType!="Wald"&AlphaSpendType!="power-type"){stop("Set AlphaSpendType= 'Wald' or AlphaSpendType= 'power-type'.",call. =FALSE)}
if(AlphaSpendType=="power-type"&is.numeric(pho)!=TRUE){stop("Symbols and texts are not applicable for 'rho'. It must be a positive number.",call. =FALSE)}
if(pho<=0&AlphaSpendType=="power-type"){stop("rho must be greater than zero or equal to the default (rho='n')",call. =FALSE)}

if(pho=="n"){pho<- 0}

if(AlphaSpendType=="Wald"){pho<- 0}

safedir<- getwd()
if(title== "n"){title<- 0}
name1<- name

# address<- choose.dir(default = "", caption = paste("Select the folder where the file '",name,"' is going to be saved.")) ## Old version of choose dir

if(SampleSizeType=="Events"){SampleSize<- K}else{SampleSize<- T}

address1<- tempdir()
y<- substr(address1,1:4,4) ; i<- 2
while(i<=nchar(address1)-3&y!="Temp"&y!="TEMP"){y<- substr(address1,i:(i+3),i+3);i<- i+1}
address1<- substr(address1,1:(i+3),i+3)

if(paste(address)=="temp"){address<- address1}
address2<- data.frame(c(0))
address2[1,1]<- address
setwd(address1)
write.table(address2,paste(name,"address.txt",sep=""),sep=";")
setwd(address)
name<- paste(name,".","txt",sep="")
if(file.exists(name)==TRUE){
stop(c("There already exists a file called"," ",name1,".
","You may want check if some test has been performed for this monitoring before. 
If you really want to overwrite the existent file, please, go to ",address," 
to delete the file '",name,"'. Alternatively, you can delete that file by using the 
following commands: ", "setwd(","'",getwd(),"'",")","; ", "file.remove(","'",name,"'","), 
then try 'AnalyzeSetUp.Binomial' again."),call. =FALSE)
                           }
MinCases<- M

if( sum(is.numeric(alpha))!=1){stop("Symbols and texts are not applicable for 'alpha'. It must be a number in the (0,0.5) interval.",call. =FALSE)}
if( sum(is.numeric(MinCases))!=1){stop("Symbols and texts are not applicable for 'M'. It must be an integer greater than zero.",call. =FALSE)}
if(is.numeric(SampleSize)==FALSE){stop("Symbols and texts are not applicable for 'SampleSize'. It must be an integer greater than zero.",call. =FALSE)}

if(SampleSize<=0){stop("'SampleSize' must be an integer greater than zero.",call. =FALSE)}

if(alpha<=0|alpha>0.5||is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(MinCases>SampleSize||is.numeric(MinCases)==FALSE){stop("'M' must be an integer smaller than or equal to 'SampleSize'.",call. =FALSE)}
if(MinCases<1){stop("'M' must be an integer greater than zero.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("'M' must be an integer.",call. =FALSE)}


alpha1<- alpha

posi<- 2

rejt<- 0

# -------------------------------------------------------------------------
# Function produces alpha spending associated to flat critical values - continuous conditional Poisson MaxSPRT
# -------------------------------------------------------------------------

SalphafLAtcv <- function(SampleSize,alpha,MinCases) {

#########################################
### APPROACH LIBERAL  ==> Inference="liberal"
########################################## 

alphar1<- alpha
#SampleSizeType="Events"
#K<- SampleSize
Late<- 0
####### Three functions for the integrated Poisson distribution with respect to the exponential distribution
# ------------------------------------------------------------
cond.ppois<- function(x,tt){if(x<0){return(0)};y<- seq(0,x); return(sum(exp(y*log(tt)+lfactorial(cc+y-1)-lfactorial(y)-lfactorial(cc-1)-(cc+y)*log(tt+1))))}
cond.dpois<- function(x,tt){if(x<0){return(0)};return(exp(x*log(tt)+lfactorial(cc+x-1)-lfactorial(x)-lfactorial(cc-1)-(cc+x)*log(tt+1)))}
cond.dpois_mult<- function(x,tt){if(x<0){return(0)}else{return(exp(x*log(tt)-lfactorial(x) ))}}
####### Tests to verify the validity of the chosen parameters
teste1<- 0
if(T<=0){teste1<- 1; out<- c("SampleSize must be > 0")}
if(teste1==0){if(alpha>0.5 | alpha<(10^(-7))){teste1<- 1; out<- c("alpha must be a number in the (1e-7,0.5] interval")}}
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer in the range [1,100]")}
#if(teste1==0 & T>1000){teste1<- 1; out<- c("Use T<=1000")}


if(SampleSizeType=="PersonTimeRatio"){ # start here 
if(teste1==0){if(1-cond.ppois(MinCases-1,T)<alpha){
                                teste1<- 1
                                T_min<- min(seq(T,M,0.01)[1-cond.ppois(M-1,seq(T,M,0.01))>=alpha])
                                out<- list("Does not have solution. For this M and alpha, SampleSize must be >=",T_min,".")
                                                  }
             }
if(Late>T & teste1==0){teste1<- 1; out<- c("D must be <= SampleSize") }
                  }# close here
if(Late<0 & teste1==0){teste1<- 1; out<- c("Negative values for D does not make sense. Use 0<=D<=SampleSize.") }
if(M<1 & teste1==0){teste1<- 1; out<- c("M must be a positive integer in the range[1,100].") }

if(teste1==1){stop(out,call.=FALSE)}

####### end parameters validity tests
# ------------------------------------------------------------


####### Function to calculate cMaxSPRT
# ------------------------------------------------------------
cLLR<- function(k,cc,tal)
{
if(k/cc<=tal){return(0)}else{return(cc*log((cc*(1+tal)/(cc+k)))+k*log((k*(1+tal)/(tal*(cc+k)))))}
}

####### Function to find the thresholds in the 'tal' scale for a given cv 
# ------------------------------------------------------------
cv_tal<- function(k,cc,cv)
{
t1<- 0
t2<- k/cc
cvt<- 0
while(max(abs(cvt-cv),abs(t2-t1))>0.00000000001){
tm<- (t1+t2)/2; cvt<- cLLR(k,cc,tm); if(cvt>cv){t1<- tm}else{t2<- tm}
                        }
return(tm)
}

if(SampleSizeType=="Events"){ks<- matrix(seq(1,K),K,1)}

#----------------------------------------------------------------------------------------------
# Function that calculates the probability of type I error for a given set of IMPUT parameters
#----------------------------------------------------------------------------------------------
Perror_I<- function(cvv,ans=0){

mu<- cv_tal(1,cc,cvv)
i<- 1
if(SampleSizeType=="PersonTimeRatio"){while(mu[i]<T){i<- i+1;mu<- cbind(mu,cv_tal(i,cc,cvv))}; imax=i ; mu[imax]<- T}else{mu = apply(ks,1,cv_tal,cc,cvv);imax=K; T<- max(mu);L<- T}

# imax is the maximum number of Events that will generate a signal.
		
mtemp = c(0,mu)
mmu = diff(mtemp) 		#The marginal difference of the mu[] vector



imin=MinCases
while (mu[imin] < Late&imin<length(mu)) imin=imin+1
if(imin>MinCases) { 
	mu[imin-1]=Late
	mmu[imin]=mu[imin]-Late
	} #END if imin>MinCases

if(imin==length(mu)&mu[imin]<Late){stop(c("For this value of 'K', 'D' must be in the [0, ", round(mu[imin],2), ") interval."),call. =FALSE)}

   		            
# NOTE: If imax=1, this code will not work


if(imin<imax){

# Defining the p[][] matrix
# -------------------------

p<- matrix(0,imax,imax+1)
#p = seq(length=(imax-1)*imax, from=0, by=0)				# p[i,j] is the probability of having j-1 Events at time mu[i]
#dim(p) = c(imax-1,imax)								# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row p[imin][] in the matrix for which there is a chance to reject H0
# When MinCases=1, there is no skipping, and it is the first row in the matrix (p[1][]).
# ------------------------------------------------------------------------------------------

if(imin==MinCases){
	for(s in 1:imin) p[imin,s] = cond.dpois_mult(s-1,mu[imin])			# Probability of having s-1 Events at time mu[imin], not rejectinh H0
	p[imin,imin+1] = 1-cond.ppois(imin-1,mu[imin])				# Probability of having s+ Events at time mu[imin], rejecting H0
	            } # end if

if(imin>MinCases){
	for(s in 1:imin) p[imin-1,s]=cond.dpois_mult(s-1,mu[imin-1])		# Probability of having s-1 Events at time mu[imin-1], not rejecting H0
	p[imin-1,imin+1] = 1-cond.ppois(imin-1,mu[imin-1])				# Probability of having s+ Events at time mu[imin-1], rejecting H0
	for(s in 1:imin) {								            # Probability of having s-1 Events at time mu[imin], not rejectinh H0
		 for(k in 1:s){ 
			p[imin,s]=p[imin,s]+p[imin-1,k]*cond.dpois_mult(s-k,mmu[imin])	
	       for(k in 1:imin){
            partial<- 0
            for(xkm1 in (k-1):(imin)){partial<- partial + exp( log(cond.dpois_mult(xkm1-(k-1),mmu[imin])) + log(p[imin-1,k]) +     lfactorial(cc+xkm1-1) - (cc+xkm1)*log(1+mu[imin]) - lfactorial(cc-1) ) }		
            p[imin,imin+1] = p[imin,imin+1] + p[imin-1,k]*exp( lfactorial(cc+(k-1)-1) - (cc+k-1)*log(1+mu[imin-1]) - lfactorial(cc-1) )  - partial
                             }
                          }
                      }
                 } # end if 

alpharef<- 0
if(imin>MinCases){alpharef=p[imin-1,imin+1]}
alpharef<- alpharef+p[imin,imin+1]

# Calculating the remaining rows in the p[][] matrix
# -------------------------------------------------

if(MinCases+1<=imax-1&((imin+1)<=(imax-1))){
i<- (imin+1)
while(alpharef<alphar1 &i<=(imax-1)){

	for(j in 1:(i-1)){							# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:j) 
			p[i,j]=p[i,j]+p[i-1,k]*cond.dpois_mult(j-k,mmu[i])	# Calculates the standard p[][] cell values
                       }


	for(k in 1:(i-1)){
		p[i,i]=p[i,i]+p[i-1,k]*cond.dpois_mult(i-k,mmu[i])		# Calculates the diagonal under the absorbing states, which requires a unique formula
                       }
  

	for(k in 1:(i-1)){ 
            partial<- 0
            for(xkm1 in (k-1):(i-1)){partial<- partial + exp( log(cond.dpois_mult(xkm1-(k-1),mmu[i])) + log(p[i-1,k]) +     lfactorial(cc+xkm1-1) - (cc+xkm1)*log(1+mu[i]) - lfactorial(cc-1) ) }
            p[i,i+1]= p[i,i+1] + exp(log(p[i-1,k]) +lfactorial(cc+(k-1)-1) - (cc+k-1)*log(1+mu[i-1]) - lfactorial(cc-1) )  - partial   # Calculates the diagonal absorbing states where H0 is rejected
                       }
                alpharef<- alpharef + p[i,i+1]
                i<- i+1

} # end for i

}	

pp=0
if(imax>imin){
for(k in 1:(imax-1)){
partial<- 0
            for(xkm1 in (k-1):(imax-1)){partial<- partial + exp( log(cond.dpois_mult(xkm1-(k-1),mmu[imax])) + log(p[imax-1,k]) +     lfactorial(cc+xkm1-1) - (cc+xkm1)*log(1+mu[imax]) - lfactorial(cc-1) ) }
            pp<- pp + exp(log(p[imax-1,k])+ lfactorial(cc+(k-1)-1) - (cc+k-1)*log(1+mu[imax-1]) - lfactorial(cc-1) )  - partial #Calculates the last probability to signal before time T
                    }
             }

# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

Salpha<- rep(0,imax) 
alpha=0
if(imin>MinCases) alpha=p[imin-1,imin+1]
for(i in imin:(imax-1)){ alpha=alpha+p[i,i+1]; Salpha[i]<- p[i,i+1]}					
alpha=alpha+pp					
Salpha[imax]<- Salpha[imax]+pp 


}else{alpha<- 1-cond.ppois(imax-1,mu[imax])} # end if(imin<imax)

if(ans==0){return(alpha)}else{return(Salpha)}

                      } # end Perror_I
#############################################################################################




###### Here starts the numerical procedure to find CV as a solution for a type I error equal to alpha 
#----------------------------------------------------------------------------------------------------

cv1<- 0
cv2<- 10
alphar<- 0
cont<- 1
ni<- log((cv2-cv1)/0.000000001)/log(2)+2
while(abs(alphar-alpha)>0.000000001&cont<=ni){
                                   cvm<- (cv1+cv2)/2 ; alphar<- Perror_I(cvv=cvm) ; if(alphar>alpha){cv1<- cvm}else{cv2<- cvm}; cont<- cont+1
                                 }

Salpha<- Perror_I(cvv=cvm,ans=1)

} # CLOSES FUNCTION THAT OBTAINS ALPHA SPENDING FOR THE POISSON MAXSPRT 




##############################################################################################################
## HERE THE TARGET ALPHA SPENDING IS DEFINED WHEN AlphaSpendTyp=Wald
##############################################################################################################

if(AlphaSpendType=="Wald"){
sa<- SalphafLAtcv(SampleSize,alpha=alpha1,MinCases)
#sa<- resE[[1]]
#mut<- resE[[2]]
j<- length(sa)
if(sum(sa)==0){stop("Choose larger SampleSize. It is not possible to find a solution for the desired alpha with the current SampleSize choice.",call. =FALSE)}
sum_sa<- sa%*%(upper.tri(matrix(0,length(sa),length(sa)),diag=TRUE))
                          }else{j<- 0}

#############################################################################################################
##   HERE WE SAVE THE KEY CONTENT TO SETUP THE SURVEILLANCE. THE CONTENT IS SAVED IN THE MATRIX  CALLED inputSetUp 
#############################################################################################################
## inputSetUp matrix contains:
# line 1: (C11) the index for the order of the test (zero entry if we did not have a first test yet), (C12) SampleSize defined by "SampleSizeType" of AnalyzeSetUp.CondPoisson function, (C13) alpha, (C14) M, (C15) base(the line of p where the looping will start in the next test), (C16) title, (C17) reject (the index indicating if and when H0 was rejected), (C18) pho (zero if Wald is used), (C19) j (the sample size in the scale defined by "SampleSizeType" if rho=0, and j=0 otherwise), C(1,10) the number of events in the historic data, C(1,11) has the SampleSizeType 
# line 2: says if the analysis has been started or not. 0 for not started. 1 for already started.
# line 3: critical values in the scale of LLR for each test
# line 4: observed events k
# line 5: actual alpha spent
# line 6: observed PersonTimeRatio values, test by test
# line 7: has the target alpha spending actually used until the (test-1)th look.
# line 8: has the critical values tau0 in the scale of the ratio Pk/V.

if(AlphaSpendType=="Wald"){k<- length(sa)}
inputSetUp<- as.data.frame(matrix(0,8,11))

if(SampleSizeType=="Events"){SampleSizeType<-1}else{SampleSizeType<- 0}

inputSetUp[1,]<- 0
inputSetUp[1,1:11]<- c(0,SampleSize,alpha,M,1,0,0,pho,j,cc,SampleSizeType) 
inputSetUp[2,]<- 0
inputSetUp[2,1]<- 0 # says if the surveillance was started or not.
inputSetUp[3,]<- 0
inputSetUp[4,]<- 0
inputSetUp[5,]<- 0
inputSetUp[6,]<- 0
inputSetUp[7,]<- 0
inputSetUp[8,]<- 0

if(AlphaSpendType=="Wald"){alphaspend<- sum_sa} #Target alpha spending to be spent event by event.


write.table(inputSetUp,name)
titlecheck<- data.frame(matrix(0,1,1))
if(title!=0){titlecheck[1,1]<- title}
 
message(c("The parameters were successfully set at '",address,"'."),domain = NULL, appendLF = TRUE)
message(c("The temporary directory of your computer has the address of the directory where the settings information of this sequential analysis is saved.
Thus, do not clean the temporary directory before finishing this sequential analysis."),domain = NULL, appendLF = TRUE)

if(AlphaSpendType=="Wald"&phoref!="n"){message(c("The value of 'rho' is ignored, as it is not used when AlphaSpendType='Wald'."),domain = NULL, appendLF = TRUE)}

write.table(titlecheck,paste(name1,"title.txt",sep=""))

if(AlphaSpendType=="Wald"){write.table(alphaspend, paste(name1,"alphaspend.txt",sep=""))}

setwd(safedir)

} ## end function AnalyzeSetUp.CondPoisson

# AnalyzeSetUp.CondPoisson(name="TestA",SampleSizeType="PersonTimeRatio",T=2,K="n",cc=20,alpha=0.05,M=1,AlphaSpendType="Wald",rho="n",title="n",address="C:/Users/Visitante/Ivair/POST-DOC/Material para construcao do pacote Sequential/PASTA PARA TREINO")






