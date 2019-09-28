CV.CondPoisson<-
function(Inference="exact",StopType="Cases",T="n",K="n",cc,D=0,M=1,alpha=0.05,Tailed="upper"){

if(Tailed!="upper"){stop("For this version of the Sequential package, CV.CondPoisson works only for 'Tailed=upper'.",call. =FALSE)}

# ------------------- INPUT VARIABLE ----------------------------------------------------------
# Inference="exact, conservative, liberal"
# StopType = "Tal" or "Cases". With StopType="Tal", the upper limit on the time of surveillance
# is given in the scale of Pk/V, and the parameter is "T". With StopType="Cases", the upper limit is given in the scale of K  
# T = maximum length of surveillance, defined in terms of the relative person-time observed during the surveillance period
# cc = number of adverse events observed in the historical period
# K = maximum length of surveillance, defined in terms of the number of events in the surveillance period
# alpha = desired alpha level
# D = Time < T for first look at the data, defined in the person-time scale used in combination with StopType="Tal" 
# M = The minimum number of cases for which a signal is allowed to occur used in combination with StopType="Cases"

if(length(Inference)>1|length(StopType)>1|length(T)>1|length(K)>1|length(cc)>1|length(D)>1|length(M)>1|length(alpha)>1){stop("The inputs cannot be vectors.",call. =FALSE)}
if(sum(StopType==c("Tal","Cases"))==0){stop("StopType must contain a label equal to 'Tal' and 'Cases'.",call. =FALSE)}
if(sum(Inference==c("conservative","exact","liberal"))==0){stop("Inference must contain a label among 'conservative', 'exact' and 'liberal'.",call. =FALSE)}
if(StopType=="Tal"){if(D>T){stop("'D' must be a number in the '[0,T]' interval.",call. =FALSE)}
                    if(is.numeric(T)==FALSE){stop("T must be a positive number.",call. =FALSE)}else{if(T<=0){stop("T must be a positive number.",call. =FALSE)}}
                    L<- T
                   }

if(StopType=="Cases"){if(M>K){stop("M must be a positive integer in the '[0,K]' interval.",call. =FALSE)}
                    if(is.numeric(K)==FALSE){stop("K must be a positive integer.",call. =FALSE)}else{if(K<=0){stop("K must be a positive integer.",call. =FALSE)}}
                   }

#########################################
### APPROACH CONSERVATIVE  ==> Inference="conservative"
##########################################

if(Inference=="conservative"){

alphar1<- alpha

Late<- D
MinCases<- M

####### Two functions for the integrated Poisson distribution with respect to the exponential distribution
# ------------------------------------------------------------
cond.ppois<- function(x,tt){if(x<0){return(0)};y<- seq(0,x); rest<- apply(matrix(y,ncol=1),1,A,tt) ; return(sum(exp( y*log(tt)-lfactorial(y)+ log(rest) ) ) )}
cond.dpois<- function(x,tt){if(x<0){return(0)};return(exp( x*log(tt)-lfactorial(x)+ log( A(x,tt) ) ) )}
cond.dpois_mult<- function(x,tt){if(x<0){return(0)}else{return(exp(x*log(tt)-lfactorial(x) ))}}

# Another auxiliar function to be used inside term A
B<- function(x,ll,tt){
                      return(1/((tt+1)^(x-ll+1)))
                     }

A<- function(x,tt){

if(x<0){return(0)}else{
                      if(x==0){return( exp( -cc*log(tt+1)+log(B(x,0,tt)) ) )}else{  # the term A in the general expression
                                 D<- function(ll){auxB<- B(x,ll,tt);return( (-1)^sum(auxB<0)* exp(lchoose(x,ll)+lfactorial(ll+cc-1)+lfactorial(x-ll)+log(abs(auxB)) -lfactorial(cc-1)-(ll+cc)*log(tt+1) ))}
                                                       ll<- matrix(seq(0,x,1),ncol=1)
                                                       return(sum(apply(ll,1,D)))}                                                                                                                 
                                }
                                 }    


####### Tests to verify the validity of the chosen parameters

teste1<- 0

if(T<=0){teste1<- 1; out<- c("SampleSize must be > 0")}
if(teste1==0){if(alpha>0.5 | alpha<(10^(-7))){teste1<- 1; out<- c("alpha must be a number in the (1e-7,0.5] interval")}}
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer in the range [1,100]")}
#if(teste1==0 & T>1000){teste1<- 1; out<- c("Use T<=1000")}


if(StopType=="Tal"){ # start here 
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

if(StopType=="Cases"){ks<- matrix(seq(1,K),K,1)}

#----------------------------------------------------------------------------------------------
# Function that calculates the probability of type I error for a given set of IMPUT parameters
#----------------------------------------------------------------------------------------------
Perror_I<- function(cvv){

mu<- cv_tal(1,cc,cvv)
i<- 1
if(StopType=="Tal"){while(mu[i]<T){i<- i+1;mu<- cbind(mu,cv_tal(i,cc,cvv))}; imax=i ; mu[imax]<- T}else{mu = apply(ks,1,cv_tal,cc,cvv);imax=K; T<- max(mu);L<- T}

# imax is the maximum number of cases that will generate a signal.
		
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
#p = seq(length=(imax-1)*imax, from=0, by=0)				# p[i,j] is the probability of having j-1 cases at time mu[i]
#dim(p) = c(imax-1,imax)								# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row p[imin][] in the matrix for which there is a chance to reject H0
# When MinCases=1, there is no skipping, and it is the first row in the matrix (p[1][]).
# ------------------------------------------------------------------------------------------

if(imin==MinCases){
	for(s in 1:imin) p[imin,s] = cond.dpois_mult(s-1,mu[imin])			# Probability of having s-1 cases at time mu[imin], not rejectinh H0
	p[imin,imin+1] = 1-cond.ppois(imin-1,mu[imin])				# Probability of having s+ cases at time mu[imin], rejectinh H0
	            } # end if

if(imin>MinCases){
	for(s in 1:imin) p[imin-1,s]=cond.dpois_mult(s-1,mu[imin-1])		# Probability of having s-1 cases at time mu[imin-1], not rejecting H0
	p[imin-1,imin+1] = 1-cond.ppois(imin-1,mu[imin-1])				# Probability of having s+ cases at time mu[imin-1], rejecting H0
	for(s in 1:imin) {								            # Probability of having s-1 cases at time mu[imin], not rejectinh H0
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
while(alpharef<=alphar1 &i<=(imax-1)){

	for(j in 1:(i-1)){							# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:j) 
			p[i,j]=p[i,j]+p[i-1,k]*cond.dpois_mult(j-k,mmu[i])	# Calculates the standard p[][] cell values
                       }


	for(k in 1:(i-1)){
		p[i,i]=p[i,i]+p[i-1,k]*cond.dpois_mult(i-k,mmu[i])		# Calculates the diagonal under the absorbing states, which requires a unique formula
                       }
  

	for(k in 1:(i-1)){ 
            partial<- 0
            #for(xkm1 in (k-1):(i-1)){partial<- partial + exp( log(cond.dpois_mult(xkm1-(k-1),mmu[i])) + log(p[i-1,k]) +     lfactorial(cc+xkm1-1) - (cc+xkm1)*log(1+mu[i]) - lfactorial(cc-1) ) }
            for(xkm1 in (k-1):(i-1)){partial<- partial + exp( log(cond.dpois_mult(xkm1-(k-1),mmu[i])) + log(p[i-1,k]) +  log(A(xkm1,mu[i])) ) }

     #p[i,i+1]= p[i,i+1] + exp(log(p[i-1,k]) +lfactorial(cc+(k-1)-1) - (cc+k-1)*log(1+mu[i-1]) - lfactorial(cc-1) )  - partial   # Calculates the diagonal absorbing states where H0 is rejected
            p[i,i+1]= p[i,i+1] + exp(log(p[i-1,k]) + log(A(k-1,mu[i-1])) )  - partial   # Calculates the diagonal absorbing states where H0 is rejected

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

alpha=0
if(imin>MinCases) alpha=p[imin-1,imin+1]
for(i in imin:(imax-1)) alpha=alpha+p[i,i+1]					
alpha=alpha+pp					



}else{alpha<- 1-cond.ppois(imax-1,mu[imax])} # end if(imin<imax)

return(alpha)

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

if(StopType=="Cases"&T!="n"){message("T was ignored since 'StopType=Cases'")}
if(StopType=="Tal"&K!="n"){message("K was ignored since 'StopType=Tal'")}

res<- list(alphar,cvm)
names(res)<- c("Type_I_Error","cv")
return(res)                    

}  ##### Closes Inference=="conservative"


#########################################
### APPROACH EXACT  ==> Inference="exact"
########################################## 

if(Inference=="exact"){
alphar1<- alpha

Late<- D
MinCases<- M

####### Two functions for the integrated Poisson distribution with respect to the exponential distribution
# ------------------------------------------------------------
cond.ppois<- function(x,tt){if(x<0){return(0)};y<- seq(0,x); rest<- apply(matrix(y,ncol=1),1,A,tt) ; return(sum(exp( y*log(tt)-lfactorial(y)+ log(rest) ) ) )}
cond.dpois<- function(x,tt){if(x<0){return(0)};return(exp( x*log(tt)-lfactorial(x)+ log( A(x,tt) ) ) )}
cond.dpois_mult<- function(x,tt){if(x<0){return(0)}else{return(exp(x*log(tt)-lfactorial(x) ))}}

# Another auxiliar function to be used inside term A
B<- function(x,ll,tt){
                     if(x==ll){return(log(tt+1)/tt)}else{                     
                     C<- function(gg){auxmod<- sum(gg-x-ll<0); return( ((-1)^auxmod)*((-1)^(x-ll-gg))*exp(-(x-ll+1)*log(tt)-log(abs(gg-x+ll))+lchoose(x-ll,gg) + (gg-x+ll)*log(tt+1) ) - ((-1)^auxmod)*(-1)^(x-ll-gg)*exp(-(x-ll+1)*log(tt)-log(abs(gg-x+ll))+lchoose(x-ll,gg) ) )}
                     g<- matrix(seq(0,x-ll-1,1),ncol=1)
                     return( exp( -(x-ll+1)*log(tt)+ log(log(tt+1)) )+ sum(apply(g,1,C))  )
                                                        }
                     }

A<- function(x,tt){

if(x<0){return(0)}else{
                      if(x==0){return( exp( -cc*log(tt+1)+log(B(x,0,tt)) ) )}else{  # the term A in the general expression
                                 D<- function(ll){auxB<- B(x,ll,tt);return( (-1)^sum(auxB<0)* exp(lchoose(x,ll)+lfactorial(ll+cc-1)+lfactorial(x-ll)+log(abs(auxB)) -lfactorial(cc-1)-(ll+cc)*log(tt+1) ))}
                                                       ll<- matrix(seq(0,x,1),ncol=1)
                                                       return(sum(apply(ll,1,D)))}                                                                                                                 
                                }
                                 }    


####### Tests to verify the validity of the chosen parameters

teste1<- 0

if(T<=0){teste1<- 1; out<- c("SampleSize must be > 0")}
if(teste1==0){if(alpha>0.5 | alpha<(10^(-7))){teste1<- 1; out<- c("alpha must be a number in the (1e-7,0.5] interval")}}
if(teste1==0 & M>100){teste1<- 1; out<- c("M must be a positive integer in the range [1,100]")}
#if(teste1==0 & T>1000){teste1<- 1; out<- c("Use T<=1000")}

if(StopType=="Tal"){ # start here 
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

if(StopType=="Cases"){ks<- matrix(seq(1,K),K,1)}

#----------------------------------------------------------------------------------------------
# Function that calculates the probability of type I error for a given set of IMPUT parameters
#----------------------------------------------------------------------------------------------
Perror_I<- function(cvv){

mu<- cv_tal(1,cc,cvv)
i<- 1
if(StopType=="Tal"){while(mu[i]<T){i<- i+1;mu<- cbind(mu,cv_tal(i,cc,cvv))}; imax=i ; mu[imax]<- T}else{mu = apply(ks,1,cv_tal,cc,cvv);imax=K; T<- max(mu);L<- T}

# imax is the maximum number of cases that will generate a signal.
		
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
#p = seq(length=(imax-1)*imax, from=0, by=0)				# p[i,j] is the probability of having j-1 cases at time mu[i]
#dim(p) = c(imax-1,imax)								# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row p[imin][] in the matrix for which there is a chance to reject H0
# When MinCases=1, there is no skipping, and it is the first row in the matrix (p[1][]).
# ------------------------------------------------------------------------------------------

if(imin==MinCases){
	for(s in 1:imin) p[imin,s] = cond.dpois_mult(s-1,mu[imin])			# Probability of having s-1 cases at time mu[imin], not rejectinh H0
	p[imin,imin+1] = 1-cond.ppois(imin-1,mu[imin])				# Probability of having s+ cases at time mu[imin], rejectinh H0
	            } # end if

if(imin>MinCases){
	for(s in 1:imin) p[imin-1,s]=cond.dpois_mult(s-1,mu[imin-1])		# Probability of having s-1 cases at time mu[imin-1], not rejecting H0
	p[imin-1,imin+1] = 1-cond.ppois(imin-1,mu[imin-1])				# Probability of having s+ cases at time mu[imin-1], rejecting H0
	for(s in 1:imin) {								            # Probability of having s-1 cases at time mu[imin], not rejectinh H0
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

# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------

if(MinCases+1<=imax-1&((imin+1)<=(imax-1))){
i<- (imin+1)
while(alpharef<=alphar1 &i<=(imax-1)){

	for(j in 1:(i-1)){							# This loop calculates the p[][] matix, one column at a time, from left to right
		for(k in 1:j) 
			p[i,j]=p[i,j]+p[i-1,k]*cond.dpois_mult(j-k,mmu[i])	# Calculates the standard p[][] cell values
                       }


	for(k in 1:(i-1)){
		p[i,i]=p[i,i]+p[i-1,k]*cond.dpois_mult(i-k,mmu[i])		# Calculates the diagonal under the absorbing states, which requires a unique formula
                       }
  

	for(k in 1:(i-1)){ 
            partial<- 0
            #for(xkm1 in (k-1):(i-1)){partial<- partial + exp( log(cond.dpois_mult(xkm1-(k-1),mmu[i])) + log(p[i-1,k]) +     lfactorial(cc+xkm1-1) - (cc+xkm1)*log(1+mu[i]) - lfactorial(cc-1) ) }
            for(xkm1 in (k-1):(i-1)){partial<- partial + exp( log(cond.dpois_mult(xkm1-(k-1),mmu[i])) + log(p[i-1,k]) +  log(A(xkm1,mu[i])) ) }

     #p[i,i+1]= p[i,i+1] + exp(log(p[i-1,k]) +lfactorial(cc+(k-1)-1) - (cc+k-1)*log(1+mu[i-1]) - lfactorial(cc-1) )  - partial   # Calculates the diagonal absorbing states where H0 is rejected
            p[i,i+1]= p[i,i+1] + exp(log(p[i-1,k]) + log(A(k-1,mu[i-1])) )  - partial   # Calculates the diagonal absorbing states where H0 is rejected

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

alpha=0
if(imin>MinCases) alpha=p[imin-1,imin+1]
for(i in imin:(imax-1)) alpha=alpha+p[i,i+1]					
alpha=alpha+pp					



}else{alpha<- 1-cond.ppois(imax-1,mu[imax])} # end if(imin<imax)

return(alpha)

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

if(StopType=="Cases"&T!="n"){message("T was ignored since 'StopType=Cases'")}
if(StopType=="Tal"&K!="n"){message("K was ignored since 'StopType=Tal'")}


res<- list(alphar,cvm)
names(res)<- c("Type_I_Error","cv")
return(res)                    


}  #### Closes Inference="exact"


#########################################
### APPROACH LIBERAL  ==> Inference="liberal"
########################################## 

if(Inference=="liberal"){

alphar1<- alpha

Late<- D
MinCases<- M

####### Two functions for the integrated Poisson distribution with respect to the exponential distribution
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


if(StopType=="Tal"){ # start here 
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

if(StopType=="Cases"){ks<- matrix(seq(1,K),K,1)}

#----------------------------------------------------------------------------------------------
# Function that calculates the probability of type I error for a given set of IMPUT parameters
#----------------------------------------------------------------------------------------------
Perror_I<- function(cvv){

mu<- cv_tal(1,cc,cvv)
i<- 1
if(StopType=="Tal"){while(mu[i]<T){i<- i+1;mu<- cbind(mu,cv_tal(i,cc,cvv))}; imax=i ; mu[imax]<- T}else{mu = apply(ks,1,cv_tal,cc,cvv);imax=K; T<- max(mu);L<- T}

# imax is the maximum number of cases that will generate a signal.
		
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
#p = seq(length=(imax-1)*imax, from=0, by=0)				# p[i,j] is the probability of having j-1 cases at time mu[i]
#dim(p) = c(imax-1,imax)								# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row p[imin][] in the matrix for which there is a chance to reject H0
# When MinCases=1, there is no skipping, and it is the first row in the matrix (p[1][]).
# ------------------------------------------------------------------------------------------

if(imin==MinCases){
	for(s in 1:imin) p[imin,s] = cond.dpois_mult(s-1,mu[imin])			# Probability of having s-1 cases at time mu[imin], not rejectinh H0
	p[imin,imin+1] = 1-cond.ppois(imin-1,mu[imin])				# Probability of having s+ cases at time mu[imin], rejecting H0
	            } # end if

if(imin>MinCases){
	for(s in 1:imin) p[imin-1,s]=cond.dpois_mult(s-1,mu[imin-1])		# Probability of having s-1 cases at time mu[imin-1], not rejecting H0
	p[imin-1,imin+1] = 1-cond.ppois(imin-1,mu[imin-1])				# Probability of having s+ cases at time mu[imin-1], rejecting H0
	for(s in 1:imin) {								            # Probability of having s-1 cases at time mu[imin], not rejectinh H0
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

alpha=0
if(imin>MinCases) alpha=p[imin-1,imin+1]
for(i in imin:(imax-1)) alpha=alpha+p[i,i+1]					
alpha=alpha+pp					



}else{alpha<- 1-cond.ppois(imax-1,mu[imax])} # end if(imin<imax)

return(alpha)

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


if(StopType=="Cases"&T!="n"){message("T was ignored since 'StopType=Cases'")}
if(StopType=="Tal"&K!="n"){message("K was ignored since 'StopType=Tal'")}

res<- list(alphar,cvm)
names(res)<- c("Type_I_Error","cv")
return(res)                    


}  ### Closes Inference="liberal"

} # End CV.CondPoisson function


###  Example
# system.time(res<- CV.CondPoisson(Inference="exact",StopType="Cases",T=1,K=15,cc=10,D=0,M=1,alpha=0.05))




