# -------------------------------------------------------------------------
# Function to perform CMaxSPRT surveillance
# -------------------------------------------------------------------------


Analyze.CondPoisson<- function(name,test,events,PersonTimeRatio,AlphaSpend="n")
{

tau<- PersonTimeRatio
k<- events
name1<- name

safedir<- getwd()
address1<- tempdir()
y<- substr(address1,1:4,4) ; i<- 2
while(i<=nchar(address1)-3&y!="Temp"&y!="TEMP"){y<- substr(address1,i:(i+3),i+3);i<- i+1}
address1<- substr(address1,1:(i+3),i+3)
setwd(address1)

if(file.exists(paste(name,"address.txt",sep=""))==FALSE){
stop(c("There is no file called"," '",name,"' yet."," You must run the function 'AnalyzeSetUp.Poisson' first."),call. =FALSE)
                                                        }

address<- as.character(read.table(paste(name,"address.txt",sep=""),sep=";")[1,1])
setwd(address)

name<- paste(name,".","txt",sep="")

if(file.exists(name)==FALSE){
stop(c("There is no file called"," '",name,"' yet."," You must run the function 'AnalyzeSetUp.Poisson' first."),call. =FALSE)
                            }


titlecheck<- paste(name1,"title.txt",sep="")
title<- read.table(titlecheck)
title<- title[1,1]
if(title==0){title<- " "}else{title<- as.character(title)}

if( sum(is.numeric(events))!=1){stop("Symbols and texts are not applicable for 'events'. It must be an integer number or zero.",call. =FALSE)}

if( sum(is.numeric(test))!=1|length(test)>1){stop("Symbols and texts are not applicable for 'test'. It must be an integer greater than zero.",call. =FALSE)}

if(test<0){stop("'test' must be an integer greater than zero.",call. =FALSE)}

if(sum(events<=0)>0){stop("The count 'events' must be an integer greater than or equal to zero.",call. =FALSE)}

if(sum(is.numeric(tau))!=1){stop("Symbols and texts are not applicable for 'tau'. It must be a number greater than zero.",call. =FALSE)}

if(sum(tau<=0)>0){stop("The entry of 'tau' must be a number greater than zero.",call. =FALSE)}

if(length(events)>1|length(tau)>1){stop("'events' and 'tau' must be single values, not vectors.",call. =FALSE)}

####
## Uploading information from previous tests
####

inputSetUp<- read.table(name)
if(inputSetUp[1,1]!=test-1){stop(c("The current test should be"," ",inputSetUp[1,1]+1,". ",
"If you do not have information about previous tests, see the user manual for more details."),call. =FALSE)}


##################################################################
## inputSetUp is a data.frame containing:
#=> inputSetUp[1,2] has the 'SampleSize' defined by "StopType" of AnalyzeSetUp.CondPoisson function, and inputSetUp[1,3] has the overall alpha level.  
#=> inputSetUp[1,8] has 'rho', which is zero for 'Wald' alpha spending.
#=> inputSetUp[1,9] has the sample size in the scale defined by "StopType" of AnalyzeSetUp.CondPoisson function. If rho>0, then inputSetUp[1,9] is settled equal to zero.
#=> inputSetUp[1,10] has the number of events in the historic data
#=> inputSetUp[1,11] has the StopType                   
#=> inputSetUp[2,1] says if the surveillance was started and, if so, when it has ocurred. 
#=> inputSetUp[3,]  has the critical values in the scale of the events 
#=> inputSetUp[4,]  has the observed number of events events, look by look, until the (test-1)th look.
#=> inputSetUp[5,]  has the actual alpha spent until the (test-1)th look.
#=> inputSetUp[6,]  has the observed tau values, test by test, until the (test-1)th look.
#=> inputSetUp[7,]  has the target alpha spending until the (test-1)th look.
#=> inputSetUp[8,]  has the critical values tau0 in the scale of the ratio Pk/V, test by test, until the (test-1)th look.

#### 

SampleSize<- inputSetUp[1,2]
alpha<- inputSetUp[1,3]
M<- inputSetUp[1,4]
start<- inputSetUp[2,1]
reject<- inputSetUp[1,7]
rho<- inputSetUp[1,8]
cc<- inputSetUp[1,10]
StopType<- inputSetUp[1,11] ; if(StopType==1){StopType<- "Cases"}else{StopType<- "Tau"}
if(test>1){CVs_old<- inputSetUp[3,1:(test-1)]; events_old<- inputSetUp[4,1:(test-1)]; tau_old<- inputSetUp[6,1:(test-1)]; target_alpha_old<- inputSetUp[7,1:(test-1)]; actual_alpha_old<- inputSetUp[5,1:(test-1)]; tau0_old<- inputSetUp[8,1:(test-1)]}else{
events_old<-0; tau_old<- 0; tau0_old=0
}



#### More checks

if( sum(is.numeric(AlphaSpend))!=1&AlphaSpend!="n"){stop("Symbols and texts are not applicable for 'AlphaSpend'. If you want to use the default, use 'n'. Otherwise,  'AlphaSpend' must be a positive number smaller than or equal to 'alpha'.",call. =FALSE)}

if( sum(length(AlphaSpend))!=1){stop("'AlphaSpend' must be a single value, not a vector.",call. =FALSE)}

if(AlphaSpend<0&AlphaSpend!="n"){stop("'AlphaSpend' must be a positive number smaller than 'alpha'.",call.=FALSE)}

if(AlphaSpend>alpha&AlphaSpend!="n"){stop(c("'AlphaSpend' must be smaller than or equal to ",alpha,"."),call. =FALSE)}


if(rho==0){
name2<- paste(name1,"alphaspend.txt",sep="")
alphaspend<- read.table(name2) #has the target alpha spending planed for future looks and settled until the (test-1)th look.
          }




##################################################################
## Setting up the target alpha spending.

if(StopType=="Cases"){CurrentSample<- sum(events_old)+k}else{CurrentSample<- sum(tau_old)+tau}

if(SampleSize<=CurrentSample){current_alpha<- alpha}else{
if(max(inputSetUp[5,])<alpha-0.00000001&inputSetUp[1,7]==0&k+sum(events_old)>=M){

if(AlphaSpend=="n"){
  if(inputSetUp[1,8]==0){ # for Wald type of alpha spending
                         if(inputSetUp[1,9]<as.numeric(k+sum(inputSetUp[4,]))){
                         current_alpha<- alpha
                                                                                   }else{current_alpha<- as.numeric(alphaspend[k+sum(inputSetUp[4,])]) }
                        }else{
                         current_alpha<- alpha*( CurrentSample /SampleSize)^rho
                             }
                   }else{
                         if(AlphaSpend<=max(inputSetUp[5,])&test>1){stop(c("For this test, 'AlphaSpend' must be selected in the (", round(max(inputSetUp[5,]),6), ",", alpha,"] interval because it has already been spent up to ",round(max(inputSetUp[5,]),6)," until the previous test."),call. =FALSE)}
                         current_alpha<- AlphaSpend
                        }

                                                                            }
                                                         }



############################################################
###### INTERNAL AUXILIARY FUNCTIONS
############################################################



####### Function to calculate cMaxSPRT
# ------------------------------------------------------------
cLLR<- function(k,cc,tal)
{
if(k/cc<=tal){return(0)}else{return(cc*log((cc*(1+tal)/(cc+k)))+k*log((k*(1+tal)/(tal*(cc+k)))))}
}

####### Function to find the thresholds in the 'tau' scale for a given cv 
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



#####
##### AUXILIAR FUNCTIONS FOR CALCULATION OF CRITICAL VALUES
#####

if(sum(events_old)+k<=20){Inference="exact"} ; if(20<sum(events_old)+k&sum(events_old)+k<=50){Inference="conservative"}; if(50<sum(events_old)+k){Inference="liberal"}

if(Inference=="exact"){

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

                      }



####<=====
if(Inference=="conservative"){

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

                             }


####<=====
if(Inference=="liberal"){

cond.ppois<- function(x,tt){if(x<0){return(0)};y<- seq(0,x); return(sum(exp(y*log(tt)+lfactorial(cc+y-1)-lfactorial(y)-lfactorial(cc-1)-(cc+y)*log(tt+1))))}
cond.dpois<- function(x,tt){if(x<0){return(0)};return(exp(x*log(tt)+lfactorial(cc+x-1)-lfactorial(x)-lfactorial(cc-1)-(cc+x)*log(tt+1)))}
cond.dpois_mult<- function(x,tt){if(x<0){return(0)}else{return(exp(x*log(tt)-lfactorial(x) ))}}
                        }







#----- Function that calculates critical values

critical_value<- function(pold,current_alpha)
{

alphas<- current_alpha-max(actual_alpha_old)
ks<- as.numeric(c(events_old,k))
tau0_old<- as.numeric(tau0_old)

mu1<- tau0_old[test-1]
mu2<- sum(ks)
tau0<- (mu1+mu2)/2
perror<- 0
scape<- 0

while(abs(perror-alphas)>0.00000001&scape==0){
   for(s in 0:(sum(ks[1:(test-1)])-1)){ perror<- perror+ (1-cond.ppois(sum(ks)-1-s,tau0-tau0_old[test-1]))*pold[s+1,1]}
    if(perror>alphas){mu2<- tau0}else{mu1<- tau0}; tau0<- (mu1+mu2)/2
    if(abs(perror-alphas)>0.00000001){perror<- 0}else{scape<- 1} 
                                     }

## Updating pold for future tests, here denoted by pf
tals<- cbind(tau0_old,tau0)

pf<- rep(0,sum(ks))
for(ki in 0:(sum(ks)-1)){
for(s in 0:min((sum(ks[1:(test-1)])-1),ki)){if(ki>0){pf[ki+1]<- pf[ki+1]+(cond.ppois(ki-s,tals[test]-tals[test-1])-cond.ppois(ki-s-1,tals[test]-tals[test-1]))*pold[s+1,1]}else{
                                                   pf[ki+1]<- cond.ppois(ki-s,tals[test]-tals[test-1])*pold[s+1,1] 
                                                                                                                                    }
                                           } # pf[s+1]: probability of having s events at time tau0
                        }

CVf<- cLLR(sum(ks),cc,tau0)

return(list(CVf,perror,pf,tau0))
}









##########################################################
###### CALCULATING CRITICAL VALUE AND ACTUAL ALPHA SPENT FOR THE CURRENT TEST
##########################################################

# Finding critical value for the 'current_alpha'

if(events+sum(inputSetUp[4,]) >= M & reject==0 & max(inputSetUp[5,])<alpha-0.00000001 ){

if(test==1|start==0){# open 1
alphas<- current_alpha

cv1<- 0
cv2<- 10
cvm<- (cv1+cv2)/2
kn<- sum(events_old+k)
tau0<- cv_tal(kn,cc,cvm)
perror<- 0
# tau value at the very first chunk of data
while(abs(perror-alphas)>0.00000001){
    perror<- 1-cond.ppois(kn-1,tau0)
    if(perror>alphas){cv1<- cvm}else{cv2<- cvm}; cvm<- (cv1+cv2)/2; tau0<- cv_tal(kn,cc,cvm)
                                     }

actualspent<- perror
CV<- cvm

## Updating pold for future tests, here denoted by p
p<- rep(0,kn)
for(s in 0:(kn-1)){if(s>0){p[s+1]<- cond.ppois(s,tau0)-cond.ppois(s-1,tau0)}else{p[s+1]<- cond.ppois(s,tau0)}} # p[s+1]: probability of having s cases at time mu1

                    }# close 1



if(test>1&start>0){# open 2
pold<- read.table(paste(name1,"p.txt",sep=""))
res<- critical_value(pold,current_alpha)
CV<- res[[1]]
actualspent<- res[[2]]+max(actual_alpha_old)
p<- res[[3]]
tau0<- res[[4]] 
                 }# close 2



# Surveillance started?
if(start==0){start<- test}

# H0 rejected?
if(test>1){ks<- as.numeric(c(events_old,k))}else{ks<- k}
llr<- cLLR(sum(ks),cc,sum(tau_old)+tau) #<= observed likelihood statistic
if(llr>=CV){reject_new<- test}else{reject_new<- 0}

                                                                                      }else{reject_new<- max(0,reject)}

if(M>events+sum(inputSetUp[4,]) |reject==1| max(inputSetUp[5,])>=alpha-0.00000001  ){actualspent<- 0; CV<- "NA"}

if(reject>0){actualspent<- max(actual_alpha_old)}
















##########################################################
###### PRINTING TABLES AND GRAPHS WITH RESULTS
##########################################################

###############
### Situation 1: SampleSize not achieved and surveillance not started because events are still smaller than M

if(start==0){ # OPEN

result<- data.frame(matrix(0,test+1,10))

colnames(result)<- c(" "," "," ", "----------","Cumulative----"," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Person-timeR","events","Person-timeR","events","LLR","target","actual","CV","Reject H0")

actualspent<- 0
CV<- "NA"
current_alpha<- 0

result[test+1,1]<- test
result[test+1,2]<- round(tau,2)
result[test+1,3]<- events
result[test+1,4]<- tau+sum(tau_old)
result[test+1,5]<- events+sum(events_old)
result[test+1,6]<- round(cLLR(events+sum(events_old),cc,tau+sum(tau_old)),2)
result[test+1,7]<- round(current_alpha,2)
result[test+1,8]<- round(actualspent,2)
result[test+1,9]<- round(CV,6)
result[test+1,10]<- paste("No")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(tau_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(tau_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i])
result[i+1,6]<- round(cLLR(sum(events_old[1:i]),cc,sum(tau_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,7]<- round(target_alpha_old[i],4)}else{result[i+1,7]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,8]<- round(actual_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
result[i+1,9]<- round(CVs_old[i],6)
result[i+1,10]<- paste("No")

                    }
          }

message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE) 
                                              message("=>    H0 cannot be rejected yet because the cumulative events is still smaller than M.",domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", Historical number of events= ",cc," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)


           } # CLOSE




###############
### Situation 2: H0 not rejected yet and sample size not achieved

if(reject==0&reject_new==0&start>0&SampleSize>=CurrentSample){# OPEN

result<- data.frame(matrix(0,test+1,10))

colnames(result)<- c(" "," "," ", "----------","Cumulative----"," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Person-timeR","events","Person-timeR","events","LLR","target","actual","CV","Reject H0")


result[test+1,1]<- test
result[test+1,2]<- round(tau,2)
result[test+1,3]<- events
result[test+1,4]<- tau+sum(tau_old)
result[test+1,5]<- events+sum(events_old)
result[test+1,6]<- round(cLLR(events+sum(events_old),cc,tau+sum(tau_old)),2)
result[test+1,7]<- round(current_alpha,4)
result[test+1,8]<- round(actualspent,4)
result[test+1,9]<- round(CV,6)
result[test+1,10]<- paste("No")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(tau_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(tau_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i]) 
result[i+1,6]<- round(cLLR(sum(events_old[1:i]),cc,sum(tau_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,7]<- round(target_alpha_old[i],4)}else{result[i+1,7]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,8]<- round(actual_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
result[i+1,9]<- round(CVs_old[i],6)
result[i+1,10]<- paste("No")

                    }
          }

message(  " ",domain = NULL, appendLF = TRUE)
message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
message("=>   Do not reject H0. Proceed to a new test as soon as you have more data.", domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", Historical number of events= ",cc," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)



# Graphic with critical values and observed test statistic

x<- seq(start,test,1) ; critical_values<- as.numeric(result[(start+1):(test+1),9]); loglike<- as.numeric(result[2:(test+1),6])

plot(seq(1,test,1),rep(max(loglike)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Log-likelihood ratio"),main="Critical value (CV) versus Observed Log-likelihood Ratio (LLR)",ylim=c(0,max(max(loglike)+1,max(critical_values+1))))
sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)
points(seq(1,test,1),loglike,col="blue",pch=20)
lines(seq(1,test,1),loglike,col="blue",lty=1)                                                                                      
points(x,critical_values,col="red",pch=20)
lines(x,critical_values,col="red",lty=2)
legend("topleft",c("Needed to reject H0 (CV)","Observed LLR"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")


                                                                                }# CLOSE





### Situation 3: SampleSize achieved with remaining alpha spending, H0 not rejected in previous tests (that is "reject==0"), and "start>0"

if(start>0&CurrentSample>=SampleSize&alpha-actualspent>0.00000001&reject==0&reject_new==0){ # OPEN

result<- data.frame(matrix(0,test+1,10))

colnames(result)<- c(" "," "," ", "----------","Cumulative----"," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Person-timeR","events","Person-timeR","events","LLR","target","actual","CV","Reject H0")


result[test+1,1]<- test
result[test+1,2]<- round(tau,2)
result[test+1,3]<- events
result[test+1,4]<- tau+sum(tau_old)
result[test+1,5]<- events+sum(events_old)
result[test+1,6]<- round(cLLR(events+sum(events_old),cc,tau+sum(tau_old)),2)
result[test+1,7]<- round(current_alpha,4)
result[test+1,8]<- round(actualspent,4)
result[test+1,9]<- round(CV,6)
result[test+1,10]<- paste("No")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(tau_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(tau_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i]) 
result[i+1,6]<- round(cLLR(sum(events_old[1:i]),cc,sum(tau_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,7]<- round(target_alpha_old[i],4)}else{result[i+1,7]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,8]<- round(actual_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
result[i+1,9]<- round(CVs_old[i],6)
result[i+1,10]<- paste("No")

                    }
          }




message(  " ",domain = NULL, appendLF = TRUE)
message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)
message(c("The upper limit on the length of surveillance has been reached."), domain = NULL, appendLF = TRUE) 
message("Then, the 'AlphaSpend' input is no longer used, and by default the target is alpha.")                                                      
message(c("You may now end the sequential analysis without rejecting H0."), domain = NULL, appendLF = TRUE)
message(c("There is still ",round(alpha-actualspent,6)," alpha to spend if you wish continue with more analyses."), domain = NULL, appendLF = TRUE)
message(" ", domain = NULL, appendLF = TRUE)  
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)                                                         
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", Historical number of events= ",cc," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)
                                                        

# Graphic with critical values and observed test statistic

x<- seq(start,test,1) ; critical_values<- as.numeric(result[(start+1):(test+1),9]); loglike<- as.numeric(result[2:(test+1),6])

plot(seq(1,test,1),rep(max(loglike)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Log-likelihood ratio"),main="Critical value (CV) versus Observed Log-likelihood Ratio (LLR)",ylim=c(0,max(max(loglike)+1,max(critical_values+1))))
sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)
points(seq(1,test,1),loglike,col="blue",pch=20)
lines(seq(1,test,1),loglike,col="blue",lty=1)                                                                                      
points(x,critical_values,col="red",pch=20)
lines(x,critical_values,col="red",lty=2)
legend("topleft",c("Needed to reject H0 (CV)","Observed LLR"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")                                              


                                                                  } # CLOSE





###############
### Situation 4: SampleSize achieved without remaining alpha spending

if(start>0&CurrentSample>=SampleSize&alpha-actualspent<=0.00000001&start>0&reject==0&reject_new==0){ # OPEN


result<- data.frame(matrix(0,test+1,10))

colnames(result)<- c(" "," "," ", "----------","Cumulative----"," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Person-timeR","events","Person-timeR","events","LLR","target","actual","CV","Reject H0")


result[test+1,1]<- test
result[test+1,2]<- round(tau,2)
result[test+1,3]<- events
result[test+1,4]<- tau+sum(tau_old)
result[test+1,5]<- events+sum(events_old)
result[test+1,6]<- round(cLLR(events+sum(events_old),cc,tau+sum(tau_old)),2)
result[test+1,7]<- round(current_alpha,4)
result[test+1,8]<- round(actualspent,4)
result[test+1,9]<- round(CV,6)
result[test+1,10]<- paste("No")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(tau_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(tau_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i]) 
result[i+1,6]<- round(cLLR(sum(events_old[1:i]),cc,sum(tau_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,7]<- round(target_alpha_old[i],4)}else{result[i+1,7]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,8]<- round(actual_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
result[i+1,9]<- round(CVs_old[i],6)
result[i+1,10]<- paste("No")

                    }
          }



message(  " ",domain = NULL, appendLF = TRUE)
message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)
message(c("The upper limit on the length of surveillance has been reached."), domain = NULL, appendLF = TRUE)                                                       
message(c("You should end the sequential analysis without rejecting H0."), domain = NULL, appendLF = TRUE)
message(c("There is no remaining alpha to spend in futures tests."), domain = NULL, appendLF = TRUE)
message(" ", domain = NULL, appendLF = TRUE)  
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)

options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", Historical number of events= ",cc," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)
                                                         

# Graphic with critical values and observed test statistic

x<- seq(start,test,1) ; critical_values<- as.numeric(result[(start+1):(test+1),9]); loglike<- as.numeric(result[2:(test+1),6])

plot(seq(1,test,1),rep(max(loglike)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Log-likelihood ratio"),main="Critical value (CV) versus Observed Log-likelihood Ratio (LLR)",ylim=c(0,max(max(loglike)+1,max(critical_values+1))))
sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)
points(seq(1,test,1),loglike,col="blue",pch=20)
lines(seq(1,test,1),loglike,col="blue",lty=1)                                                                                      
points(x,critical_values,col="red",pch=20)
lines(x,critical_values,col="red",lty=2)
legend("topleft",c("Needed to reject H0 (CV)","Observed LLR"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")                                               


                                                                                                   } # CLOSE




###############
### Situation 5: H0 rejected in the current test

if(reject==0&reject_new>0){# OPEN

result<- data.frame(matrix(0,test+1,10))

colnames(result)<- c(" "," "," ", "----------","Cumulative----"," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Person-timeR","events","Person-timeR","events","LLR","target","actual","CV","Reject H0")


result[test+1,1]<- test
result[test+1,2]<- round(tau,2)
result[test+1,3]<- events
result[test+1,4]<- tau+sum(tau_old)
result[test+1,5]<- events+sum(events_old)
result[test+1,6]<- round(cLLR(events+sum(events_old),cc,tau+sum(tau_old)),2)
result[test+1,7]<- round(current_alpha,4)
result[test+1,8]<- round(actualspent,4)
result[test+1,9]<- round(CV,6)
result[test+1,10]<- paste("Yes")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(tau_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(tau_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i]) 
result[i+1,6]<- round(cLLR(sum(events_old[1:i]),cc,sum(tau_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,7]<- round(target_alpha_old[i],4)}else{result[i+1,7]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,8]<- round(actual_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
result[i+1,9]<- round(CVs_old[i],6)
result[i+1,10]<- paste("No")

                    }
          }


message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE) 
                                              message("=>    Reject H0. No further sequential analyses are needed.",domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)

options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", Historical number of events= ",cc," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

# Graphic with critical values and observed test statistic

x<- seq(start,test,1) ; critical_values<- as.numeric(result[(start+1):(test+1),9]); loglike<- as.numeric(result[2:(test+1),6])

plot(seq(1,test,1),rep(max(loglike)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Log-likelihood ratio"),main="Critical value (CV) versus Observed Log-likelihood Ratio (LLR)",ylim=c(0,max(max(loglike)+1,max(critical_values+1))))
sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)
points(seq(1,test,1),loglike,col="blue",pch=20)
lines(seq(1,test,1),loglike,col="blue",lty=1)                                                                                      
points(x,critical_values,col="red",pch=20)
lines(x,critical_values,col="red",lty=2)
legend("topleft",c("Needed to reject H0 (CV)","Observed LLR"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")                                              


            }# CLOSE






###############
### Situation 6: H0 rejected in previous tests

if(reject>0){# OPEN

result<- data.frame(matrix(0,test+1,10))

colnames(result)<- c(" "," "," ", "----------","Cumulative----"," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Person-timeR","events","Person-timeR","events","LLR","target","actual","CV","Reject H0")


result[test+1,1]<- test
result[test+1,2]<- round(tau,2)
result[test+1,3]<- events
result[test+1,4]<- tau+sum(tau_old)
result[test+1,5]<- events+sum(events_old)
result[test+1,6]<- round(cLLR(events+sum(events_old),cc,tau+sum(tau_old)),2)
result[test+1,7]<- paste("NA")
result[test+1,8]<- paste("NA")
result[test+1,9]<- paste("NA")
result[test+1,10]<- paste("Yes")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(tau_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(tau_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i]) 
result[i+1,6]<- round(cLLR(sum(events_old[1:i]),cc,sum(tau_old[1:i])),2)

if(i>reject){result[i+1,c(7,8,9)]<- paste("NA")}else{
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,7]<- round(target_alpha_old[i],4)}else{result[i+1,7]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,8]<- round(actual_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
result[i+1,9]<- round(CVs_old[i],6)
                                                       }
if(i<reject){result[i+1,10]<- paste("No")}else{result[i+1,10]<- paste("Yes")}

                    }
          }





message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)      
          message(paste(c("=>    H0 was rejected on test"," ",reject,". ","No further sequential analyses are needed.")),domain = NULL, appendLF = TRUE)  
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)

options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", Historical number of events= ",cc," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

            }# CLOSE

############################################################
## UPDATING INFORMATION FOR FUTURE TESTES
############################################################

### For alpha spending
if(AlphaSpend!="n"&max(inputSetUp[5,])<alpha-0.00000001&inputSetUp[1,7]==0&inputSetUp[1,9]>=events+sum(inputSetUp[4,])&start==1&rho==0){
j<- events+sum(inputSetUp[4,])
while(j<=length(alphaspend)&actualspent>alphaspend[j]){alphaspend[j]<- actualspent; j<- j+1}
                                                                                                                                       }

### For decision matters
if(test> ncol(inputSetUp) ){inputSetUp<- cbind(inputSetUp,matrix(0,nrow(inputSetUp),1))}
inputSetUp[1,1]<- test
inputSetUp[2,1]<- start
if(start>0){inputSetUp[1,7]<- reject_new}
if(reject==0){inputSetUp[3,test]<- CV}
inputSetUp[4,test]<- k
inputSetUp[5,test]<- actualspent
inputSetUp[6,test]<- tau
if(reject==0){inputSetUp[7,test]<- current_alpha;inputSetUp[8,test]<- tau0}else{inputSetUp[7,test]<- alpha}




############################################################
## SAVING INFORMATION FOR FUTURE TESTES
############################################################

if(rho==0){write.table(alphaspend,paste(name1,"alphaspend.txt",sep=""))}

write.table(inputSetUp,name)

if(start>0&reject==0){write.table(p,paste(name1,"p.txt",sep=""))}

result2<- result[2:(test+1),]
colnames(result2)<- c("Test","Person-timeR","events","Cum. Person-timeR","Cum. events","LLR","target alpha","actual alpha","CV","Reject H0")
invisible(result2)

#####################################
}##### Close function Analyze.Poisson
#####################################

#AnalyzeSetUp.CondPoisson(name="TestA",SampleSizeType="Events",K=100,cc=20,alpha=0.05,M=1,AlphaSpendType="Wald",rho="n",title="n",address="C:/Users/Visitante/Ivair/POST-DOC/Material para construcao do pacote Sequential/PASTA PARA TREINO")


#Analyze.CondPoisson(name="TestA",test=1,events=5,PersonTimeRatio=0.5,AlphaSpend="n")
#Analyze.CondPoisson(name="TestA",test=2,events=6,PersonTimeRatio=0.3,AlphaSpend="n")
#Analyze.CondPoisson(name="TestA",test=3,events=10,PersonTimeRatio=0.1,AlphaSpend="n")



