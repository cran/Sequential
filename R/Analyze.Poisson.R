# -------------------------------------------------------------------------
# Function to perform the unpredictable binomial MaxSPRT surveillance - Version edited at Jan-14-2015
# -------------------------------------------------------------------------


Analyze.Poisson<- function(name,test,mu0,events,AlphaSpend="n")
{

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

if(sum(events<0)>0){stop("The count 'events' must be an integer greater than or equal to zero.",call. =FALSE)}

if(sum(is.numeric(mu0))!=1){stop("Symbols and texts are not applicable for 'mu0'. It must be a number greater than zero.",call. =FALSE)}

if(sum(mu0<=0)>0){stop("The entry of 'mu0' must be a number greater than zero.",call. =FALSE)}

if(length(events)>1|length(mu0)>1){stop("'events' and 'mu0' must be single values, not vectors.",call. =FALSE)}

####
## Uploading information from previous tests
####

inputSetUp<- read.table(name)
if(inputSetUp[1,1]!=test-1){stop(c("The current test should be"," ",inputSetUp[1,1]+1,". ",
"If you do not have information about previous tests, see the user manual for more details."),call. =FALSE)}


##################################################################
## inputSetUp is a data.frame containing:
#=> inputSetUp[1,2] has the 'SampleSize', and inputSetUp[1,3] has the overall alpha level.  
#=> inputSetUp[1,8] has 'rho', which is zero for 'Wald' alpha spending.
#=> inputSetUp[1,9] has the sample size in the scale of the events. If rho>0, then inputSetUp[1,9] is settled equal to zero.
#=> inputSetUp[2,1] says if the surveillance was started and, if so, when it has ocurred. 
#=> inputSetUp[3,]  has the critical values in the scale of the events 
#=> inputSetUp[4,]  has the observed number of events, look by look, until the (test-1)th look.
#=> inputSetUp[5,]  has the actual alpha spent until the (test-1)th look.
#=> inputSetUp[6,]  has the expected number of events under H0, mu0, look by look, until the (test-1)th look.
#=> inputSetUp[7,]  has the target alpha spending until the (test-1)th look.


#### 

SampleSize<- inputSetUp[1,2]
alpha<- inputSetUp[1,3]
M<- inputSetUp[1,4]
start<- inputSetUp[2,1]
reject<- inputSetUp[1,7]
rho<- inputSetUp[1,8]
if(test>1){CVs_old<- inputSetUp[3,1:(test-1)]; events_old<- inputSetUp[4,1:(test-1)]; mu0_old<- inputSetUp[6,1:(test-1)]; target_alpha_old<- inputSetUp[7,1:(test-1)]; actual_alpha_old<- inputSetUp[5,1:(test-1)]}else{
events_old<-0; mu0_old<- 0
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

if(SampleSize<=mu0+sum(inputSetUp[6,1:(test-1)])){current_alpha<- alpha}else{
if(max(inputSetUp[5,])<alpha-0.00000001&inputSetUp[1,7]==0&events+sum(events_old)>=M){

if(AlphaSpend=="n"){
  if(inputSetUp[1,8]==0){ # for Wald type of alpha spending
                         if(inputSetUp[1,9]<as.numeric(events+sum(inputSetUp[4,]))){
                         current_alpha<- alpha
                                                                                   }else{current_alpha<- as.numeric(alphaspend[events+sum(inputSetUp[4,])]) }
                        }else{
                         current_alpha<- alpha*( (mu0+sum(inputSetUp[6,]) ) /inputSetUp[1,2])^rho
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


#----- THE MAXSPRT STATISTIC

LLR <- function(cc,uu) {
	if(cc<=uu) x=0
	if(cc>uu) x = (uu-cc) + cc*log(cc/uu)
	x
	}
#--------------------------





#----- Function that calculates critical values

critical_value<- function(pold,current_alpha,CVold)
{

alphas<- current_alpha-max(actual_alpha_old)

CV1<- CVold-1 ; CV2<- CV1+ qpois(1-alphas,mu0)+1
alphat<- 1
count<- 0
while(abs(alphat-alphas)>10^(-6)&count<30){
count<- count+1
CVm<- ceiling((CV1+CV2)/2)

p<- rep(0,CVm+1) # p[y+1] is the probability of having y events in the current test

for(y in 0:(CVm-1)){for(x in 0:min(CVold-1,y)){p[y+1]<- p[y+1]+dpois(y-x,mu0)*pold[x+1,1]}}

for(xx in 0:(CVold-1)){p[CVm+1]<- p[CVm+1]+ (1-ppois(CVm-1-xx,mu0))*pold[xx+1,1]}
alphat<- p[CVm+1]
if(p[CVm+1]>alphas){CV1<- CVm}else{CV2<- CVm; alphahere<- p[CVm+1];CVf<- CVm; pf<- p}

                                  }
                
return(list(CVf,alphahere,pf))
}
#--------------------------




##########################################################
###### CALCULATING CRITICAL VALUE AND ACTUAL ALPHA SPENT FOR THE CURRENT TEST
##########################################################

# Finding critical value for the 'current_alpha'

if(events+sum(inputSetUp[4,]) >= M & reject==0 & max(inputSetUp[5,])<alpha-0.00000001 ){

if(test==1){
alphas<- current_alpha
CV<- qpois(1-alphas,mu0)+1
p<- rep(0,CV)
for(x in 0:(CV-1)){p[x+1]<- dpois(x,mu0)}
actualspent<- 1-ppois(CV-1,mu0)
           }

if(test>1&start==0){
mu0h<- sum(mu0_old)+mu0
alphas<- current_alpha
CV<- qpois(1-alphas,mu0h)+1
p<- rep(0,CV)
for(x in 0:(CV-1)){p[x+1]<- dpois(x,mu0)}
actualspent<- 1-ppois(CV-1,mu0h)
                   }

if(test>1&start>0){
pold<- read.table(paste(name1,"p.txt",sep=""))
CVold<- as.numeric(CVs_old[test-1])

res<- critical_value(pold,current_alpha,CVold)

CV<- res[[1]]
actualspent<- res[[2]]+max(actual_alpha_old)
p<- res[[3]]

                 }

# Surveillance started?
if(start==0){start<- test}

# H0 rejected?
if(CV<=events+sum(inputSetUp[4,])){reject_new<- test}else{reject_new<- 0}

                                                                                      }else{reject_new<- max(0,reject)}

if(M>events+sum(inputSetUp[4,]) |reject==1| max(inputSetUp[5,])>=alpha-0.00000001  ){actualspent<- 0; CV<- "NA"}

if(reject>0){actualspent<- max(actual_alpha_old)}

##########################################################
###### PRINTING TABLES AND GRAPHS WITH RESULTS
##########################################################

###############
### Situation 1: SampleSize not achieved and surveillance not started because events are still smaller than M

if(start==0){ # OPEN

result<- data.frame(matrix(0,test+1,11))

colnames(result)<- c(" "," "," ", "-----","Cumulative----"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","mu0","Events","mu0","Events","RR estimate","LLR","target","actual","CV","Reject H0")

actualspent<- 0
CV<- "NA"
current_alpha<- 0

result[test+1,1]<- test
result[test+1,2]<- round(mu0,2)
result[test+1,3]<- events
result[test+1,4]<- mu0+sum(mu0_old)
result[test+1,5]<- events+sum(events_old)
if(events+sum(events_old)>=mu0+sum(mu0_old)){result[test+1,6]<- round((events+sum(events_old))/(mu0+sum(mu0_old)),2) }else{result[test+1,6]<- 1}
result[test+1,7]<- round(LLR(events+sum(events_old),mu0+sum(mu0_old)),2)
result[test+1,8]<- round(current_alpha,2)
result[test+1,9]<- round(actualspent,2)
result[test+1,10]<- CV
result[test+1,11]<- paste("No")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(mu0_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(mu0_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i])
if(sum(events_old[1:i])>=sum(mu0_old[1:i])){result[i+1,6]<- round(sum(events_old[1:i])/sum(mu0_old[1:i]),2) }else{result[i+1,6]<- 1} 
result[i+1,7]<- round(LLR(sum(events_old[1:i]),sum(mu0_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,8]<- round(target_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,9]<- round(actual_alpha_old[i],4)}else{result[i+1,9]<- "NA"}
result[i+1,10]<- CVs_old[i]
result[i+1,11]<- paste("No")

                    }
          }

message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE) 
                                              message("=>    H0 cannot be rejected yet because the cumulative events is still smaller than M.",domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", mu0= ",mu0," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)


           } # CLOSE




###############
### Situation 2: H0 not rejected yet and sample size not achieved

if(reject==0&reject_new==0&start>0&mu0+sum(mu0_old)<SampleSize){# OPEN

result<- data.frame(matrix(0,test+1,11))

colnames(result)<- c(" "," "," ", "----","Cumulative----"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","mu0","Events","mu0","Events","RR estimate","LLR","target","actual","CV","Reject H0")


result[test+1,1]<- test
result[test+1,2]<- round(mu0,2)
result[test+1,3]<- events
result[test+1,4]<- round(mu0+sum(mu0_old),2)
result[test+1,5]<- events+sum(events_old)
if(events+sum(events_old)>=mu0+sum(mu0_old)){result[test+1,6]<- round((events+sum(events_old))/(mu0+sum(mu0_old)),2) }else{result[test+1,6]<- 1}
result[test+1,7]<- round(LLR(events+sum(events_old),mu0+sum(mu0_old)),2)
result[test+1,8]<- round(current_alpha,4)
result[test+1,9]<- round(actualspent,4)
result[test+1,10]<- CV
result[test+1,11]<- paste("No")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(mu0_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(mu0_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i])
if(sum(events_old[1:i])>=sum(mu0_old[1:i])){result[i+1,6]<- round(sum(events_old[1:i])/sum(mu0_old[1:i]),2) }else{result[i+1,6]<- 1} 
result[i+1,7]<- round(LLR(sum(events_old[1:i]),sum(mu0_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,8]<- round(target_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,9]<- round(actual_alpha_old[i],4)}else{result[i+1,9]<- "NA"}
result[i+1,10]<- CVs_old[i]
result[i+1,11]<- paste("No")

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
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", mu0= ",mu0," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)


x<- seq(start,test,1) ; observed<- as.numeric(result[2:(test+1),5]) ; critical_values<- as.numeric(result[(start+1):(test+1),10])
target<- as.numeric(result[(start+1):(test+1),8]) ; actual<- as.numeric(result[(start+1):(test+1),9])
loglike<- as.numeric(result[2:(test+1),7]); RRest<- as.numeric(result[2:(test+1),6])

#########>>>>>>>>>>>>>> critical_values in the scale of MaxSPRT

MaxSPRT_critical_values<- rep(0,length(critical_values))
if(test==1){mu0h<- mu0}else{mu0h<- c(as.numeric(mu0_old),mu0)}
for(i in 1:length(critical_values)){MaxSPRT_critical_values[i]<- LLR(critical_values[i],sum(mu0h[1:i]))}


# Graphic 1
par(mfrow=c(2,2))
plot(seq(1,test,1),rep(max(observed,critical_values),test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Cumulative Events"),main="Number of events scale",sub=title,ylim=c(0,5+max(observed,critical_values)))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),observed,col="blue",pch=20)
lines(seq(1,test,1),observed,col="blue",lty=1)                                              
points(x,critical_values,col="red",pch=20)
lines(x,critical_values,col="red",lty=2)
legend("topleft",c("Needed to reject H0 (CV)","Observed"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

# Graphic 2
plot(seq(1,test,1),rep(alpha+0.02,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Alpha spending"),main="Alpha Spending",ylim=c(0,alpha+0.02))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(x,actual,col="blue",pch=20)
lines(x,actual,col="blue",lty=1)                                              
points(x,target,col="red",pch=20)
lines(x,target,col="red",lty=2)
legend("topleft",c("Target","Actual spent"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

# Graphic 3
plot(seq(1,test,1),rep(max(RRest)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Observed relative risk"),main="Observed Relative Risk",ylim=c(0,max(RRest)+1))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),RRest,col="blue",pch=20)
lines(seq(1,test,1),RRest,col="blue",lty=1)                                              

# Graphic 4
plot(seq(1,test,1),rep(max(loglike)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Log-likelihood ratio"),main="MaxSPRT scale",ylim=c(0,max(MaxSPRT_critical_values)+5))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),loglike,col="blue",pch=20)
lines(seq(1,test,1),loglike,col="blue",lty=1)  

points(seq(start,test,1),MaxSPRT_critical_values,col="red",pch=20)
lines(seq(start,test,1),MaxSPRT_critical_values,col="red",lty=2)

legend("topleft",c("Needed to reject H0 (CV)","Observed"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")



                                                                                }# CLOSE




###############
### Situation 3: SampleSize achieved with remaining alpha spending, H0 not rejected in previous tests (that is "reject==0"), and "start>0"

if(start>0&mu0+sum(mu0_old)>=SampleSize&alpha-actualspent>0.00000001&reject==0&reject_new==0){ # OPEN

result<- data.frame(matrix(0,test+1,11))

colnames(result)<- c(" "," "," ", "-----","Cumulative----"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","mu0","Events","mu0","Events","RR estimate","LLR","target","actual","CV","Reject H0")


result[test+1,1]<- test
result[test+1,2]<- round(mu0,2)
result[test+1,3]<- events
result[test+1,4]<- round(mu0+sum(mu0_old),2)
result[test+1,5]<- events+sum(events_old)
if(events+sum(events_old)>=mu0+sum(mu0_old)){result[test+1,6]<- round((events+sum(events_old))/(mu0+sum(mu0_old)),2) }else{result[test+1,6]<- 1}
result[test+1,7]<- round(LLR(events+sum(events_old),mu0+sum(mu0_old)),2)
result[test+1,8]<- round(current_alpha,4)
result[test+1,9]<- round(actualspent,4)
result[test+1,10]<- CV
result[test+1,11]<- paste("No")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(mu0_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(mu0_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i])
if(sum(events_old[1:i])>=sum(mu0_old[1:i])){result[i+1,6]<- round(sum(events_old[1:i])/sum(mu0_old[1:i]),2) }else{result[i+1,6]<- 1} 
result[i+1,7]<- round(LLR(sum(events_old[1:i]),sum(mu0_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,8]<- round(target_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,9]<- round(actual_alpha_old[i],4)}else{result[i+1,9]<- "NA"}
result[i+1,10]<- CVs_old[i]
result[i+1,11]<- paste("No")

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
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", mu0= ",mu0," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)                                                         

x<- seq(start,test,1) ; observed<- as.numeric(result[2:(test+1),5]) ; critical_values<- as.numeric(result[(start+1):(test+1),10])
target<- result[(start+1):(test+1),8] ; actual<- result[(start+1):(test+1),9]
loglike<- as.numeric(result[2:(test+1),7]); RRest<- as.numeric(result[2:(test+1),6])

#########>>>>>>>>>>>>>> critical_values in the scale of MaxSPRT

MaxSPRT_critical_values<- rep(0,length(critical_values))
if(test==1){mu0h<- mu0}else{mu0h<- c(as.numeric(mu0_old),mu0)}
for(i in 1:length(critical_values)){MaxSPRT_critical_values[i]<- LLR(critical_values[i],sum(mu0h[1:i]))}


# Graphic 1
par(mfrow=c(2,2))
plot(seq(1,test,1),rep(max(observed,critical_values),test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Cumulative Events"),main="Number of events scale",sub=title,ylim=c(0,5+max(observed,critical_values)))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),observed,col="blue",pch=20)
lines(seq(1,test,1),observed,col="blue",lty=1)                                              
points(x,critical_values,col="red",pch=20)
lines(x,critical_values,col="red",lty=2)
legend("topleft",c("Needed to reject H0 (CV)","Observed"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

# Graphic 2
plot(seq(1,test,1),rep(alpha+0.02,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Alpha spending"),main="Alpha Spending",ylim=c(0,alpha+0.02))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(x,actual,col="blue",pch=20)
lines(x,actual,col="blue",lty=1)                                              
points(x,target,col="red",pch=20)
lines(x,target,col="red",lty=2)
legend("topleft",c("Target","Actual spent"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

# Graphic 3
plot(seq(1,test,1),rep(max(RRest)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Observed relative risk"),main="Observed Relative Risk",ylim=c(0,max(RRest)+1))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),RRest,col="blue",pch=20)
lines(seq(1,test,1),RRest,col="blue",lty=1)                                              

# Graphic 4
plot(seq(1,test,1),rep(max(loglike)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Observed log-likelihood"),main="MaxSPRT scale",ylim=c(0,max(MaxSPRT_critical_values)+5))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),loglike,col="blue",pch=20)
lines(seq(1,test,1),loglike,col="blue",lty=1)                                              

points(seq(start,test,1),MaxSPRT_critical_values,col="red",pch=20)
lines(seq(start,test,1),MaxSPRT_critical_values,col="red",lty=2)

legend("topleft",c("Needed to reject H0 (CV)","Observed"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")


                                                                  } # CLOSE


###############
### Situation 4: SampleSize achieved without remaining alpha spending

if(start>0&mu0+sum(mu0_old)>=SampleSize&alpha-actualspent<=0.00000001&start>0&reject==0&reject_new==0){ # OPEN

result<- data.frame(matrix(0,test+1,11))

colnames(result)<- c(" "," "," ", "----","Cumulative----"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","mu0","Events","mu0","Events","RR estimate","LLR","target","actual","CV","Reject H0")

result[test+1,1]<- test
result[test+1,2]<- round(mu0,2)
result[test+1,3]<- events
result[test+1,4]<- round(mu0+sum(mu0_old),2)
result[test+1,5]<- events+sum(events_old)
if(events+sum(events_old)>=mu0+sum(mu0_old)){result[test+1,6]<- round((events+sum(events_old))/(mu0+sum(mu0_old)),2) }else{result[test+1,6]<- 1}
result[test+1,7]<- round(LLR(events+sum(events_old),mu0+sum(mu0_old)),2)
result[test+1,8]<- round(current_alpha,4)
result[test+1,9]<- round(actualspent,4)
result[test+1,10]<- CV
result[test+1,11]<- paste("No")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(mu0_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(mu0_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i])
if(sum(events_old[1:i])>=sum(mu0_old[1:i])){result[i+1,6]<- round(sum(events_old[1:i])/sum(mu0_old[1:i]),2) }else{result[i+1,6]<- 1} 
result[i+1,7]<- round(LLR(sum(events_old[1:i]),sum(mu0_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,8]<- round(target_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,9]<- round(actual_alpha_old[i],4)}else{result[i+1,9]<- "NA"}
result[i+1,10]<- CVs_old[i]
result[i+1,11]<- paste("No")

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
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", mu0= ",mu0," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)                                                         

x<- seq(start,test,1) ; observed<- as.numeric(result[2:(test+1),5]) ; critical_values<- as.numeric(result[(start+1):(test+1),10])
target<- result[(start+1):(test+1),8] ; actual<- result[(start+1):(test+1),9]
loglike<- as.numeric(result[2:(test+1),7]); RRest<- as.numeric(result[2:(test+1),6])

#########>>>>>>>>>>>>>> critical_values in the scale of MaxSPRT

MaxSPRT_critical_values<- rep(0,length(critical_values))
if(test==1){mu0h<- mu0}else{mu0h<- c(as.numeric(mu0_old),mu0)}
for(i in 1:length(critical_values)){MaxSPRT_critical_values[i]<- LLR(critical_values[i],sum(mu0h[1:i]))}


# Graphic 1
par(mfrow=c(2,2))
plot(seq(1,test,1),rep(max(observed,critical_values),test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Cumulative Events"),main="Number of events scale",sub=title,ylim=c(0,5+max(observed,critical_values)))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),observed,col="blue",pch=20)
lines(seq(1,test,1),observed,col="blue",lty=1)                                              
points(x,critical_values,col="red",pch=20)
lines(x,critical_values,col="red",lty=2)
legend("topleft",c("Needed to reject H0 (CV)","Observed"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

# Graphic 2
plot(seq(1,test,1),rep(alpha+0.02,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Alpha spending"),main="Alpha Spending",ylim=c(0,alpha+0.02))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(x,actual,col="blue",pch=20)
lines(x,actual,col="blue",lty=1)                                              
points(x,target,col="red",pch=20)
lines(x,target,col="red",lty=2)
legend("topleft",c("Target","Actual spent"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

# Graphic 3
plot(seq(1,test,1),rep(max(RRest)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Observed relative risk"),main="Observed Relative Risk",ylim=c(0,max(RRest)+1))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),RRest,col="blue",pch=20)
lines(seq(1,test,1),RRest,col="blue",lty=1)                                              

# Graphic 4
plot(seq(1,test,1),rep(max(loglike)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Observed log-likelihood"),main="MaxSPRT scale",ylim=c(0,max(MaxSPRT_critical_values)+5))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),loglike,col="blue",pch=20)
lines(seq(1,test,1),loglike,col="blue",lty=1)                                              

points(seq(start,test,1),MaxSPRT_critical_values,col="red",pch=20)
lines(seq(start,test,1),MaxSPRT_critical_values,col="red",lty=2)

legend("topleft",c("Needed to reject H0 (CV)","Observed"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

                                                                                                   } # CLOSE


###############
### Situation 5: H0 rejected in the current test

if(reject==0&reject_new>0){# OPEN

result<- data.frame(matrix(0,test+1,11))

colnames(result)<- c(" "," "," ", "----","Cumulative----"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","mu0","Events","mu0","Events","RR estimate","LLR","target","actual","CV","Reject H0")

result[test+1,1]<- test
result[test+1,2]<- round(mu0,2)
result[test+1,3]<- events
result[test+1,4]<- round(mu0+sum(mu0_old),2)
result[test+1,5]<- events+sum(events_old)
if(events+sum(events_old)>=mu0+sum(mu0_old)){result[test+1,6]<- round((events+sum(events_old))/(mu0+sum(mu0_old)),2) }else{result[test+1,6]<- 1}
result[test+1,7]<- round(LLR(events+sum(events_old),mu0+sum(mu0_old)),2)
result[test+1,8]<- round(current_alpha,4)
result[test+1,9]<- round(actualspent,4)
result[test+1,10]<- CV
result[test+1,11]<- paste("Yes")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(mu0_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(mu0_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i])
if(sum(events_old[1:i])>=sum(mu0_old[1:i])){result[i+1,6]<- round(sum(events_old[1:i])/sum(mu0_old[1:i]),2) }else{result[i+1,6]<- 1} 
result[i+1,7]<- round(LLR(sum(events_old[1:i]),sum(mu0_old[1:i])),2)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,8]<- round(target_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,9]<- round(actual_alpha_old[i],4)}else{result[i+1,9]<- "NA"}
result[i+1,10]<- CVs_old[i]
result[i+1,11]<- paste("No")

                    }
          }


message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE) 
                                              message("=>    Reject H0. No further sequential analyses are needed.",domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", mu0= ",mu0," and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

x<- seq(start,test,1) ; observed<- as.numeric(result[2:(test+1),5]) ; critical_values<- as.numeric(result[(start+1):(test+1),10])
target<- result[(start+1):(test+1),8] ; actual<- result[(start+1):(test+1),9]
loglike<- as.numeric(result[2:(test+1),7]); RRest<- as.numeric(result[2:(test+1),6])

#########>>>>>>>>>>>>>> critical_values in the scale of MaxSPRT

MaxSPRT_critical_values<- rep(0,length(critical_values))
if(test==1){mu0h<- mu0}else{mu0h<- c(as.numeric(mu0_old),mu0)}
for(i in 1:length(critical_values)){MaxSPRT_critical_values[i]<- LLR(critical_values[i],sum(mu0h[1:i]))}


# Graphic 1
par(mfrow=c(2,2))
plot(seq(1,test,1),rep(max(observed,critical_values),test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Cumulative Events"),main="Number of events scale",sub=title,ylim=c(0,5+max(observed,critical_values)))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),observed,col="blue",pch=20)
lines(seq(1,test,1),observed,col="blue",lty=1)                                              
points(x,critical_values,col="red",pch=20)
lines(x,critical_values,col="red",lty=2)
legend("topleft",c("Needed to reject H0 (CV)","Observed"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

# Graphic 2
plot(seq(1,test,1),rep(alpha+0.02,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Alpha spending"),main="Alpha Spending",ylim=c(0,alpha+0.02))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(x,actual,col="blue",pch=20)
lines(x,actual,col="blue",lty=1)                                              
points(x,target,col="red",pch=20)
lines(x,target,col="red",lty=2)
legend("topleft",c("Target","Actual spent"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")

# Graphic 3
plot(seq(1,test,1),rep(max(RRest)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Observed relative risk"),main="Observed Relative Risk",ylim=c(0,max(RRest)+1))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),RRest,col="blue",pch=20)
lines(seq(1,test,1),RRest,col="blue",lty=1)                                              

# Graphic 4
plot(seq(1,test,1),rep(max(loglike)+1,test),col="white",pch=18,xlab="Test",cex.main=1.3,cex.lab=1.3,cex.main=1.5,xaxt="n",
ylab=c("Observed log-likelihood"),main="MaxSPRT scale",ylim=c(0,max(MaxSPRT_critical_values)+5))

sequencia<- seq(1,test,1)
rotulos<- seq(1,test,1)
axis(1, at=sequencia, labels=rotulos,las=1,cex.axis=1.2)

points(seq(1,test,1),loglike,col="blue",pch=20)
lines(seq(1,test,1),loglike,col="blue",lty=1)                                              

points(seq(start,test,1),MaxSPRT_critical_values,col="red",pch=20)
lines(seq(start,test,1),MaxSPRT_critical_values,col="red",lty=2)

legend("topleft",c("Needed to reject H0 (CV)","Observed"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")


            }# CLOSE



###############
### Situation 6: H0 rejected in previous tests

if(reject>0){# OPEN

result<- data.frame(matrix(0,test+1,11))

colnames(result)<- c(" "," "," ", "----","Cumulative----"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","mu0","Events","mu0","Events","RR estimate","LLR","target","actual","CV","Reject H0")


result[test+1,1]<- test
result[test+1,2]<- round(mu0,2)
result[test+1,3]<- events
result[test+1,4]<- round(mu0+sum(mu0_old),2)
result[test+1,5]<- events+sum(events_old)
if(events+sum(events_old)>=mu0+sum(mu0_old)){result[test+1,6]<- round((events+sum(events_old))/(mu0+sum(mu0_old)),2) }else{result[test+1,6]<- 1}
result[test+1,7]<- round(LLR(events+sum(events_old),mu0+sum(mu0_old)),2)
result[test+1,c(8,9,10)]<- paste("NA")
result[test+1,11]<- paste("Yes")

if(test>1){
for(i in 1:(test-1)){

result[i+1,1]<- i
result[i+1,2]<- round(mu0_old[i],2)
result[i+1,3]<- events_old[i]
result[i+1,4]<- round(sum(mu0_old[1:i]),2) 
result[i+1,5]<- sum(events_old[1:i])
if(sum(events_old[1:i])>=sum(mu0_old[1:i])){result[i+1,6]<- round(sum(events_old[1:i])/sum(mu0_old[1:i]),2) }else{result[i+1,6]<- 1} 
result[i+1,7]<- round(LLR(sum(events_old[1:i]),sum(mu0_old[1:i])),2)
if(i>reject){result[i+1,c(8,9,10)]<- paste("NA")}else{
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,8]<- round(target_alpha_old[i],4)}else{result[i+1,8]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,9]<- round(actual_alpha_old[i],4)}else{result[i+1,9]<- "NA"}
result[i+1,10]<- CVs_old[i]
                                                     }
if(i<reject){result[i+1,11]<- paste("No")}else{result[i+1,11]<- paste("Yes")}

                    }
          }


message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)      
          message(paste(c("=>    H0 was rejected on test"," ",reject,". ","No further sequential analyses are needed.")),domain = NULL, appendLF = TRUE)  
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: Sample size= ",SampleSize,", alpha= ",alpha,", mu0= ",mu0," and M= ",M,"."),domain = NULL, appendLF = TRUE)
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
inputSetUp[4,test]<- events
inputSetUp[5,test]<- actualspent
inputSetUp[6,test]<- mu0
if(reject==0){inputSetUp[7,test]<- current_alpha}else{inputSetUp[7,test]<- alpha}

############################################################
## SAVING INFORMATION FOR FUTURE TESTES
############################################################

if(rho==0){write.table(alphaspend,paste(name1,"alphaspend.txt",sep=""))}

write.table(inputSetUp,name)

if(start>0&reject==0&actualspent>0){write.table(p,paste(name1,"p.txt",sep=""))}

result2<- result[2:(test+1),]
colnames(result2)<- c("Test","mu0","Events","Cum. mu0","Cum. Events","RR estimate","LLR","target alpha","actual alpha","CV","Reject H0")
invisible(result2)

#####################################
}##### Close function Analyze.Poisson
#####################################


# AnalyzeSetUp.Poisson(name="VaccineA", SampleSize=100, alpha=0.05,M=1,AlphaSpendType="power-type",rho=0.5,title="n",address="C:/Users/Ivair/Documents")
# AnalyzeSetUp.Poisson(name="CvH MMR",SampleSize=25,alpha=0.05,M=3,AlphaSpendType="Wald",rho="n",title="CvH MMR_MMRV",address="C:/Users/Visitante/Ivair")
# Analyze.Poisson(name="VaccineA",test=1,mu0=1,events=1,AlphaSpend="n")


