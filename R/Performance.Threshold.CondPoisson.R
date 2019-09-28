
Performance.Threshold.CondPoisson<- function(K,cc,CV.upper="n", Person_timeRatioH0="n",GroupSizes="n",Tailed="upper",RR)
{

Threshold<- CV.upper

# Person_timeRatioH0 = critical value in the scale tau_k= P_k/V notation of Silva et al.(2019) instead of the "Threshold" CMaxSPRT scale 

if(Tailed!="upper"){stop("For this version of the Sequential package, AlphaSpend.CondPoisson works only for 'Tailed=upper'.",call. =FALSE)}

if(length(GroupSizes)==1){if(GroupSizes=="n"){GroupSizes<- rep(1,K)};if(is.numeric(GroupSizes)==T){GroupSizes<- rep(GroupSizes,K/GroupSizes)}}

RR<- RR[order(RR)]
AlphaSpend<- rep(0,length(GroupSizes))

if(min(Threshold)==max(Threshold)){Threshold<- rep(Threshold,length(GroupSizes))}


## More verifications
Groups<- GroupSizes
if(length(Groups)==1){
if(is.numeric(Groups)==FALSE){stop("'GroupSizes' must be a single integer or a vector of positive integers summing up 'K'.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'GroupSizes' must be a positive integer or a vector of positive integers summing up 'K'.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'GroupSizes' must be a positive integer or a vector of positive integers summing up 'K'.",call. =FALSE)}

if(Groups==0){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}
if(K/Groups!=round(K/Groups)){stop("The maximum length of surveillance, 'K', must be a multiple of 'GroupsSizes'.",call. =FALSE)}
if(Groups>K){stop("The maximum length of surveillance, 'K', must be a multiple of 'GroupSizes'.",call.=FALSE)}
}

if(length(Groups)>1){
if(sum(is.numeric(Groups))==0){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}else{
if(is.numeric(K)==FALSE){stop("The maximum length of surveillance, 'K', must be a positive integer.",call. =FALSE)}
if(K!=round(K)){stop("The maximum length of surveillance, 'K', must be a positive integer.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}
if(sum(Groups)!=K){stop("'GroupSizes' must sums up 'K'.",call. =FALSE)}
}
}



### IMPORTANT VARIABLE
if(length(Groups)==1){if(Groups=="n"|Groups==1){GroupSizes<- rep(1,K)}else{GroupSizes<- rep(Groups,K/Groups)}}
ks<- GroupSizes

if((is.numeric(K)==FALSE)){stop("'K' must be a positive integer.",call. =FALSE)}
if(K<=0|round(K)!=K){stop("'K' must be a positive integer.",call. =FALSE)}





## AUXILIAR FUNCTIONS
cond.ppois<- function(x,tt){if(x<0){return(0)};y<- seq(0,x); return(sum(exp(y*log(tt)+lfactorial(cc+y-1)-lfactorial(y)-lfactorial(cc-1)-(cc+y)*log(tt+1))))}
cond.dpois<- function(x,tt){if(x<0){return(0)};return(exp(x*log(tt)+lfactorial(cc+x-1)-lfactorial(x)-lfactorial(cc-1)-(cc+x)*log(tt+1)))}
cond.dpois_mult<- function(x,tt){if(x<0){return(0)}else{return(exp(x*log(tt)-lfactorial(x) ))}}
######################     


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



####### Function to find the thresholds in the 'tau' scale for a given cv 
# ------------------------------------------------------------

tau0_old<- rep(0,length(GroupSizes))

for(ii in 1:length(GroupSizes)){
if(is.numeric(Person_timeRatioH0)==FALSE){
ksa<- sum(ks[1:ii])
tau0_old[ii]<- cv_tal(ksa,cc,Threshold[ii])
                                         }else{tau0_old[ii]<- Person_timeRatioH0[ii]}
                               }

tals<- tau0_old
tauR_old<- matrix(tau0_old*RR[1],nrow=1)
if(length(RR)>1){
for(ii in 2:length(RR)){
tauR_old<- rbind(tauR_old,matrix(tau0_old*RR[ii],nrow=1))
                       } 
                }
tauR<- tauR_old




#----- Function that calculates critical values

critical_value<- function(pold,poldR)
{

alpha<- 0
ksa<- sum(ks[1:test])

   for(s in 0:(sum(ks[1:(test-1)])-1)){ alpha<- alpha+ (1-cond.ppois(ksa-1-s,tau0_old[test]-tau0_old[test-1]))*pold[s+1]}

## Updating pold for future tests, here denoted by pf

pf<- rep(0,ksa) ; pfR<- matrix(rep(0,ksa*length(RR)),length(RR),ksa)
for(ki in 0:(ksa-1)){
for(s in 0:min((sum(ks[1:(test-1)])-1),ki)){
if(ki>0){pf[ki+1]<- pf[ki+1]+(cond.ppois(ki-s,tals[test]-tals[test-1])-cond.ppois(ki-s-1,tals[test]-tals[test-1]))*pold[s+1]
for(jj in 1:length(RR)){pfR[jj,ki+1]<- pfR[jj,ki+1]+(cond.ppois(ki-s,tauR_old[jj,test]-tauR_old[jj,test-1])-cond.ppois(ki-s-1,tauR_old[jj,test]-tauR_old[jj,test- 1]))*poldR[jj,s+1]}
                                            }else{
                                                   pf[ki+1]<- cond.ppois(ki-s,tals[test]-tals[test-1])*pold[s+1] 
                                                   
                                                for(jj in 1:length(RR)){pfR[jj,ki+1]<- cond.ppois(ki-s,tauR_old[jj,test]-tauR_old[jj,test-1])*poldR[jj,s+1]}
                                                 }
                                           } # pf[s+1]: probability of having s events at time tau0
                        }


power<- rep(0,length(RR))
for(jj in 1:length(RR)){
 for(s in 0:(sum(ks[1:(test-1)])-1)){ power[jj]<- power[jj]+ (1-cond.ppois(ksa-1-s,tauR_old[jj,test]-tauR_old[jj,test-1]))*poldR[jj,s+1]}
                       }

return(list(pf,power,pfR,alpha))
}



### RUNING THE CALCULATION OF ALPHA SPENDING AND PERFORMANCE MEASURES FOR EACH GROUP


Power<- matrix(rep(0,length(GroupSizes)*length(RR)), length(RR),length(GroupSizes))

for(i in 1:length(GroupSizes)){

test<- i

if(i==1){# open 1

kn<- sum(GroupSizes[1:i])


AlphaSpend[i]<- 1-cond.ppois(kn-1,tau0_old[i])

for(jj in 1:length(RR)){
Power[jj,i]<- 1-cond.ppois(kn-1,tauR_old[jj,1])
                       }
tau0<- tau0_old[i]
## Updating pold for future tests, here denoted by p
p<- rep(0,kn); poldR<- matrix(rep(0,kn*length(RR)),length(RR),kn)
for(s in 0:(kn-1)){if(s>0){p[s+1]<- cond.ppois(s,tau0)-cond.ppois(s-1,tau0)}else{p[s+1]<- cond.ppois(s,tau0)}} # p[s+1]: probability of having s cases attime mu1 
for(jj in 1:length(RR)){
for(s in 0:(kn-1)){if(s>0){poldR[jj,s+1]<- cond.ppois(s,tau0*RR[jj])-cond.ppois(s-1,tau0*RR[jj])}else{poldR[jj,s+1]<- cond.ppois(s,tau0*RR[jj])}}
                       }

        }# close 1



if(i>1){# open 2

res<- critical_value(pold,poldR)
p<- res[[1]];  poldR<- res[[3]]
poweaux<- as.numeric(res[[2]])
AlphaSpend[i]<- res[[4]]
for(jj in 1:length(RR)){ Power[jj,i]<- poweaux[jj]}



       }# close 2

pold<- p

                           }## CLOSE FOR i 


##### PERFORMANCE
Times<- as.numeric(GroupSizes%*%upper.tri(matrix(0,length(GroupSizes),length(GroupSizes)),diag=T))
perf<- matrix(,length(RR),4)
colnames(perf)<- c("RR","power","ESignalTime","ESampleSize")
for(jj in 1:length(RR)){
power<- sum(Power[jj,])
perf[jj,1]<- RR[jj]
perf[jj,2]<- power
perf[jj,3]<- sum(Power[jj,]*Times)/power
perf[jj,4]<- sum(Power[jj,]*Times)+(1-power)*sum(GroupSizes)
                      }

AlphaSpend<- AlphaSpend%*%upper.tri(matrix(0,length(GroupSizes),length(GroupSizes)),diag=T)

results<- list(AlphaSpend,perf)
names(results)<- c("AlphaSpend","Performance")

return(results)
##########################################
}# CLOSE Performance.Threshold.CondPoisson

# Performance.Threshold.CondPoisson(K=30,cc=10,CV.upper=3.6, Person_timeRatioH0="n",GroupSizes=10,RR=2)

# Performance.Threshold.CondPoisson(K=30,cc=10,Threshold="n", Person_timeRatioH0=c(0.1,0.5,1),GroupSizes=10,RR=2)


