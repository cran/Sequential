
ConfidenceInterval.Binomial<- function(Gamma=0.9, CV.lower="n", CV.upper="n",GroupSizes,z="n",p="n",Cum.cases,Tailed="upper")
{

# Gamma: Confidence coefficient for the interval estimator of the relative risk.
# CV.lower: lower signaling threshold given in the scale of the cases for all tests so far performed until termination of the analysis. Put "NA" for initial looks at the data where a test is not applicable.
# CV.upper: lower signaling threshold given in the scale of the cases for all tests so far performed until termination of the analysis. Put "NA" for initial looks at the data where a test is not applicable.
# GroupSizes: test specific number of events (cases+controls) for all tests so far performed until termination of the analysis.
# z: vector of matching ratios for all tests so far performed until termination of the analysis.
# p: probability of having heads from the Bernoulli experiment.
# Cum.cases: total number of cases on termination of the analysis.

if(length(z)==1){if(z=="n"){z<- 1/p-1}}

if(length(z)==1){z<- rep(z,length(GroupSizes))}

Cum.controls<- sum(GroupSizes)-Cum.cases
GroupSizes<- as.numeric(GroupSizes%*%upper.tri(matrix(0,length(GroupSizes),length(GroupSizes)),diag=TRUE))

if(length(CV.lower)==1){if(CV.lower=="n"){CV.lower<- rep(NA,length(GroupSizes))}}
if(length(CV.upper)==1){if(CV.upper=="n"){CV.upper<- rep(NA,length(GroupSizes))}}
for(ii in 1: length(CV.lower)){if(is.na(CV.lower[ii])==TRUE){CV.lower[ii]<- -1}}
for(ii in 1: length(CV.upper)){if(is.na(CV.upper[ii])==TRUE){CV.upper[ii]<- GroupSizes[ii]+1}}



##############################################################
### Finding the confidence interval by the upper signaling threshold:

if(Tailed=="upper"|Tailed=="two"){ # open 1

#### Auxiliar function for calculating the probability of crossing the upper limit throshold
####
Dist_Cases<- function(RRh){ # opens Dist_Cases

if(length(GroupSizes)==1){
  prob_reject<- 1-pbinom(Cum.cases-1,GroupSizes[1],1/(1+z[1]/RRh))
  return(prob_reject)
                         }



if(length(GroupSizes)>1){ # opens GroupSizes

for(i in 1:(length(GroupSizes)-1)){ # opens i
   probs<- rep(0,CV.upper[i]-1-(CV.lower[i]+1)+1)
   k<-0
for(j in (CV.lower[i]+1):(CV.upper[i]-1)){ # opens j
k<- k+1
if(i==1){probs[k]<- dbinom(j,GroupSizes[i],1/(1+(z[i]/RRh)))}else{
kk<- 0;   for(l in (CV.lower[i-1]+1):(CV.upper[i-1]-1)){
kk<- kk+1
    probs[k]<- probs[k]+dbinom(j-l,GroupSizes[i]-GroupSizes[i-1],1/(1+(z[i]/RRh)))*probs_o[kk] 
                                                     }
                                                            }
                                          } # closes j

if(i==1){prob_reject<- 1-pbinom(CV.upper[i]-1,GroupSizes[i],1/(1+(z[i]/RRh)))}
if(i>1){
kk<- 0
  for(l in (CV.lower[i-1]+1):(CV.upper[i-1]-1)){
kk<- kk+1
prob_reject<- prob_reject+ (1-pbinom(CV.upper[i]-l-1,GroupSizes[i]-GroupSizes[i-1],1/(1+(z[i]/RRh))))*probs_o[kk]
                                              }
       }

probs_o<- probs             
                               } # closes i

kk<- 0
for(k in (CV.lower[i]+1):(CV.upper[i]-1)){kk<- kk+1; prob_reject<- prob_reject+(1-pbinom(Cum.cases-k-1,GroupSizes[i+1]-GroupSizes[i],1/(1+(z[length(z)]/RRh))))*probs[kk]}
return(prob_reject)
                       }  # closes GroupSizes


                        } # closes Dist_Cases 
####
####


##### Finding the lower limit:
RR1<- 0; RR2<- 1000; RRl<- (RR1+RR2)/2

pp<- 0
while(abs(pp-(1-Gamma)/2)>0.001){
pp<- Dist_Cases(RRl)
if(pp>(1-Gamma)/2){RR2<- RRl}else{RR1<- RRl}
RRl<- (RR1+RR2)/2
                                }

##### Finding the upper limit:
RR1<- RRl; RR2<- 1000; RRu<- (RR1+RR2)/2
pp<- 0
while(abs(pp-(1+Gamma)/2)>0.001){
pp<- Dist_Cases(RRu)
if(pp>(1+Gamma)/2){RR2<- RRu}else{RR1<- RRu}
RRu<- (RR1+RR2)/2
                                }


            } # close 1
###################################################################

##############################################################
### Finding the confidence interval by the lower signaling threshold:

if(Tailed=="lower"){ # open 2

Dist_Cases<- function(RRh)
{


if(length(GroupSizes)==1){
  prob_reject<- pbinom(Cum.cases,GroupSizes[1],1/(1+z[1]/RRh))
  return(prob_reject)
                         }

if(length(GroupSizes)>1){ # opens GroupSizes

for(i in 1:(length(GroupSizes)-1)){
   probs<- rep(0,CV.upper[i]-1-(CV.lower[i]+1)+1)
   k<-0
for(j in (CV.lower[i]+1):(CV.upper[i]-1)){
k<- k+1
if(i==1){probs[k]<- dbinom(j,GroupSizes[i],1/(1+(z[i]/RRh)))}else{
kk<- 0;   for(l in (CV.lower[i-1]+1):(CV.upper[i-1]-1)){
kk<- kk+1
    probs[k]<- probs[k]+dbinom(j-l,GroupSizes[i]-GroupSizes[i-1],1/(1+(z[i]/RRh)))*probs_o[kk] 
                                                     }
                                                            }
                                       } # close j

if(i==1){prob_reject<- pbinom(CV.lower[i],GroupSizes[i],1/(1+(z[i]/RRh)))}
if(i>1){
kk<- 0
  for(l in (CV.lower[i-1]+1):(CV.upper[i-1]-1)){
kk<- kk+1
prob_reject<- prob_reject+ pbinom(CV.lower[i]-l,GroupSizes[i]-GroupSizes[i-1],1/(1+(z[i]/RRh)))*probs_o[kk]
                                              }
       }

probs_o<- probs             
                      } # close i

kk<- 0
for(k in (CV.lower[i]+1):(CV.upper[i]-1)){kk<- kk+1; prob_reject<- prob_reject+pbinom(Cum.cases-k,GroupSizes[i+1]-GroupSizes[i],1/(1+(z[length(z)]/RRh)))*probs[kk]}
return(prob_reject)
                     } # closes GroupSizes

                        } # closes Dist_Cases 
####
####

##### Finding the lower limit:
RR1<- 0; RR2<- 1000; RRl<- (RR1+RR2)/2

pp<- 0
while(abs(pp-(1-Gamma)/2)>0.001){
pp<- Dist_Cases(RRl)
if(pp>(1-Gamma)/2){RR1<- RRl}else{RR2<- RRl}
RRl<- (RR1+RR2)/2
                                }

##### Finding the upper limit:
RR1<- RRl; RR2<- 1000; RRu<- (RR1+RR2)/2

while(abs(pp-(1+Gamma)/2)>0.001){
pp<- Dist_Cases(RRu)
if(pp>(1+Gamma)/2){RR1<- RRu}else{RR2<- RRu}
RRu<- (RR1+RR2)/2
                                }


            } # close 2
###################################################################

res<- list(RRl,RRu)
names(res)<- c("RRl","RRu")
return(res)
}

