


Performance.AlphaSpend.CondPoisson<- function(K,cc,alpha,AlphaSpend="n",GroupSizes="n",rho=0.5,gamma="n",Tailed="upper",RR)
{

if(Tailed!="upper"){stop("For this version of the Sequential package, Threshold.CondPoisson works only for 'Tailed=upper'.",call. =FALSE)}

if(length(GroupSizes)==1){if(GroupSizes=="n"){GroupSizes<- rep(1,K)}}

RR<- RR[order(RR)]
gama<- gamma
############################
### Four classical alpha spending shapes

# AlphaSpend=1: power-type t^rho, Kim and DeMets (1987a), and Jennison & Turnball(1989,1990), noted that this function produces Pocock and O'Brien & Fleming tests.
# Values of 1,2 and 3 for mimicking Wang & Tsiatis(1987) test.

alpha_spendT1<- function(alpha,rho){
x<- as.numeric(GroupSizes%*%((upper.tri(matrix(0,length(GroupSizes),length(GroupSizes)),diag=T)*1))/sum(GroupSizes))
sum_sa<- alpha*(x^rho)
return(sum_sa)
}


# AlphaSpend= 2: Gaussian-type, Lan & DeMets(1983) suggested this one for mimicking O'Brien & Fleming's test.  

alpha_spendT2<- function(alpha){
x<- as.numeric(GroupSizes%*%((upper.tri(matrix(0,length(GroupSizes),length(GroupSizes)),diag=T)*1))/sum(GroupSizes))
za<- qnorm(1-alpha/2)
sum_sa<- 2-2*pnorm(za/sqrt(x))
return(sum_sa)
}


# AlphaSpend= 3: LogExp-type, Lan & DeMets(1983) indicated this one for mimicking Pocock's test.

alpha_spendT3<- function(alpha){
x<- as.numeric(GroupSizes%*%((upper.tri(matrix(0,length(GroupSizes),length(GroupSizes)),diag=T)*1))/sum(GroupSizes))
sum_sa<- alpha*log(1+(exp(1)-1)*x)
return(sum_sa)
}


# AlphaSpend= 4: Gamma-type, Hwang, Shih & DeCani(1990). 

alpha_spendT4<- function(alpha,gamma){
x<- as.numeric(GroupSizes%*%((upper.tri(matrix(0,length(GroupSizes),length(GroupSizes)),diag=T)*1))/sum(GroupSizes))
if(gamma==0){sum_sa<- alpha*x}else{
sum_sa<- alpha*(1-exp(-gamma*x))/(1-exp(-gamma))
                                  }
return(sum_sa)
}


## More verifications
Groups<- GroupSizes
if(length(Groups)==1){
if(is.numeric(Groups)==FALSE){stop("'GroupSizes' must be an integer smaller than or equal to 'K'.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'GroupSizes' must be a positive integer or a vector of positive integers.",call. =FALSE)}

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
if(sum(Groups)!=K){stop("'GroupSizes' must sum up equal to 'K'.",call. =FALSE)}
}
}



### IMPORTANT VARIABLE
if(length(Groups)==1){if(Groups=="n"|Groups==1){GroupSizes<- rep(1,K)}else{GroupSizes<- rep(Groups,K/Groups)}}
ks<- GroupSizes


if((is.numeric(alpha)==FALSE)){stop("'alpha' must be a number in the '(0,0.5]' interval.",call. =FALSE)}
if((alpha<=0|alpha>0.5)){stop("'alpha' must be a number in the '(0,0.5]' interval.",call. =FALSE)}

if((is.numeric(K)==FALSE)){stop("'K' must be a positive integer.",call. =FALSE)}
if(K<=0|round(K)!=K){stop("'K' must be a positive integer.",call. =FALSE)}

AlphaSpend2<- AlphaSpend

if(length(AlphaSpend2)==1){ 
                          if( AlphaSpend2=="n" ){stop(" 'AlphaSpend' must be a vector containing a non-decreasing and positive sequence of numbers with maximum at alpha, with same length as 'GroupSizes', or a single integer among 1 to 4.",call. =FALSE)}
                          if(sum(AlphaSpend2)!=alpha&sum(AlphaSpend2==c(1,2,3,4))==0){stop(" 'AlphaSpend' must be a vector containing a non-decreasing and positive sequence of numbers with maximum at alpha, with same length as 'GroupSizes', or a single integer among 1 to 4.",call. =FALSE)}
                          if(AlphaSpend2==alpha){if(length(GroupSizes)>1){stop(" 'AlphaSpend' and 'GroupSizes' must have the same length.",call. =FALSE)}}
                          if( AlphaSpend2==1 ){
                                              if(length(rho)>1){stop(" 'rho' must be a single number greater than zero and smaller than 5.",call. =FALSE)}               
                                              if(rho=="n"){stop(" Please, choose 'rho' between 0 and 5.",call. =FALSE)} 
                                              if(sum(is.numeric(rho))!=1){stop("Symbols and texts are not applicable for 'rho'. It must be a number greater than zero and smaller than 5.",call. =FALSE)}
                                              if(rho<=0|rho>5){stop(" 'rho' must be a single number greater than zero and smaller than 5.",call. =FALSE)}
                                              AlphaSpend<- alpha_spendT1(alpha,rho)
                                             }
                          if( AlphaSpend2==2 ){AlphaSpend<- alpha_spendT2(alpha)}
                          if( AlphaSpend2==3 ){AlphaSpend<- alpha_spendT3(alpha)}
                          if( AlphaSpend2==4 ){
                                              if(length(gama)>1){stop(" 'gama' must be a single number greater than -10 and smaller than 10.",call. =FALSE)}               
                                              if(gama=="n"){stop(" Please, choose 'gama' between -10 and 10.",call. =FALSE)} 
                                              if(sum(is.numeric(gama))!=1){stop("Symbols and texts are not applicable for 'gama'. It must be a number greater than -10 and smaller than 10.",call. =FALSE)}
                                              if(gama<=0|gama>5){stop(" 'gama' must be a single number greater than -10 and smaller than 10.",call. =FALSE)}
                                              AlphaSpend<- alpha_spendT4(alpha,gama)
                                             }
                         }else{
                               if(is.numeric(AlphaSpend2)==FALSE){stop(" 'AlphaSpend' must be a vector containing a non-decreasing and positive sequence of numbers with maximum at alpha or a single integer among 1 to 4.",call. =FALSE)}
                               if(min(AlphaSpend2)<0|max(AlphaSpend)!=alpha){stop(" 'AlphaSpend' must be a vector containing a non-decreasing and positive sequence of numbers with maximum at alpha or a single integer among 1 to 4.",call. =FALSE)}
                              }








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






#----- Function that calculates critical values

critical_value<- function(pold,current_alpha,poldR)
{

alphas<- current_alpha-actualspent

mu1<- tau0_old[test-1]
mu2<- sum(ks)
tau0<- (mu1+mu2)/2
perror<- 0
scape<- 0
ksa<- sum(ks[1:test])

while(abs(perror-alphas)>0.00000001&scape==0){
   for(s in 0:(sum(ks[1:(test-1)])-1)){ perror<- perror+ (1-cond.ppois(ksa-1-s,tau0-tau0_old[test-1]))*pold[s+1]}
    if(perror>alphas){mu2<- tau0}else{mu1<- tau0}; tau0<- (mu1+mu2)/2
    if(abs(perror-alphas)>0.00000001){perror<- 0}else{scape<- 1} 
                                     }

## Updating pold for future tests, here denoted by pf
tau0_old<- c(tau0_old,tau0)
tauR<- matrix(tau0*RR,,1) ; tauR_old<- cbind(tauR_old,tauR)
tals<- tau0_old



pf<- rep(0,ksa) ; pfR<- matrix(rep(0,ksa*length(RR)),length(RR),ksa)
for(ki in 0:(ksa-1)){
for(s in 0:min((sum(ks[1:(test-1)])-1),ki)){
if(ki>0){pf[ki+1]<- pf[ki+1]+(cond.ppois(ki-s,tals[test]-tals[test-1])-cond.ppois(ki-s-1,tals[test]-tals[test-1]))*pold[s+1]
for(jj in 1:length(RR)){pfR[jj,ki+1]<- pfR[jj,ki+1]+(cond.ppois(ki-s,tauR_old[jj,test]-tauR_old[jj,test-1])-cond.ppois(ki-s-1,tauR_old[jj,test]-tauR_old[jj,test-1]))*poldR[jj,s+1]}
                                            }else{
                                                   pf[ki+1]<- cond.ppois(ki-s,tals[test]-tals[test-1])*pold[s+1] 
                                                   
                                                for(jj in 1:length(RR)){pfR[jj,ki+1]<- cond.ppois(ki-s,tauR_old[jj,test]-tauR_old[jj,test-1])*poldR[jj,s+1]}
                                                 }
                                           } # pf[s+1]: probability of having s events at time tau0
                        }

CVf<- cLLR(ksa,cc,tau0)

power<- rep(0,length(RR))
for(jj in 1:length(RR)){
 for(s in 0:(sum(ks[1:(test-1)])-1)){ power[jj]<- power[jj]+ (1-cond.ppois(ksa-1-s,tauR[jj,1]-tauR_old[jj,test-1]))*poldR[jj,s+1]}
                       }

return(list(CVf,perror,pf,tau0_old,power,pfR,tauR_old))
}



### RUNING THE CALCULATION OF CRITICAL VALUE AND PERFORMANCE MEASURES FOR EACH GROUP

CV<- rep(0,length(GroupSizes))
Power<- matrix(rep(0,length(GroupSizes)*length(RR)), length(RR),length(GroupSizes))

for(i in 1:length(GroupSizes)){

test<- i
current_alpha<- AlphaSpend[i]

if(i==1){# open 1

alphas<- current_alpha

cv1<- 0
cv2<- 10
cvm<- (cv1+cv2)/2
kn<- sum(GroupSizes[1:i])
tau0<- cv_tal(kn,cc,cvm)
perror<- 0
# tau value at the very first chunk of data
while(abs(perror-alphas)>0.00000001){
    perror<- 1-cond.ppois(kn-1,tau0)
    if(perror>alphas){cv1<- cvm}else{cv2<- cvm}; cvm<- (cv1+cv2)/2; tau0<- cv_tal(kn,cc,cvm)
                                     }
tauR<- matrix(tau0*RR,,1)

for(jj in 1:length(RR)){
Power[jj,i]<- 1-cond.ppois(kn-1,tauR[jj,1])
                       }
actualspent<- perror
CV[i]<- cvm

## Updating pold for future tests, here denoted by p
p<- rep(0,kn); poldR<- matrix(rep(0,kn*length(RR)),length(RR),kn)
for(s in 0:(kn-1)){if(s>0){p[s+1]<- cond.ppois(s,tau0)-cond.ppois(s-1,tau0)}else{p[s+1]<- cond.ppois(s,tau0)}} # p[s+1]: probability of having s cases attime mu1 
for(jj in 1:length(RR)){
for(s in 0:(kn-1)){if(s>0){poldR[jj,s+1]<- cond.ppois(s,tau0*RR[jj])-cond.ppois(s-1,tau0*RR[jj])}else{poldR[jj,s+1]<- cond.ppois(s,tau0*RR[jj])}}
                       }

tau0_old<- tau0
tauR_old<- tauR
        }# close 1



if(i>1){# open 2

res<- critical_value(pold,current_alpha,poldR)
CV[i]<- res[[1]]
actualspent<- res[[2]]+actualspent 
p<- res[[3]];  poldR<- res[[6]]
tau0_old<- as.numeric(res[[4]]); tauR_old<- res[[7]]

poweaux<- as.numeric(res[[5]])
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

results<- list(CV,tau0_old,perf)
names(results)<- c("CV","Person-timeRatioH0","Performance")

return(results)
##########################################
}# CLOSE Performance.AlphaSpend.CondPoisson


#res<- Performance.AlphaSpend.CondPoisson(K=30,cc=10,alpha=0.05,AlphaSpend=c(0.01,0.02,0.05),GroupSizes=c(10,10,10),rho=0.5,RR=c(1.2,1.5,2))


