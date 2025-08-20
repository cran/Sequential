# -------------------------------------------------------------------------
# Function to perform the unpredictable multinomial marginal MaxSPRT surveillance - Version 4
# -------------------------------------------------------------------------

Analyze.Multinomial<- function(name,test,cases,controls,N_exposures,N_controls,AlphaSpend="n")
{

# name: name to be used in each analysis to read the information saved from previus test.
# test: the number of tests already performed plus the current test
# cases: vector with the number of events per multinomial entry
# controls: the number of events in the control window.
# N_exposures: vector with the number of individuals in the same risk window per multinomial entry excluding the entry for the control group. Must have the same dimension as cases. 
# N_controls: number of individuals in the control group 
# AlphaSpend: Acummulative alpha spending up to the current test. The default is the power-type with rho defined in the AnalyzeSetUp.Multinomial function.



#### Setting the local where the files from previous tests are located

name1<- name

safedir<- getwd()
address1<- tempdir()
y<- substr(address1,1:4,4) ; i<- 2
while(i<=nchar(address1)-3&y!="Temp"&y!="TEMP"){y<- substr(address1,i:(i+3),i+3);i<- i+1}
address1<- substr(address1,1:(i+3),i+3)
setwd(address1)

if(file.exists(paste(name,"address.txt",sep=""))==FALSE){
stop(c("There is no file called"," '",name,"' yet."," You must run the function 'AnalyzeSetUp.Binomial' first."),call. =FALSE)
                                                        }

address<- as.character(read.table(paste(name,"address.txt",sep=""),sep=";")[1,1])
setwd(address)

name<- paste(name,".","txt",sep="")

if(file.exists(name)==FALSE){
stop(c("There is no file called"," '",name,"' yet."," You must run the function 'AnalyzeSetUp.Binomial' first."),call. =FALSE)
                            }


titlecheck<- paste(name1,"title.txt",sep="")
title<- read.table(titlecheck)
title<- title[1,1]
if(title==0){title<- " "}else{title<- as.character(title)}





NN<- sum(cases)+sum(controls) # number of events in the current test.





####
## Uploading information from previous tests
####

inputSetUp<- read.table(name)
if(inputSetUp[1,1]!=test-1){stop(c("The current test should be"," ",inputSetUp[1,1]+1,". ",
"If you do not have information about previous tests, see the user manual for more details."),call. =FALSE)}

k<- as.numeric(inputSetUp[1,4]) # number of exposures
ExposuresNames<- read.table("ExposuresNames.txt") # TESTAR AQUI
ExposuresNames<- ExposuresNames[1:k,1]

## inputSetUp matrix contains:
# line 1: (C11) the index for the order of the test (zero entry if we did not have a first test yet), (C12) Maximum SampleSize, (C13) alpha, (C14) k, (C15) m, (C16) title, (C17) reject (the index indicating if and when H0 was rejected), (C18) rho, (C19) vazio, (C1,10) vazio
# line 2: says if the analysis has been started or not. 0 for not started. 1 for already started.
# line 3: critical values in the scale of MaxSPRT
# line 4: sum of the cases per test  
# line 5: AlphaOld, actual alpha spent up to the last previous test
# line 6: cumulative expected value of Sn E[sum(w*C)] under H0, test by test
# line 7: has the target alpha spending until the (test-1)th look.
# line 8: observed controls
# line 9: rejection time. It has zero for those exposures not rejected.
# line 10: the number of rows (different weights/z ) per test
# line 11: matched case-control
# line 12: Target alpha spending defined with AnalyzeSetUp
# line 13: weight for each observation
# line 14: test per weight
# line 15: vector with test margins R0.
# line 16: AlphaSpendType
# line 17: active, vector informing which entries are still elegible for alpha spending in the current test. Entries must be settled equal to 1 for active(elegible), or equal to zero if the null hypothesis related to it has been rejected in some previous test.
# line 18: Ns, vector with test-specific number of events
# line 19: N_old, total cumulative number of events in previous tests



AlphaSpendType<- inputSetUp[16,1]
pmin<- inputSetUp[16,2]
pmax<- inputSetUp[16,3]
target_power<- inputSetUp[16,4] 
R1<- inputSetUp[16,5] # relative risk under the alternative hypothesis given the target_power
N<- as.numeric(inputSetUp[1,2])
alpha<- as.numeric(inputSetUp[1,3])
m<- as.numeric(inputSetUp[1,5])
rho<- as.numeric(inputSetUp[1,8])
Rmin<- as.numeric(inputSetUp[16,6])
Rmax<- as.numeric(inputSetUp[16,7])
gamma<- as.numeric(inputSetUp[16,8])

R0<- as.numeric(inputSetUp[15,1])
active<- as.numeric(inputSetUp[17,1:k])
if(test>1){Ns<- as.numeric(inputSetUp[18,1:(test-1)])}; if(test==1){Ns<- NN}else{Ns<- c(Ns,NN)}
if(test>1){N_old<- as.numeric(inputSetUp[19,test-1])}else{N_old<- 0}
if(test>1){cvs_old<- inputSetUp[20:(20+k),1:(test-1)]; cv_llr_old<- inputSetUp[3,test-1]}
if(test>2){CumCases<- inputSetUp[(20+k+1):(20+2*k),1:(test-1)]}; if(test==2){CumCases<- matrix(inputSetUp[(20+k+1):(20+2*k),1:(test-1)],k,1)}
if(test>1){CumControls<- as.numeric(inputSetUp[20+3*k+1,1:(test-1)])}
if(test==1){AlphaOld<- 0}else{AlphaOld<- as.numeric(inputSetUp[5,test-1])} 
if(test>2){ps<- inputSetUp[(20+3*k+2):(20+3*k+2+k),1:(test-1)]}; if(test==2){ps<- matrix(inputSetUp[(20+3*k+2):(20+3*k+2+k),1:(test-1)],k+1,1)}



# NN: number of events in the current test.
# alpha: overall significance level
# N: maximum length of surveillance
# probs_old: probabilities of each state under the non-reject H0 status in the previous tests.
# cvs_old: critical values used in the previous tests.
# N_old: total cumulative number of events in previous tests
# Ns: vector with test-specific number of events 
# R0: test margin under the null hypothesis
# m: Monte Carlo replications of the multinomial for critical values calculations
# active: vector informing which entries are still elegible for alpha spending in the current test. Entries must be settled equal to 1 for active(elegible), or equal to zero if the null hypothesis related to it has been rejected in some previous test. 
# ExposuresNames: Optional. This is to inform the name of the exposures related to each entry of the multinomial vector. For example, it can be c("A","B","AB") indicating that the count entries in the "cases" vector are related to populations exposed to vaccines A, B, and AB, respectively.
# rho: parameter for the power-type alpha spending function



#### vector with probabilities under H0 per multinomial entry excluding the entry for the control group. Must have the same dimension as cases.
p<- N_exposures/(sum(N_exposures)+N_controls) 

p_h0<- c(p,1-sum(p))  # probabilities of events under H0

#### Updating the multinomial probability vector for the active entries
activeh<- c(active,1) 
p0a<- p_h0
p0a[activeh==0]<- 0 
for(i in 1:(k+1)){p_h0[i]<- p0a[i]/sum(p0a)}

if(test==1){ps<- matrix(p_h0,k+1,1)}else{ps<- cbind(ps,matrix(p_h0,k+1,1))}









#### Setting the robust alpha spending function

if(AlphaSpendType==1){alpha1<- alpha/2; alpha2<- alpha/2}else{alpha1<- alpha/(k+1); alpha2<- k*alpha/(k+1)}



NNs<- matrix(0,k,test)

for(i in 1:test){if(i!=test){NNs[,i]<- Ns[i]}else{NNs[,i]<- NN}}


if(AlphaSpend=="n"){

Mm<- c(as.numeric(inputSetUp[9,1:k]),1)
Rejection<- matrix(0,k+1,test); Rejection[k+1,]<- 1

for(i in 1:k){if(Mm[i]>0){Rejection[i,Mm[i]:test]<- 1}}





powers<- matrix(0,k+1,test); powers[k+1,]<- 1
 
if(test>1){
for(i in 1:k){
   for(j in 1:test){
pss<- as.numeric(ps[i,1:j])
Nss<- Ns[1:j][pss>0]
      if(as.numeric(ps[i,j])==0|sum(Nss)<=2){if(j>1){powers[i,j]<- powers[i,j-1]}}else{


pss<- pss[pss>0]

     if(AlphaOld>0){res<- Performance.AlphaSpend.Binomial(N=sum(Nss), alpha= AlphaOld/k,
AlphaSpend=1,p=pss,GroupSizes=Nss,Tailed="upper",RR=R1,
Statistic="MaxSPRT",rho)
     powers[i,j]<- res$Performance[2]}else{powers[i,j]<- 0}
                                                                          }
                   }
             }
           }




I<- 1*( (ps>=pmin & ps<=pmax) | Rejection==1 | powers>= target_power)

if(test>1){MI<- (I[1:k,]*NNs[1:k,])}else{MI<- matrix((I[1:k,]*NNs[1:k,]),,1)}



contribution<- 0 
for(l in 1:k){contribution<- contribution + min( alpha2/k, (alpha2/k)*(sum(MI[l,])/N)^(rho) ) }

AlphaSpend = min( alpha1,alpha1*((N_old+NN)/N)^rho ) + contribution

                       
                   }




if(AlphaSpend<=AlphaOld){
stop(c("Choose AlphaSpend greater than"," ",AlphaOld),call. =FALSE)
                        }



alphah<- AlphaSpend - AlphaOld










#### Actual lower bound for power

power<- matrix(0,k,test)
  
for(i in 1:k){
   for(j in 1:test){
pss<- as.numeric(ps[i,1:j])
Nss<- Ns[1:j][pss>0]
      if(as.numeric(ps[i,j])==0|sum(Nss)<=2){
                                                if(j>1){power[i,j]<- power[i,j-1]}
                                                if(j==1){power[i,j]<- 0}
                                                                         }else{

pss<- pss[pss>0]

     if(AlphaSpend>0){res<- Performance.AlphaSpend.Binomial(N=sum(Nss), alpha= AlphaSpend/k,
AlphaSpend=1,p=pss,GroupSizes=Nss,Tailed="upper",RR=R1,
Statistic="MaxSPRT",rho)
     power[i,j]<- res$Performance[2]}else{power[i,j]<- 0}
                                                                              }
                   }
             }





#### Maximum likelihood estimator for the vector of relative risks


if(test>1){xxx<- CumCases[,ncol(CumCases)]+cases ; n<- CumControls[length(CumControls)] + controls + sum(xxx) }else{xxx<- cases; n<- controls+sum(xxx)}
hp<- xxx/n  # MLE for the unknown p


valid_entries<- seq(1,k)[N_exposures>0]

R_partial<- matrix(0,length(valid_entries),length(valid_entries))

for(i in 1:length(valid_entries)){R_partial[i,]<-
 -N_exposures[valid_entries]*hp[valid_entries[i]]; R_partial[i,i]<- N_exposures[valid_entries[i]]*(1-hp[valid_entries[i]])}

hR<-  solve(R_partial)%*%matrix(N_controls*hp[valid_entries],length(valid_entries),1)  # MLE for the relative risk vector, R.


hR1<- rep(0,k); hR1[valid_entries]<- hR; hR<- matrix(hR1,k,1) 





#### Confidence Interval for RR


CI_RR<- matrix(0,k,2)

for(i in 1:k){

# finding the lower bound of the confidence interval for Ri
 
RR<- rep(Rmin,k)
RR1<- 0.001; RR2<- 100; RRm<- (RR1+RR2)/2; RR[i]<- RRm
prob<- 0
while(abs(prob-(1-gamma)/2)>10^(-6) & RRm>0.001){

pi<- N_exposures[i]*RRm/(sum(RR*N_exposures)+N_controls)

prob<- 1-pbinom(xxx[i]-1,n,pi)

if(prob>(1-gamma)/2){RR2<- RRm}else{RR1<- RRm}; RRm<- (RR1+RR2)/2; RR[i]<- RRm 

                                     }

CI_RR[i,1]<- RRm


# finding the upper bound of the confidence interval for Ri
 
RR<- rep(Rmax,k)
RR1<- 0.001; RR2<- 100; RRm<- (RR1+RR2)/2; RR[i]<- RRm
prob<- 0
while(abs(prob-(1-gamma)/2)>10^(-6) &RRm<100){

pi<- N_exposures[i]*RRm/(sum(RR*N_exposures)+N_controls)

prob<- pbinom(xxx[i],n,pi)

if(prob>(1-gamma)/2){RR1<- RRm}else{RR2<- RRm}; RRm<- (RR1+RR2)/2; RR[i]<- RRm 

                                     }

CI_RR[i,2]<- RRm

            }
 





############################################################
###### INTERNAL AUXILIARY FUNCTIONS
############################################################


###########################################
#----- THE MAXSPRT STATISTIC

LLR <- function(cc,n,z){

       if(cc==n){x = n*log(1+z/R0)}else{
         if((z/R0)*cc/(n-cc)<=1){x=0}else{
	       x = cc*log(cc/n)+(n-cc)*log((n-cc)/n)  -cc*log(1/(z/R0+1))-(n-cc)*log((z/R0)/(z/R0+1))
                                    }
                                  } 	
      	x
	                 }
#--------------------------







#-------------------------------------------------
## Finding the critical values for the current test

#### SAMPLE SPACE IN THE LLR SCALE
aux3<- 0
for(i in 1:k){
if(p[i]>0&p[i]<1){
  for(x in 0:(NN+N_old)){
             z<- 1/p[i]-1
             if(aux3==0){aux3<- 1; llrs<-  LLR(x,NN+N_old,z)}else{
             llrh<- LLR(x,NN+N_old,z)
             llrs<-  c(llrs, llrh ) 
                                                                }  
                       }
                 }
             }

llrs<- unique(llrs); llrs<- llrs[order(llrs)]

#------





#library(pmultinom)

if(test==1){

cv_new<- rep(NN+1,k+1 )

    ii<- length(llrs)
    probE1<- 0
       while(probE1<= alphah& ii>0){ 
           cv_llr<- llrs[ii] 
           for(j in 1:k){
                         if(p[j]>0 & p[j]<1){
                         z<- 1/p[j]-1
                         x<- NN; llr2<- LLR(x,NN,z) ; while(llr2>=cv_llr){cv_new[j]<- x; x<- x-1; llr2<- LLR(x,NN,z) }
                                            }
                        }      
          probE1<- max(0,1-pmultinom(upper = cv_new-1, size=NN, probs= p_h0, method="exact"))
          if(probE1<= alphah){ ii<- ii-1; probE1ref<- probE1; cvs<- cv_new}   
                                  }
  if(ii==0|probE1>alphah){cv_new<- rep(N_old+NN+1,k+1); probE1ref<- 0}else{cv_new<- cvs}
           } 





### VERSION 4



if(test>1){


cv_new<- rep(N_old+NN+1,k+1 )

    ii<- sum(llrs<cv_llr_old/2)
    
    probE1<- 1
       while(probE1> alphah & ii<=length(llrs)){ # first while
           cv_llr<- llrs[ii] 
           for(j in 1:k){
                        if(p[j]>0 & p[j]<1){
                        z<- 1/p[j]-1
 ### MUDEI O while(llr2>=cv_llr) ABAIXO POR while(llr2>cv_llr)
                        x<- NN+N_old; llr2<- LLR(x,NN+N_old,z) ; while(llr2>cv_llr){cv_new[j]<- x; x<- x-1; llr2<- LLR(x,NN+N_old,z) }
                                           }
                        } 
             
cvs<- cbind(cvs_old,cv_new)

ys<- matrix(0,k+1,test)
G<- 0
 for(i in 1:m){
         ys[,1]<- rmultinom(1,as.numeric(Ns[1]),p_h0); aux<- sum(sum(ys[,1]<cvs[,1])==(k+1))
         for(j in 2:test){ys[,j]<- ys[,j-1]+rmultinom(1,as.numeric(Ns[j]),p_h0); if(j==test){aux<- aux+sum(sum(ys[1:k,j]>=cvs[1:k,j])>0)}else{aux<- aux+sum(sum(ys[,j]<cvs[,j])==(k+1))}} 
         G<- G + sum(aux==test)        
              }
         
         
     
          probE1<- G/m
          if(probE1> alphah){ ii<- ii+1}else{ probE1ref<- probE1}   
                                  }# close first while
  cv_new<- cvs ; if(ii== (length(llrs)+1) ){cv_new<- rep(N_old+NN+1,k+1)}
 

          }







############################################################
## UPDATING INFORMATION FOR FUTURE TESTES
############################################################

if(ncol(inputSetUp)<test){inputSetUp<- cbind(inputSetUp,rep(0,nrow(inputSetUp)))}
N_old<- N_old+NN
AlphaOld<- AlphaOld + probE1ref
cvs_old<- cv_new
cv_llr_old<- cv_llr


inputSetUp[1,1]<- test
inputSetUp[3,test]<- cv_llr_old
inputSetUp[4,test]<- sum(cases)
inputSetUp[5,test]<- AlphaOld
inputSetUp[7,test]<- AlphaSpend
inputSetUp[8,test]<- controls
inputSetUp[18,test]<- NN
inputSetUp[19,test]<- N_old
if(test==1){inputSetUp[20:(20+k),test]<- cvs_old}else{
inputSetUp[20:(20+k),test]<- cvs_old[,test] # critical values for this test
                                                     }
# CumCases
if(test==1){inputSetUp[(20+k+1):(20+2*k),test]<- cases}else{inputSetUp[(20+k+1):(20+2*k),test]<- cases + as.numeric(inputSetUp[(20+k+1):(20+2*k),test-1])}   

if(test==1){new_active<- 1*( cv_new[1:k]> as.numeric(inputSetUp[(20+k+1):(20+2*k),test]) )}else{
new_active<- active*1*( cv_new[1:k,test]> as.numeric(inputSetUp[(20+k+1):(20+2*k),test]) )
                                                                                    }






inputSetUp[17,1:k]<- new_active

inputSetUp[(20+2*k+1):(20+3*k),test]<- hR  # MLE for the relative risk vector 

# CumControls
if(test>1){inputSetUp[20+3*k+1,test]<- inputSetUp[20+3*k+1,test-1]+controls}else{inputSetUp[20+3*k+1,test]<- controls} 

# ps

inputSetUp[(20+3*k+2):(20+3*k+2+k),test]<- p_h0


if(sum(new_active==0&active==1)>0){inputSetUp[9,new_active==0&active==1]<- test}  
Reject_Test_Time<- inputSetUp[9,1:k]
if(sum(Reject_Test_Time==0)>0){Reject_Test_Time[Reject_Test_Time==0]<- "na"}
names(Reject_Test_Time)<- ExposuresNames
rownames(Reject_Test_Time)<- " "


# Confidence intervals for the relative risk

inputSetUp[(20+4*k+3):(20+4*k+3+k-1),test]<- CI_RR[,1]
inputSetUp[(20+4*k+3+k):(20+5*k+3+k-1),test]<- CI_RR[,2]



############################################################
## SAVING INFORMATION FOR FUTURE TESTS
############################################################

write.table(inputSetUp,name)


##########################################################
## PRINTING OUTPUT TABLES
##########################################################

linhas<- rep(0,test)
for(i in 1:test){linhas[i]<- paste("test",i)}

Reject_H0<- rep("No",k); Reject_H0[new_active==0]<- "Yes"
names(Reject_H0)<- ExposuresNames



ps_under_H0<- round(t(ps[1:k,]),4)
colnames(ps_under_H0)<- ExposuresNames
rownames(ps_under_H0)<- linhas


power<- round(t(power),4)
colnames(power)<- ExposuresNames
rownames(power)<- linhas



Critical_Values<- t(inputSetUp[20:(20+k-1),1:test])
colnames(Critical_Values)<- ExposuresNames
rownames(Critical_Values)<- linhas

Cumulative_Cases<-  t(inputSetUp[(20+k+1):(20+2*k),1:test])
Cumulative_Cases<- cbind(Cumulative_Cases, matrix(c(inputSetUp[20+3*k+1,1:test]),,1))
colnames(Cumulative_Cases)<- c(ExposuresNames,"Controls")
rownames(Cumulative_Cases)<- linhas

Relative_Risk_estimates<-  round(t(inputSetUp[(20+2*k+1):(20+3*k),1:test]),2)
colnames(Relative_Risk_estimates)<- ExposuresNames
rownames(Relative_Risk_estimates)<- linhas

Critical_Values_LLR<- matrix(round(as.numeric(inputSetUp[3,1:test]),6),1,test)
colnames(Critical_Values_LLR)<- linhas
rownames(Critical_Values_LLR)<- " "

Alpha_spending<- matrix(0,2,test); rownames(Alpha_spending)<- c("Target","Actual")
Alpha_spending[1,]<- round(as.numeric(inputSetUp[7,1:test]),6)
Alpha_spending[2,]<- round(as.numeric(inputSetUp[5,1:test]),6) 
colnames(Alpha_spending)<- linhas

Lower_bound_CI<- round(t(inputSetUp[(20+4*k+3):(20+4*k+3+k-1),1:test]),2)
colnames(Lower_bound_CI)<- ExposuresNames
rownames(Lower_bound_CI)<- linhas


Upper_bound_CI<- round(t(inputSetUp[(20+4*k+3+k):(20+5*k+3+k-1),1:test]),2)
colnames(Upper_bound_CI)<- ExposuresNames
rownames(Upper_bound_CI)<- linhas


result<- list(Reject_H0,Reject_Test_Time,ps_under_H0,power,Critical_Values,Cumulative_Cases,Relative_Risk_estimates,Critical_Values_LLR,Alpha_spending,Lower_bound_CI,Upper_bound_CI)
names(result)<- c("Reject_H0","Rejection_time","ps_under_H0","power","Critical_values_in_cumulative_cases_scale","Cumulative_cases","Relative_risk_estimates","Critical_values_in_MaxSPRT_scale","Alpha_spending","Lower_bound_CI","Upper_bound_CI")

invisible(result)

#####################################
}##### Close function Analyze.Multinomial
#####################################




## Illustrative example 

#library(Sequential)

#AnalyzeSetUp.Multinomial(name="testee",N=200,alpha=0.05,k=7,R0=1,rho=1,m=10000,
#title="tituto da tabela",ExposuresNames=c("A","B","C","AB","AC","AD","ABC"),
#address="C:/Users/User/Documents/Viagens a Boston/2024/BACKUP DEVIDO AO PROBLEMA DE VIRUS/TRABALHO V2/MULTIPLE VACCINES/TESTES")


# Test 1
#res1<- Analyze.Multinomial(name="testee",test=1,cases=c(8,6,1,4,1,0,0),controls= 3,N_exposures=c(1000,1500,300,300,320,280,240),N_controls=500,AlphaSpend="n")
#res1

# Test 2
#res2<- Analyze.Multinomial(name="testee",test=2,cases=c(13,1,4,6,4,2,1),controls= 5,N_exposures=c(900,1200,800,280,250,250,275),N_controls=1500,AlphaSpend="n")
#res2

# Test 3
#res3<- Analyze.Multinomial(name="testee",test=3,cases=c(8,5,3,3,3,1,7),controls= 3,N_exposures=c(800,1000,1000,230,280,310,300),N_controls=2000,AlphaSpend="n")
#res3

# Test 4
#res4<- Analyze.Multinomial(name="testee",test=4,cases=c(6,6,5,3,6,3,8),controls= 8,N_exposures=c(700,900,1000,200,300,350,310),N_controls=3000,AlphaSpend="n")
#res4

# Test 5
#res5<- Analyze.Multinomial(name="testee",test=5,cases=c(4,2,3,4,6,2,5),controls= 7,N_exposures=c(800, 850, 900, 200, 250, 300, 200),N_controls=1000,AlphaSpend="n")
#res5

# Test 6
#res6<- Analyze.Multinomial(name="testee",test=6,cases=c(6,3,3,6,7,3,6),controls= 8,N_exposures=c(815, 790, 880, 210, 240, 310, 208),N_controls=950,AlphaSpend="n")
#res6


#### Graph with the alpha spending of this example:

#target<- as.numeric(res6$Alpha_spending[1,])
#actual<- as.numeric(res6$Alpha_spending[2,])
#ns<- c( sum(as.numeric(res6$Cumulative_cases[1,])), sum(as.numeric(res6$Cumulative_cases[2,])),  sum(as.numeric(res6$Cumulative_cases[3,])) , sum(as.numeric(res6$Cumulative_cases[4,])) , sum(as.numeric(res6$Cumulative_cases[5,])), sum(as.numeric(res6$Cumulative_cases[6,])) )
#plot(ns, target, type="b", lty=1, lwd=2, xlab="Sample size", ylab= "Alpha spending",cex=1.5) 
#lines(ns, actual, type="b", lty=2, lwd= 2)
#abline(h=0.05)
#legend(550,0.02, c("Target", "Actual"),lty= c(1,2),pch=1, lwd=2, bty="n")
#text(locator(7),c("(23, 0.004527)", "(59, 0.013427)", "(92, 0.019827)", "(137, 0.032627)", "(170, 0.037327)", "(212, 0.037627)", "approx 0.036573"),cex=1.5)





