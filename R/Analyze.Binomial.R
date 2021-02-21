
Analyze.Binomial<- function(name,test,z="n",p="n",cases,controls,AlphaSpend="n")
{

CIgamma="n"; dplaces=3

if(length(CIgamma)>1){stop(" 'CIgamma' must be a single number in the '(0.5,1)' interval.",call. =FALSE)}
if(is.numeric(CIgamma)==F&CIgamma!="n"){stop(" 'CIgamma' must be a single number in the '(0.5,1)' interval.",call. =FALSE)}
if(is.numeric(CIgamma)==T){if(CIgamma>=1|CIgamma<=0.5){stop(" 'CIgamma' must be a single number in the '(0.5,1)' interval.",call. =FALSE)}}

if(length(dplaces)>1){stop(" 'dplaces' must be a positive integer.",call. =FALSE)}
if(as.numeric(dplaces)==F){stop(" 'dplaces'  must be a positive integer.",call. =FALSE)}
if(round(dplaces)!=dplaces|dplaces<=0|dplaces>6){stop(" 'dplaces' must be a positive integer in the '[1,6]' interval.",call. =FALSE)}

if(length(p)==1&length(z)==1){if(p=="n"&z=="n"){stop("Please, at least one of the inputs, z or p, must be provided.",call. =FALSE)}}

if(length(z)>1){if(sum(is.numeric(z))!=1){stop("Symbols and texts are not applicable for 'z'. It must be a vector of positive numbers.",call. =FALSE)}}
if(length(z)==1){if(z!="n"&is.numeric(z)!=TRUE){stop("Symbols and texts are not applicable for 'z'. It must be a vector of positive numbers.",call. =FALSE)}}
if(sum(z<=0)>0){stop("'z' must be a vector of positive numbers.",call. =FALSE)}

if(length(p)==1){
if(p!="n"){
if(sum(is.numeric(p))!=1){stop("Symbols and texts are not applicable for 'p'. It must contain only probability measures.",call. =FALSE)}
if(sum(p<=0)>0|sum(p>=1)>0){stop("Each entry of p must be a number greater than zero and smaller than 1.",call. =FALSE)}

if(length(z)==1){
    if(z!="n"){
if(length(z)!=length(p)){stop("'z' and 'p' are vectors that must have the same dimension.",call. =FALSE)}
if(sum(p== 1/(1+z))!=length(p)&z!="n"){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed. .",call. =FALSE)}
              }
                }
if(length(z)>1){
    
if(length(z)!=length(p)){stop("'z' and 'p' are vectors that must have the same dimension.",call. =FALSE)}
if(sum(p== 1/(1+z))!=length(p)&z!="n"){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed. .",call. =FALSE)}
              
                }

          }
                }

if(length(p)>1){
if(sum(is.numeric(p))!=1){stop("Symbols and texts are not applicable for 'p'. It must contain only probability measures.",call. =FALSE)}
if(sum(p<=0)>0|sum(p>=1)>0){stop("Each entry of p must be a number greater than zero and smaller than 1.",call. =FALSE)}
if(length(z)==1){
    if(z!="n"){
if(length(z)!=length(p)){stop("'z' and 'p' are vectors that must have the same dimension.",call. =FALSE)}
if(sum(p== 1/(1+z))!=length(p)&z!="n"){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed. .",call. =FALSE)}
              }
                }
if(length(z)>1){
   
if(length(z)!=length(p)){stop("'z' and 'p' are vectors that must have the same dimension.",call. =FALSE)}
if(sum(p== 1/(1+z))!=length(p)){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed. .",call. =FALSE)}
              
                }
                }


if(length(z)==1){if(z=="n"){z<- 1/p-1}}


ExposureA<- cases
ExposureB<- controls
w<- rep(1,length(cases)) 


# Function to specify the number of decimal places to be printed inside the output tables
specify_decimal<- function(x, k){format(round(x, k), nsmall=k)} 

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



cases<- ExposureA
controls<- ExposureB


if( sum(is.numeric(cases))!=1){stop("Symbols and texts are not applicable for 'cases'. It must be a vector of counts.",call. =FALSE)}
if( sum(is.numeric(controls))!=1){stop("Symbols and texts are not applicable for 'controls'. It must be a vector of counts.",call. =FALSE)}
#if( sum(is.numeric(w))!=1){stop("Symbols and texts are not applicable for 'w'. It must be a vector of positive numbers.",call. =FALSE)}
if( sum(is.numeric(z))!=1){stop("Symbols and texts are not applicable for 'z'. It must be a vector of positive numbers.",call. =FALSE)}


if(sum(cases<0)>0){stop("cases must contain only zeros and positive integers.",call. =FALSE)}
if(sum(controls<0)>0){stop("controls must contain only zeros and positive integers.",call. =FALSE)}
#if(sum(w<0)>0){stop("w must contain only positive numbers.",call. =FALSE)}
if(sum(z<0)>0){stop("z must contain only zeros and positive numbers.",call. =FALSE)}



if( sum(is.numeric(test))!=1|length(test)>1){stop("Symbols and texts are not applicable for 'test'. It must be an integer greater than zero.",call. =FALSE)}

if(test<0){stop("'test' must be an integer greater than zero.",call. =FALSE)}

if(sum(cases)+sum(controls)<=0){stop("It must be informed a number greater than zero for the total number of events.",call. =FALSE)}

if( length(names(table(c(length(cases),length(controls),length(z),length(w)))))>1 ){stop("ExposureA, ExposureB, z, and w must be of the same length.",call. =FALSE)}




####
## Uploading information from previous tests
####

inputSetUp<- read.table(name)
if(inputSetUp[1,1]!=test-1){stop(c("The current test should be"," ",inputSetUp[1,1]+1,". ",
"If you do not have information about previous tests, see the user manual for more details."),call. =FALSE)}



## inputSetUp matrix contains:
# line 1: (C11) the index for the order of the test (zero entry if we did not have a first test yet), (C12) SampleSize, (C13) alpha, (C14) M, (C15) base(the line of p where the looping will start in the next test), (C16) title, (C17) reject (the index indicating if and when H0 was rejected), (C18) rho (zero if Wald is used), (C19) Tailed (1 for one-sided test and 2 for two-sided), (C1,10) z (used for alpha spending=Wald)
# line 2: says if the analysis has been started or not. 0 for not started. 1 for already started.
# line 3: CVu which is the upper critical values in the scale of Sn
# line 4: observed cases  
# line 5: actual alpha spent
# line 6: cumulative expected value of Sn E[sum(w*C)] under H0, test by test
# line 7: has the target alpha spending until the (test-1)th look.
# line 8: observed controls
# line 9: CVl which is the lower critical values in the scale of Sn
# line 10: the number of rows (different weights/z ) per test
# line 11: matched case-control (old z)
# line 12: Target alpha spending defined with AnalyzeSetUp
# line 13: weight for each observation (old w)
# line 14: test associated to each weight 
# lines 15: the first column has R0

#### 

SampleSize<- inputSetUp[1,2] ; N<- SampleSize
alpha<- inputSetUp[1,3]
M<- inputSetUp[1,4]
start<- inputSetUp[2,1]
reject<- inputSetUp[1,7]
rho<- inputSetUp[1,8] ; if(rho==0){rho<- "n"}
Tailed<- inputSetUp[1,9]
maxspent<- max(inputSetUp[5,])
actual_alpha_old<- inputSetUp[5,]
target_alpha_old<- inputSetUp[7,]
row_old<- inputSetUp[10,1:test]; if(test>1){row_old<- row_old[1:(test-1)]} 
CVu_old<- as.numeric(inputSetUp[3,]) 
cases_old<- as.numeric(inputSetUp[4,1:sum(row_old)])  # <= cases_old  
CVl_old<- as.numeric(inputSetUp[9,])  # <= colocar CVl aqui
controls_old<- as.numeric(inputSetUp[8,1:sum(row_old)]) # <= controls_old  
z_old<- as.numeric(inputSetUp[11,1:sum(row_old)])        # z_old   
alphaspend<- inputSetUp[12,]        # 
w_old<- as.numeric(inputSetUp[13,1:sum(row_old)])         # w_old   
test_old<- as.numeric(inputSetUp[14,1:sum(row_old)])      # test indexes
z_setup<- inputSetUp[1,10] ; if(z_setup==0){z_setup<- "n"}
R0<- inputSetUp[15,1]

if(test>1){
cases_current<- as.numeric(c(cases_old,cases))
controls_current<- as.numeric(c(controls_old, controls))
rows_current<- as.numeric(c(row_old,length(cases))) 
z_current<- as.numeric(c(z_old,z))
          }else{
cases_current<- cases
controls_current<- controls
rows_current<- length(cases) 
z_current<- z
               }

if(CIgamma=="n"){CIgamma<- 1-alpha}


#### More checks

if( sum(is.numeric(AlphaSpend))!=1&AlphaSpend!="n"){stop("Symbols and texts are not applicable for 'AlphaSpend'. If you want to use the default, use 'n'. Otherwise,  'AlphaSpend' must be a positive number smaller than or equal to 'alpha'.",call. =FALSE)}

if( sum(length(AlphaSpend))!=1){stop("'AlphaSpend' must be a single value, not a vector.",call. =FALSE)}

if(AlphaSpend<0&AlphaSpend!="n"){stop("'AlphaSpend' must be a positive number smaller than 'alpha'.",call.=FALSE)}

if(AlphaSpend>alpha&AlphaSpend!="n"){stop(c("'AlphaSpend' must be smaller than or equal to ",alpha,"."),call. =FALSE)}


if(sum(cases_old)+sum(controls_old)>=M){
name2<- paste(name1,"pold.txt",sep="")
pold<- matrix(as.numeric(read.table(name2)[[1]]), ,1) # has the previous state probabilities 
name2<- paste(name1,"statesOld.txt",sep="")
statesOld<- as.numeric(read.table(name2)[[1]]) # has the possible state of Sn until (test-1)th analysis.
                                       }

##################################################################
## Setting up the target alpha spending.

kref<- sum(cases+controls)+sum(controls_old+cases_old)

if(SampleSize<=kref){current_alpha<- alpha}else{
if(max(inputSetUp[5,])<alpha-0.00000001&inputSetUp[1,7]==0&kref>=M){

if(AlphaSpend=="n"){
  if(inputSetUp[1,8]==0){ # for Wald type of alpha spending
                         if(N<kref){
                         current_alpha<- alpha
                                   }else{current_alpha<- as.numeric(alphaspend[kref]) }
                        }else{
                         current_alpha<- alpha*(kref/inputSetUp[1,2])^rho
                             }
                   }else{
                         if(AlphaSpend<=max(inputSetUp[5,])&test>1){stop(c("For this test, 'AlphaSpend' must be selected in the (", specify_decimal(max(inputSetUp[5,]),6), ",", alpha,"] interval because it has already been spent up to ",specify_decimal(max(inputSetUp[5,]),6)," until the previous test."),call. =FALSE)}
                         current_alpha<- AlphaSpend
                         if(kref<=N){alphaspend[kref]<- AlphaSpend; for(ii in kref:N){if(alphaspend[ii]<alphaspend[ii-1]){alphaspend[ii]<- alphaspend[ii-1]}}}
                        }

                                                                            }
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

###########################################
## Expectation of Cases under H0. j is the order of the test for which the expectation is to be calculated
Exp<- function(ord){
cum.events<- cases_current[1:sum(rows_current[1:ord])]+controls_current[1:sum(rows_current[1:ord])]
zs<- z_current[1:sum(rows_current[1:ord])]
return(sum(cum.events/(1+zs/R0))) 
} 

###########################################
## Relative risk estimation

RelRisk<- function(ord)
{
lenZs<- length(names(table(z_current[1:sum(rows_current[1:ord])])))
if(lenZs==1){ 
rrest<- min(z)*sum(cases_current[1:sum(rows_current[1:ord])])/sum(controls_current[1:sum(rows_current[1:ord])])
          }else{
counta<- 1
rrest<- 0
for(yy in 1:sum(rows_current[1:ord])){rrest<- max(rrest,max(z_current[1:yy])*sum(cases_current[1:yy])/sum(controls_current[1:yy]))}

LR<- function(Rcand){return(prod( dbinom(cases_current[1:sum(rows_current[1:ord])],cases_current[1:sum(rows_current[1:ord])]+controls_current[1:sum(rows_current[1:ord])],1/(1+(z_current[1:sum(rows_current[1:ord])])/Rcand))) )}
if(rrest<Inf){Rsup<- rrest}else{Rsup<- 200} 
if(max(controls_current[1:sum(rows_current[1:ord])])==0){rrest<- Inf}else{
tesaux<- 0 ; Rinf<- 0
while(tesaux==0){Rs<- matrix(seq(Rinf,Rsup,min(0.01,Rsup)),,1); LRs<- apply(Rs,1,LR); rrest<- Rs[LRs==max(LRs)];if(rrest< counta*200|counta>20){tesaux<- 1}else{counta<- counta+1 ; Rinf<- Rsup; Rsup<- Rsup+200}}
                                                                         }
if(counta>20){rrest<- Inf}
               }

return(rrest)
}


###########################################
## LLR for cumulative data with variable z

LLRcum<- function(ord)
{
Rest<- RelRisk(ord)
LR<- function(Rcand){return(prod( dbinom(cases_current[1:sum(rows_current[1:ord])],cases_current[1:sum(rows_current[1:ord])]+controls_current[1:sum(rows_current[1:ord])],1/(1+(z_current[1:sum(rows_current[1:ord])])/Rcand))) )}
if(max(controls_current[1:sum(rows_current[1:ord])])==0){llrcum<- -log(LR(R0))}
if(max(cases_current[1:sum(rows_current[1:ord])])==0|Rest<=R0){llrcum<- 0}else{
llrcum<- log(LR(Rest))-log(LR(R0))
                                                                             }
return(llrcum)
}


# Function to calculate the marginal distribution of Sn for an unique group of data
##########################################

sn_marginal<- function(z,w,cases,controls,R=R0)
{
z_w_levels<- names(table(paste(z,w)))
data_matrix<- matrix(0,length(z_w_levels),5)
colnames(data_matrix)<- c("z","w","Cases","Controls","Total")

for(ii in 1:length(z_w_levels)){ # 1
sum_cases<- sum(cases[z_w_levels[ii]==paste(z,w)]) 
sum_controls<- sum(controls[z_w_levels[ii]==paste(z,w)])
sum_total<- sum_cases+sum_controls

data_matrix[ii,]<- c(max(z[z_w_levels[ii]==paste(z,w)]),max(w[z_w_levels[ii]==paste(z,w)]),sum_cases,sum_controls,sum_cases+sum_controls)

pii<- 1/(1+max(z[z_w_levels[ii]==paste(z,w)])/R)
 
if(ii==1){ # 2
states=seq(0,sum_total)*max(w[z_w_levels[ii]==paste(z,w)])
probs_b<- dbinom(seq(0,sum_total),sum_total,pii)
         }else{


### List of states for the current ii

states_matrix<- matrix(0,sum_total+1,length(states_b))
colnames(states_matrix)<-  states_b
rownames(states_matrix)<- seq(0,sum_total)

for(jj in 0:sum_total){
 states_matrix[jj+1,]<- states_b+jj*max(w[z_w_levels[ii]==paste(z,w)]) 
                       }

### Obtaining the probabilities for each state

probs<- rep(0,nrow(states_matrix)*ncol(states_matrix))

for(iii in 1:nrow(states_matrix)){for(jjj in 1:ncol(states_matrix)){
probs[(iii-1)*ncol(states_matrix)+jjj]<- probs_b[jjj]*dbinom(iii-1,sum_total,pii)
if(iii==1 & jjj==1){states<- states_matrix[1,1]}else{states<- c(states,states_matrix[iii,jjj])}
                                                                   }
                                 }

probs_b<- as.numeric(by(probs, states, FUN=sum, simplify = TRUE))

         } # close 2

states_b<- as.numeric(by(states,states,max))                         
                            } # close 1

probs<- as.numeric(by(probs_b,states_b,sum))
states<- as.numeric(by(states_b,states_b,max))

return(list(states,probs))

} # close sn_marginal function
##############################




#####
#----- Function that calculates critical values

critical_value<- function(z,w,casesr,controlsr,current_alpha,maxspent)
{

if(sum(cases_old)+sum(controls_old)<M){#1
     if(sum(cases_old)+sum(controls_old)>=M){res<- sn_marginal(z,w,casesr,controlsr)}else{res<- sn_marginal(z,w,casesr+sum(cases_old),controlsr+sum(controls_old))}
     statesNew<- res[[1]]
     pnew<- res[[2]]
                                       }#1 close if test==1



if(sum(cases_old)+sum(controls_old)>=M){#2

if(sum(cases_old)+sum(controls_old)>=M){resm<- sn_marginal(z,w,casesr,controlsr)}else{resm<- sn_marginal(z,w,casesr+sum(cases_old),controlsr+sum(controls_old))}
statesMed<- resm[[1]]
pmed<- resm[[2]]

paux<- rep(0,length(statesOld)*length(statesMed))
statesaux<- paux 
cc<- 1
for(kk in 1:length(statesOld)){
   for(jj in 1:length(statesMed)){
     
                           statesaux[cc]<- statesOld[kk]+statesMed[jj]                                        
                           paux[cc]<- pold[kk]*pmed[jj]
 cc<- cc+1
                                 }
                              }
statesNew<- as.numeric(names(table(statesaux)))
pnew<- aggregate(paux, by=list(statesaux), FUN=sum)[,2]

                                       }#2 close if test>1


perrorI<- current_alpha-maxspent 

if(Tailed==1){ # 3

jj<- length(pnew) ; i<- jj ; p<- 0 ; pp<-0
while(p<perrorI){p<- sum(pnew[length(pnew):jj]); if(p<perrorI){i<- jj;jj<- jj-1;pp<- p}}


alphahere<- pp              ## actual alpha spent in the current test 
if(pp>0){CVf<- statesNew[i]}else{CVf<- NA}
statesOld<- statesNew[1:i]
pold<- pnew[1:i]

return(list(CVf,alphahere,pold,statesOld))
            } # 3

if(Tailed==2){ # 4

jj<- length(pnew) ; iu<- jj ; p<- 0 ; pp1<-0 
while(p< perrorI/2){p<- sum(pnew[length(pnew):jj]); if(p< perrorI/2){iu<- jj;jj<- jj-1;pp1<- p}}

if(pp1>0){CVfu<- statesNew[iu]}else{CVfu<- NA}

jj<- 1 ; il<- 1 ; p<- 0 ; pp2<-0
while(p< perrorI/2){p<- sum(pnew[1:jj]); if(p< perrorI/2){il<- jj;jj<- jj+1;pp2<- p}}

if(pp2>0){CVfl<- statesNew[il]}else{CVfl<- NA}

alphahere<- pp1+pp2              ## actual alpha spent in the current test 
statesOld<- statesNew[il:iu]
pold<- pnew[il:iu]

return(list(CVfu,CVfl,alphahere,pold,statesOld))
            } # 4


}
#--------------------------
#### Closes the critical value function

##########################################################
###### CALCULATING CRITICAL VALUE AND ACTUAL ALPHA SPENT FOR THE CURRENT TEST
##########################################################

# Finding critical value for the 'current_alpha'

totalevents<- kref

if(totalevents >= M&reject==0&max(inputSetUp[5,])<alpha-0.00000001 ){

res<- critical_value(z,w,casesr=cases,controlsr=controls,current_alpha,maxspent)

if(Tailed==1){
CVu<- res[[1]]
actualspent<- res[[2]]+max(actual_alpha_old)
if(actualspent==0){CVu<- NA}
pold<- res[[3]]
statesOld<- res[[4]]
             }

if(Tailed==2){
CVu<- res[[1]]
CVl<- res[[2]]
actualspent<- res[[3]]+max(actual_alpha_old)
if(actualspent==0){CVu<- NA; CVl<- NA}
pold<- res[[4]]
statesOld<- res[[5]]
             }

# Surveillance started?
if(start>0){Sn<- sum(cases*w)+sum(cases_old*w_old)}else{if(test==1){Sn<- sum(cases*w)}else{Sn<- sum(cases*w)+sum(cases_old*w_old)}}

if(start==0){if(actualspent>0){start<- test}}             

# H0 rejected?
if(Tailed==1){
if(is.na(CVu)==TRUE){reject_new<- 0}else{  if(Sn>=CVu){reject_new<- test}else{reject_new<- 0}}
             }else{
aux1<- 0 ; aux2<- 0
if(is.na(CVu)!=TRUE){if(Sn>=CVu){aux2<- 1}}
if(is.na(CVl)!=TRUE){if(Sn<=CVl){aux1<- 1}}
if(aux1==1|aux2==1){reject_new<- test}else{reject_new<- 0}
                  }
                                                                     }else{reject_new<- max(0,reject)}

if(M>totalevents |reject==1| max(inputSetUp[5,])>=alpha-0.00000001 ){actualspent<- 0; CVu<- NA;  if(Tailed==2){CVl<- NA} }



##########################################################
###### PRINTING TABLES AND GRAPHS WITH RESULTS
##########################################################


# First we have to colapse the data per w/z
z_w_levels<- names(table(paste(z,w)))
data_matrix<- matrix(0,length(z_w_levels),5)
colnames(data_matrix)<- c("z","w","Cases","Controls","Total")
for(ii in 1:length(z_w_levels)){ # 1
sum_cases<- sum(cases[z_w_levels[ii]==paste(z,w)]) 
sum_controls<- sum(controls[z_w_levels[ii]==paste(z,w)])
sum_total<- sum_cases+sum_controls
data_matrix[ii,]<- c(max(z[z_w_levels[ii]==paste(z,w)]),max(w[z_w_levels[ii]==paste(z,w)]),sum_cases,sum_controls,sum_cases+sum_controls)
                               }


# Relative risk estimate
rre<- z*cases/controls 
if(test>1){rre_old<- z_old*cases_old/controls_old}

# Weight history
if(test>1){ws<- c(w_old,w)}else{ws<- w}
ws<- as.numeric(by(ws,ws,max))
nw<- length(ws)
namesw<- rep(0,nw)
for(k in 1:nw){if(k>1){namesw[k]<- paste("w=",ws[k])}else{namesw[k]<- paste("  w=",ws[k])}}


###############
### Situation 1: SampleSize not achieved and surveillance not started because events are still smaller than M

if(start==0){ # OPEN

if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,12))
result[2:(test+1),1]<- seq(1,test,1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases|H0]","RR estimate","LLR","target","actual","CV","Reject H0")

result[2:(test+1),1]<- seq(1,test)
result[test+1,2]<- sum(cases)
result[test+1,3]<- sum(controls)
result[test+1,4]<- sum(cases)+sum(cases_old)
result[test+1,5]<- sum(controls)+sum(controls_old)
result[test+1,6]<- specify_decimal(Exp(test),2)
result[test+1,7]<- specify_decimal(RelRisk(test),2)
result[test+1,8]<- specify_decimal(LLRcum(test),6)
result[test+1,9]<- 0
result[test+1,10]<- 0
result[test+1,11]<-  NA
result[test+1,12]<- paste("No")

if(test>1){

for(i in 1:(test-1)){

if(i==1){a1<- 1; a2<- sum(rows_current[1:i])}else{a1<- sum(rows_current[1:(i-1)])+1; a2<- sum(rows_current[1:i])}

result[i+1,2]<- sum(cases_current[a1:a2])
result[i+1,3]<- sum(controls_current[a1:a2])
result[i+1,4]<- sum(cases_current[1:sum(rows_current[1:i])])
result[i+1,5]<- sum(controls_current[1:sum(rows_current[1:i])])
result[i+1,6]<- specify_decimal(Exp(i),2)
result[i+1,7]<- specify_decimal(RelRisk(i),2)
result[i+1,8]<- specify_decimal(LLRcum(i),6)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,9]<- 0}else{result[i+1,9]<- NA}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,10]<- 0}else{result[i+1,10]<- NA}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,11]<- CVu_old[i]}else{result[i+1,11]<- NA}
result[i+1,12]<- paste("No")

                    }

          }

               } # close 1


####
if(Tailed==2){ # 2

result<- data.frame(matrix(0,test+1,2*nw+7))

names<- rep(0,2*nw+7)
names[1]<- c(" ")

if(nw==1){names[2]<- c("#Events"); names[3]<- c("Relative Risk")}else{
if(nw==2){names[2]<- c("--") ; names[3]<- c("#Events --"); names[4]<- " --Relative"; names[5]<- "Risk--"}else{
if(nw==3){names[2]<- c("--");names[3]<- c("---"); names[4]<- c("#Events"); names[5]<- " --"; names[6]<- "Relative"; names[7]<- "Risk---"}else{
ac<- floor(nw/2+1);  names[2:(ac-1)]<- c("--"); names[ac]<- c("--"); names[ac+1]<- c("#Events "); names[(ac+2):(nw+1)]<- c("--")
if(nw==4){
names[(nw+2)]<- " ---"; names[nw+3]<- c("Relative"); names[nw+4]<- c("Risk"); names[nw+5]<- c("---") 
         }else{
names[(nw+2)]<- " ---"; names[(nw+3):(nw+ac)]<- "--" ;names[nw+ac+1]<- c("Relative"); names[nw+ac+2]<- c("Risk"); names[(nw+ac+3):(2*nw+1)]<- c("--"); names[2*nw+1]<- c("---")
              }          
                                                                                                                              }   
                                                                                                        }
                                                                      }
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "--Critical"; names[1+nw+nw+3]<- "Value--" ; names[1+nw+nw+4]<- "--alpha"; names[1+nw+nw+5]<- "-----" ; names[1+nw+nw+6]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "low"
result[1,1+nw+nw+3]<- "high"
result[1,1+nw+nw+4]<- "target"
result[1,1+nw+nw+5]<- "spent"
result[1,1+nw+nw+6]<- "Reject H0"

result[test+1,1]<- test

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)

for(k in 1:nw){
if(sum(ws[k]==w)>0){cum_events[k]<- sum(cases[ws[k]==w])+sum(controls[ws[k]==w]); cum_cases[k]<- sum(cases[ws[k]==w]); cum_controls[k]<- sum(controls[ws[k]==w])}
if(sum(w_old==ws[k])>0){
cum_events[k]<- cum_events[k]+sum(cases_old[w_old==ws[k]])+sum(controls_old[w_old==ws[k]])
cum_cases[k]<- cum_cases[k]+sum(cases_old[w_old==ws[k]])
cum_controls[k]<- cum_controls[k]+sum(controls_old[w_old==ws[k]]) 
                       }
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&is.na(max(rre[w==ws[k]]))!=TRUE&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- NA}
 result[test+1,1+k]<- as.numeric(cum_events[k]) 
              }

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

result[test+1,2*nw+3]<- NA  # <= low critical value

result[test+1,2*nw+4]<- NA  # <= high critical value

result[test+1,2*nw+5]<- 0
result[test+1,2*nw+6]<- 0
result[test+1,2*nw+7]<- paste("No")

if(test>1){
for(i in 1:(test-1)){
result[i+1,1]<- i

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)
for(k in 1:nw){
  if(sum(ws[k]==w_old)>0){
               cum_events[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])+sum(controls_old[ws[k]==w_old&test_old<=i])
               cum_cases[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])
               cum_controls[k]<- sum(controls_old[ws[k]==w_old&test_old<=i])
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&is.na(max(rre_old[w_old==ws[k]&test_old==i]))!=TRUE&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(rre_old[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- NA}
               result[i+1,1+k]<- as.numeric(cum_events[k])
                         }else{result[i+1,nw+1+k]<- NA}
    
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVl_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVl_old[i]/(sum(cum_events*ws)-CVl_old[i]),2)}else{result[i+1,2*nw+3]<- NA}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+4]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- NA}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+5]<- NA}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+6]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+6]<- NA}
result[i+1,2*nw+7]<- paste("No")
                    }
          }


             } # close 2


message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE) 
                                              message("=>    H0 cannot be rejected yet because the cumulative events is still smaller than M.",domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho, ", Tailed= ", Tailed,", and M= ",M, ", H0: RR<=",R0, "."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)


           } # CLOSE




###############
### Situation 2: Surveillance started, H0 not rejected yet and sample size not achieved

if(reject==0&reject_new==0&start>0&totalevents<SampleSize){# OPEN

if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,12))
result[2:(test+1),1]<- seq(1,test,1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases|H0]","RR estimate","LLR","target","actual","CV","Reject H0")

result[2:(test+1),1]<- seq(1,test)
result[test+1,2]<- sum(cases)
result[test+1,3]<- sum(controls)
result[test+1,4]<- sum(cases)+sum(cases_old)
result[test+1,5]<- sum(controls)+sum(controls_old)
result[test+1,6]<- specify_decimal(Exp(test),2)
result[test+1,7]<- specify_decimal(RelRisk(test),2)
result[test+1,8]<- specify_decimal(LLRcum(test),6)
result[test+1,9]<- specify_decimal(current_alpha,4)
result[test+1,10]<- specify_decimal(actualspent,4)
result[test+1,11]<-  CVu
result[test+1,12]<- paste("No")

if(test>1){

for(i in 1:(test-1)){

if(i==1){a1<- 1; a2<- sum(rows_current[1:i])}else{a1<- sum(rows_current[1:(i-1)])+1; a2<- sum(rows_current[1:i])}

result[i+1,2]<- sum(cases_current[a1:a2])
result[i+1,3]<- sum(controls_current[a1:a2])
result[i+1,4]<- sum(cases_current[1:sum(rows_current[1:i])])
result[i+1,5]<- sum(controls_current[1:sum(rows_current[1:i])])
result[i+1,6]<- specify_decimal(Exp(i),2)
result[i+1,7]<- specify_decimal(RelRisk(i),2)
result[i+1,8]<- specify_decimal(LLRcum(i),6)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,9]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,9]<- paste(NA)}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,10]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,10]<- paste(NA)}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,11]<- CVu_old[i]}else{result[i+1,11]<- paste(NA)}
result[i+1,12]<- paste("No")

                    }

          }

###### Graphics

par(mfrow=c(2,2))

## GRAPH FOR CV, OBSERVED, AND EXPECTED CASES

limy<-  3*(sum(cases)+sum(cases_old))/2+10
plot(seq(1,test), rep(limy,test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Cumulative Cases"),main=title,ylim= c(0,limy) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)

points(seq(1,test), result[2:(test+1),4],col="blue",pch=20)
lines(seq(1,test), result[2:(test+1),4],col="blue",lty=1)
points(seq(1,test), result[2:(test+1),6],col="black",pch=2)
lines(seq(1,test), result[2:(test+1),6],col="black",lty=1)

ini<- 1
if(is.na(result[2,11])==TRUE){
while(is.na(result[ini+1,11])==TRUE){ini<- ini+1}; if(ini==test){cvs<- result[ini+1,11]}else{cvs<- rep(0,test-ini+1) ; cvs[1]<- result[ini+1,11]; g1<- 1; for(gg in (ini+1):test){g1<- g1+1; if(is.na(result[gg+1,11])==TRUE){cvs[g1]<-result[gg,11]}else{cvs[g1]<-result[gg+1,11] }}}
cvs<- as.numeric(cvs)
points(seq(ini,test),cvs ,col="red",pch=20)
lines(seq(ini,test), cvs,col="red",lty=1) 
legend("topleft",c("Needed to reject H0 (CV)","Observed","Expected = E[Cases|H0]"),col=c("red","blue","black"),pch=c(18,20,2),lty=c(2,1,1),bty="n")
                                            }else{
legend("topleft",c("Observed","Expected = E[Cases|H0]"),col=c("blue","black"),pch=c(20,2),lty=c(1,1),bty="n")
                                                 }

## GRAPH FOR ALPHA SPEND 

if(sum(is.na(result[2:(test+1),9]))<test){
plot(seq(1,test), rep(as.numeric(result[test+1,9])+0.02,test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Alpha spending"),main="Alpha Spending",ylim= c(0,as.numeric(result[test+1,9])+0.02 ) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(ini,test),result[(ini+1):(test+1),9],lty=2,col="red")
points(seq(ini,test),result[(ini+1):(test+1),9],pch=18,col="red")
lines(seq(ini,test),result[(ini+1):(test+1),10],lty=2,col="blue")
points(seq(ini,test),result[(ini+1):(test+1),10],pch=18,col="blue")
legend("topleft",c("Target","Actual spent"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")
                                           }else{
plot(seq(1,test,1),rep(alpha,test),pch=18,col="white",ylab="Alpha spending",xlab="Test",xaxt="n",cex.lab=1.3, 
sub=" ",font.sub=1,main="Alpha spending"); legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")
                                                } 


## GRAPH FOR RR ESTIMATES

if(sum(result[2:(test+1),7]==Inf)<test){
ini<- 1; while(result[ini+1,7]==Inf){ini<- ini+1}
plot(seq(1,test,1),rep(max(as.numeric(result[(ini+1):(test+1),7])),test),pch=18,col="white",ylab="Observed relative risk",xlab="Test",xaxt="n",cex.lab=1.3, 
ylim=c(0,max(as.numeric(result[(ini+1):(test+1),7]))),sub=" ",font.sub=1,main="Observed Relative Risk")
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(ini,test),result[(ini+1):(test+1),7],lty=2,col="blue")
points(seq(ini,test),result[(ini+1):(test+1),7],pch=18,col="blue")
                                        }else{
plot(seq(1,test,1),rep(1,test),pch=18,col="white",ylab="Observed relative risk",xlab="Test",xaxt="n",cex.lab=1.3, 
sub=" ",font.sub=1,main="Observed relative risk"); legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")
                                             }


## GRAPH FOR LLR

plot(seq(1,test), rep(max(as.numeric(result[2:(test+1),8])),test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Log-likelihood ratio"),main="Log-likelihood ratio",ylim= c(0,max(as.numeric(result[2:(test+1),8])) ) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(1,test),result[2:(test+1),8],lty=2,col="blue")
points(seq(1,test),result[2:(test+1),8],pch=18,col="blue")


             } # close 1


####
if(Tailed==2){ # 2

if(test>1){ws<- c(w_old,w)}else{ws<- w}
ws<- as.numeric(by(ws,ws,max))
nw<- length(ws)
namesw<- rep(0,nw)
for(k in 1:nw){if(k>1){namesw[k]<- paste("w=",ws[k])}else{namesw[k]<- paste("  w=",ws[k])}}

result<- data.frame(matrix(0,test+1,2*nw+7))

names<- rep(0,2*nw+7)
names[1]<- c(" ")

if(nw==1){names[2]<- c("#Events"); names[3]<- c("Relative Risk")}else{
if(nw==2){names[2]<- c("--") ; names[3]<- c("#Events --"); names[4]<- " --Relative"; names[5]<- "Risk--"}else{
if(nw==3){names[2]<- c("--");names[3]<- c("---"); names[4]<- c("#Events"); names[5]<- " --"; names[6]<- "Relative"; names[7]<- "Risk---"}else{
ac<- floor(nw/2+1);  names[2:(ac-1)]<- c("--"); names[ac]<- c("--"); names[ac+1]<- c("#Events "); names[(ac+2):(nw+1)]<- c("--")
if(nw==4){
names[(nw+2)]<- " ---"; names[nw+3]<- c("Relative"); names[nw+4]<- c("Risk"); names[nw+5]<- c("---") 
         }else{
names[(nw+2)]<- " ---"; names[(nw+3):(nw+ac)]<- "--" ;names[nw+ac+1]<- c("Relative"); names[nw+ac+2]<- c("Risk"); names[(nw+ac+3):(2*nw+1)]<- c("--"); names[2*nw+1]<- c("---")
              }          
                                                                                                                              }   
                                                                                                        }
                                                                      }
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "--Critical"; names[1+nw+nw+3]<- "Value--" ; names[1+nw+nw+4]<- "--alpha"; names[1+nw+nw+5]<- "-----" ; names[1+nw+nw+6]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "low"
result[1,1+nw+nw+3]<- "high"
result[1,1+nw+nw+4]<- "target"
result[1,1+nw+nw+5]<- "spent"
result[1,1+nw+nw+6]<- "Reject H0"

result[test+1,1]<- test

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)

for(k in 1:nw){
if(sum(ws[k]==w)>0){cum_events[k]<- sum(cases[ws[k]==w])+sum(controls[ws[k]==w]); cum_cases[k]<- sum(cases[ws[k]==w]); cum_controls[k]<- sum(controls[ws[k]==w])}
if(sum(w_old==ws[k])>0){
cum_events[k]<- cum_events[k]+sum(cases_old[w_old==ws[k]])+sum(controls_old[w_old==ws[k]])
cum_cases[k]<- cum_cases[k]+sum(cases_old[w_old==ws[k]])
cum_controls[k]<- cum_controls[k]+sum(controls_old[w_old==ws[k]]) 
                       }
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&is.na(max(rre[w==ws[k]]))!=TRUE&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- NA}
 result[test+1,1+k]<- as.numeric(cum_events[k]) 
              }

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

if(is.na(CVl)!=TRUE){result[test+1,2*nw+3]<- specify_decimal(CVl/(sum(cum_events*ws)-CVl),2)}else{result[test+1,2*nw+3]<- NA}  # <= low critical value

if(is.na(CVu)!=TRUE){result[test+1,2*nw+4]<- specify_decimal(CVu/(sum(cum_events*ws)-CVu),2)}else{result[test+1,2*nw+4]<- NA}  # <= high critical value

result[test+1,2*nw+5]<- specify_decimal(current_alpha,4)
result[test+1,2*nw+6]<- specify_decimal(actualspent,4)
result[test+1,2*nw+7]<- paste("No")

if(test>1){
for(i in 1:(test-1)){
result[i+1,1]<- i

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)
for(k in 1:nw){
  if(sum(ws[k]==w_old)>0){
               cum_events[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])+sum(controls_old[ws[k]==w_old&test_old<=i])
               cum_cases[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])
               cum_controls[k]<- sum(controls_old[ws[k]==w_old&test_old<=i])
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&is.na(max(rre_old[w_old==ws[k]&test_old==i]))!=TRUE&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(rre_old[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- NA}
               result[i+1,1+k]<- as.numeric(cum_events[k])
                         }else{result[i+1,nw+1+k]<- NA}
    
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVl_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVl_old[i]/(sum(cum_events*ws)-CVl_old[i]),2)}else{result[i+1,2*nw+3]<- NA}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+4]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- NA}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+5]<- NA}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+6]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+6]<- NA}
result[i+1,2*nw+7]<- paste("No")
                    }
          }


             } # close 2

message(  " ",domain = NULL, appendLF = TRUE)
message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
message("=>   Do not reject H0. Proceed to a new test as soon as you have more data.", domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,",zp= ", z_setup,", and M= ",M, ", H0: RR<=",R0, "."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)



                                                                                }# CLOSE




###############
### Situation 3: SampleSize achieved with remaining alpha spending and "start>0"
if(reject>0){actualspent<- 0}
if(totalevents>=SampleSize&alpha-actualspent>0.00000001&start>0&reject==0&reject_new==0){ # OPEN


if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,12))
result[2:(test+1),1]<- seq(1,test,1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases|H0]","RR estimate","LLR","target","actual","CV","Reject H0")

result[2:(test+1),1]<- seq(1,test)
result[test+1,2]<- sum(cases)
result[test+1,3]<- sum(controls)
result[test+1,4]<- sum(cases)+sum(cases_old)
result[test+1,5]<- sum(controls)+sum(controls_old)
result[test+1,6]<- specify_decimal(Exp(test),2)
result[test+1,7]<- specify_decimal(RelRisk(test),2)
result[test+1,8]<- specify_decimal(LLRcum(test),6)
result[test+1,9]<- specify_decimal(current_alpha,4)
result[test+1,10]<- specify_decimal(actualspent,4)
result[test+1,11]<-  CVu
result[test+1,12]<- paste("No")

if(test>1){

for(i in 1:(test-1)){

if(i==1){a1<- 1; a2<- sum(rows_current[1:i])}else{a1<- sum(rows_current[1:(i-1)])+1; a2<- sum(rows_current[1:i])}

result[i+1,2]<- sum(cases_current[a1:a2])
result[i+1,3]<- sum(controls_current[a1:a2])
result[i+1,4]<- sum(cases_current[1:sum(rows_current[1:i])])
result[i+1,5]<- sum(controls_current[1:sum(rows_current[1:i])])
result[i+1,6]<- specify_decimal(Exp(i),2)
result[i+1,7]<- specify_decimal(RelRisk(i),2)
result[i+1,8]<- specify_decimal(LLRcum(i),6)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,9]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,9]<- NA}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,10]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,10]<- NA}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,11]<- CVu_old[i]}else{result[i+1,11]<- NA}
result[i+1,12]<- paste("No")

                    }

          }

###### Graphics

par(mfrow=c(2,2))

## GRAPH FOR CV, OBSERVED, AND EXPECTED CASES

limy<-  3*(sum(cases)+sum(cases_old))/2+10
plot(seq(1,test), rep(limy,test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Cumulative Cases"),main=title,ylim= c(0,limy) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)

points(seq(1,test), result[2:(test+1),4],col="blue",pch=20)
lines(seq(1,test), result[2:(test+1),4],col="blue",lty=1)
points(seq(1,test), result[2:(test+1),6],col="black",pch=2)
lines(seq(1,test), result[2:(test+1),6],col="black",lty=1)

ini<- 1
if(is.na(result[2,11])==TRUE){
while(is.na(result[ini+1,11])==TRUE){ini<- ini+1}; if(ini==test){cvs<- result[ini+1,11]}else{cvs<- rep(0,test-ini+1) ; cvs[1]<- result[ini+1,11]; g1<- 1; for(gg in (ini+1):test){g1<- g1+1; if(is.na(result[gg+1,11])==TRUE){cvs[g1]<-result[gg,11]}else{cvs[g1]<-result[gg+1,11] }}}
cvs<- as.numeric(cvs)
points(seq(ini,test),cvs ,col="red",pch=20)
lines(seq(ini,test), cvs,col="red",lty=1) 
legend("topleft",c("Needed to reject H0 (CV)","Observed","Expected = E[Cases|H0]"),col=c("red","blue","black"),pch=c(18,20,2),lty=c(2,1,1),bty="n")
                                            }else{
legend("topleft",c("Observed","Expected = E[Cases|H0]"),col=c("blue","black"),pch=c(20,2),lty=c(1,1),bty="n")
                                                 }

## GRAPH FOR ALPHA SPEND 

if(sum(is.na(result[2:(test+1),9]))<test){
plot(seq(1,test), rep(as.numeric(result[test+1,9])+0.02,test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Alpha spending"),main="Alpha Spending",ylim= c(0,as.numeric(result[test+1,9])+0.02 ) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(ini,test),result[(ini+1):(test+1),9],lty=2,col="red")
points(seq(ini,test),result[(ini+1):(test+1),9],pch=18,col="red")
lines(seq(ini,test),result[(ini+1):(test+1),10],lty=2,col="blue")
points(seq(ini,test),result[(ini+1):(test+1),10],pch=18,col="blue")
legend("topleft",c("Target","Actual spent"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")
                                           }else{
plot(seq(1,test,1),rep(alpha,test),pch=18,col="white",ylab="Alpha spending",xlab="Test",xaxt="n",cex.lab=1.3, 
sub=" ",font.sub=1,main="Alpha spending"); legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")
                                                } 


## GRAPH FOR RR ESTIMATES

if(sum(result[2:(test+1),7]==Inf)<test){
ini<- 1; while(result[ini+1,7]==Inf){ini<- ini+1}
plot(seq(1,test,1),rep(max(as.numeric(result[(ini+1):(test+1),7])),test),pch=18,col="white",ylab="Observed relative risk",xlab="Test",xaxt="n",cex.lab=1.3, 
ylim=c(0,max(as.numeric(result[(ini+1):(test+1),7]))),sub=" ",font.sub=1,main="Observed Relative Risk")
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(ini,test),result[(ini+1):(test+1),7],lty=2,col="blue")
points(seq(ini,test),result[(ini+1):(test+1),7],pch=18,col="blue")
                                        }else{
plot(seq(1,test,1),rep(1,test),pch=18,col="white",ylab="Observed relative risk",xlab="Test",xaxt="n",cex.lab=1.3, 
sub=" ",font.sub=1,main="Observed relative risk"); legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")
                                             }


## GRAPH FOR LLR

plot(seq(1,test), rep(max(as.numeric(result[2:(test+1),8])),test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Log-likelihood ratio"),main="Log-likelihood ratio",ylim= c(0,max(as.numeric(result[2:(test+1),8])) ) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(1,test),result[2:(test+1),8],lty=2,col="blue")
points(seq(1,test),result[2:(test+1),8],pch=18,col="blue")


             } # close 1


####
if(Tailed==2){ # 2

if(test>1){ws<- c(w_old,w)}else{ws<- w}
ws<- as.numeric(by(ws,ws,max))
nw<- length(ws)
namesw<- rep(0,nw)
for(k in 1:nw){if(k>1){namesw[k]<- paste("w=",ws[k])}else{namesw[k]<- paste("  w=",ws[k])}}

result<- data.frame(matrix(0,test+1,2*nw+7))

names<- rep(0,2*nw+7)
names[1]<- c(" ")

if(nw==1){names[2]<- c("#Events"); names[3]<- c("Relative Risk")}else{
if(nw==2){names[2]<- c("--") ; names[3]<- c("#Events --"); names[4]<- " --Relative"; names[5]<- "Risk--"}else{
if(nw==3){names[2]<- c("--");names[3]<- c("---"); names[4]<- c("#Events"); names[5]<- " --"; names[6]<- "Relative"; names[7]<- "Risk---"}else{
ac<- floor(nw/2+1);  names[2:(ac-1)]<- c("--"); names[ac]<- c("--"); names[ac+1]<- c("#Events "); names[(ac+2):(nw+1)]<- c("--")
if(nw==4){
names[(nw+2)]<- " ---"; names[nw+3]<- c("Relative"); names[nw+4]<- c("Risk"); names[nw+5]<- c("---") 
         }else{
names[(nw+2)]<- " ---"; names[(nw+3):(nw+ac)]<- "--" ;names[nw+ac+1]<- c("Relative"); names[nw+ac+2]<- c("Risk"); names[(nw+ac+3):(2*nw+1)]<- c("--"); names[2*nw+1]<- c("---")
              }          
                                                                                                                              }   
                                                                                                        }
                                                                      }
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "--Critical"; names[1+nw+nw+3]<- "Value--" ; names[1+nw+nw+4]<- "--alpha"; names[1+nw+nw+5]<- "-----" ; names[1+nw+nw+6]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "low"
result[1,1+nw+nw+3]<- "high"
result[1,1+nw+nw+4]<- "target"
result[1,1+nw+nw+5]<- "spent"
result[1,1+nw+nw+6]<- "Reject H0"

result[test+1,1]<- test

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)

for(k in 1:nw){
if(sum(ws[k]==w)>0){cum_events[k]<- sum(cases[ws[k]==w])+sum(controls[ws[k]==w]); cum_cases[k]<- sum(cases[ws[k]==w]); cum_controls[k]<- sum(controls[ws[k]==w])}
if(sum(w_old==ws[k])>0){
cum_events[k]<- cum_events[k]+sum(cases_old[w_old==ws[k]])+sum(controls_old[w_old==ws[k]])
cum_cases[k]<- cum_cases[k]+sum(cases_old[w_old==ws[k]])
cum_controls[k]<- cum_controls[k]+sum(controls_old[w_old==ws[k]]) 
                       }
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"}
 result[test+1,1+k]<- as.numeric(cum_events[k]) 
              }

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

if(CVl!="NA"){result[test+1,2*nw+3]<- specify_decimal(CVl/(sum(cum_events*ws)-CVl),2)}else{result[test+1,2*nw+3]<- "NA"}  # <= low critical value

if(CVu!="NA"){result[test+1,2*nw+4]<- specify_decimal(CVu/(sum(cum_events*ws)-CVu),2)}else{result[test+1,2*nw+4]<-"NA"}  # <= high critical value

result[test+1,2*nw+5]<- specify_decimal(current_alpha,4)
result[test+1,2*nw+6]<- specify_decimal(actualspent,4)
result[test+1,2*nw+7]<- paste("No")

if(test>1){
for(i in 1:(test-1)){
result[i+1,1]<- i

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)
for(k in 1:nw){
  if(sum(ws[k]==w_old)>0){
               cum_events[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])+sum(controls_old[ws[k]==w_old&test_old<=i])
               cum_cases[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])
               cum_controls[k]<- sum(controls_old[ws[k]==w_old&test_old<=i])
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(rre_old[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
               result[i+1,1+k]<- as.numeric(cum_events[k])
                         }else{result[i+1,nw+1+k]<- "NA"}
    
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVl_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVl_old[i]/(sum(cum_events*ws)-CVl_old[i]),2)}else{result[i+1,2*nw+3]<- "NA"}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+4]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+6]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+6]<- "NA"}
result[i+1,2*nw+7]<- paste("No")
                    }
          }


             } # close 2

message(  " ",domain = NULL, appendLF = TRUE)
message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)
message(c("The upper limit on the length of surveillance has been reached."), domain = NULL, appendLF = TRUE) 
message("Then, the 'AlphaSpend' input is no longer used, and by default the target is alpha.")                                                      
message(c("You may now end the sequential analysis without rejecting H0."), domain = NULL, appendLF = TRUE)
message(c("There is still ",specify_decimal(alpha-actualspent,6)," alpha to spend if you wish continue with more analyses."), domain = NULL, appendLF = TRUE)
message(" ", domain = NULL, appendLF = TRUE)  
message("-------------------------------------------------------------------------------------------", domain = NULL, appendLF = TRUE)                                                         
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", zp= ", z_setup,", and M= ",M, ", H0: RR<=",R0, "."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)                                                         
                                            
                                                                  } # CLOSE




###############
### Situation 4: SampleSize achieved without remaining alpha spending

if(totalevents>=SampleSize&alpha-actualspent<=0.00000001&start>0&reject==0&reject_new==0){ # OPEN

if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,12))
result[2:(test+1),1]<- seq(1,test,1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases|H0]","RR estimate","LLR","target","actual","CV","Reject H0")

result[2:(test+1),1]<- seq(1,test)
result[test+1,2]<- sum(cases)
result[test+1,3]<- sum(controls)
result[test+1,4]<- sum(cases)+sum(cases_old)
result[test+1,5]<- sum(controls)+sum(controls_old)
result[test+1,6]<- specify_decimal(Exp(test),2)
result[test+1,7]<- specify_decimal(RelRisk(test),2)
result[test+1,8]<- specify_decimal(LLRcum(test),6)
result[test+1,9]<- specify_decimal(current_alpha,4)
result[test+1,10]<- specify_decimal(actualspent,4)
result[test+1,11]<-  CVu
result[test+1,12]<- paste("No")

if(test>1){

for(i in 1:(test-1)){

if(i==1){a1<- 1; a2<- sum(rows_current[1:i])}else{a1<- sum(rows_current[1:(i-1)])+1; a2<- sum(rows_current[1:i])}

result[i+1,2]<- sum(cases_current[a1:a2])
result[i+1,3]<- sum(controls_current[a1:a2])
result[i+1,4]<- sum(cases_current[1:sum(rows_current[1:i])])
result[i+1,5]<- sum(controls_current[1:sum(rows_current[1:i])])
result[i+1,6]<- specify_decimal(Exp(i),2)
result[i+1,7]<- specify_decimal(RelRisk(i),2)
result[i+1,8]<- specify_decimal(LLRcum(i),6)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,9]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,9]<- NA}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,10]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,10]<- NA}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,11]<- CVu_old[i]}else{result[i+1,11]<- NA}
result[i+1,12]<- paste("No")

                    }

          }

###### Graphics

par(mfrow=c(2,2))

## GRAPH FOR CV, OBSERVED, AND EXPECTED CASES

limy<-  3*(sum(cases)+sum(cases_old))/2+10
plot(seq(1,test), rep(limy,test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Cumulative Cases"),main=title,ylim= c(0,limy) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)

points(seq(1,test), result[2:(test+1),4],col="blue",pch=20)
lines(seq(1,test), result[2:(test+1),4],col="blue",lty=1)
points(seq(1,test), result[2:(test+1),6],col="black",pch=2)
lines(seq(1,test), result[2:(test+1),6],col="black",lty=1)

ini<- 1
if(is.na(result[2,11])==TRUE){
while(is.na(result[ini+1,11])==TRUE){ini<- ini+1}; if(ini==test){cvs<- result[ini+1,11]}else{cvs<- rep(0,test-ini+1) ; cvs[1]<- result[ini+1,11]; g1<- 1; for(gg in (ini+1):test){g1<- g1+1; if(is.na(result[gg+1,11])==TRUE){cvs[g1]<-result[gg,11]}else{cvs[g1]<-result[gg+1,11] }}}
cvs<- as.numeric(cvs)
points(seq(ini,test),cvs ,col="red",pch=20)
lines(seq(ini,test), cvs,col="red",lty=1) 
legend("topleft",c("Needed to reject H0 (CV)","Observed","Expected = E[Cases|H0]"),col=c("red","blue","black"),pch=c(18,20,2),lty=c(2,1,1),bty="n")
                                            }else{
legend("topleft",c("Observed","Expected = E[Cases|H0]"),col=c("blue","black"),pch=c(20,2),lty=c(1,1),bty="n")
                                                 }

## GRAPH FOR ALPHA SPEND 

if(sum(is.na(result[2:(test+1),9]))<test){
plot(seq(1,test), rep(as.numeric(result[test+1,9])+0.02,test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Alpha spending"),main="Alpha Spending",ylim= c(0,as.numeric(result[test+1,9])+0.02 ) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(ini,test),result[(ini+1):(test+1),9],lty=2,col="red")
points(seq(ini,test),result[(ini+1):(test+1),9],pch=18,col="red")
lines(seq(ini,test),result[(ini+1):(test+1),10],lty=2,col="blue")
points(seq(ini,test),result[(ini+1):(test+1),10],pch=18,col="blue")
legend("topleft",c("Target","Actual spent"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")
                                           }else{
plot(seq(1,test,1),rep(alpha,test),pch=18,col="white",ylab="Alpha spending",xlab="Test",xaxt="n",cex.lab=1.3, 
sub=" ",font.sub=1,main="Alpha spending"); legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")
                                                } 



## GRAPH FOR RR ESTIMATES

if(sum(result[2:(test+1),7]==Inf)<test){
ini<- 1; while(result[ini+1,7]==Inf){ini<- ini+1}
plot(seq(1,test,1),rep(max(as.numeric(result[(ini+1):(test+1),7])),test),pch=18,col="white",ylab="Observed relative risk",xlab="Test",xaxt="n",cex.lab=1.3, 
ylim=c(0,max(as.numeric(result[(ini+1):(test+1),7]))),sub=" ",font.sub=1,main="Observed Relative Risk")
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(ini,test),result[(ini+1):(test+1),7],lty=2,col="blue")
points(seq(ini,test),result[(ini+1):(test+1),7],pch=18,col="blue")
                                        }else{
plot(seq(1,test,1),rep(1,test),pch=18,col="white",ylab="Observed relative risk",xlab="Test",xaxt="n",cex.lab=1.3, 
sub=" ",font.sub=1,main="Observed relative risk"); legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")
                                             }


## GRAPH FOR LLR

plot(seq(1,test), rep(max(as.numeric(result[2:(test+1),8])),test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Log-likelihood ratio"),main="Log-likelihood ratio",ylim= c(0,max(as.numeric(result[2:(test+1),8])) ) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(1,test),result[2:(test+1),8],lty=2,col="blue")
points(seq(1,test),result[2:(test+1),8],pch=18,col="blue")


             } # close 1


####
if(Tailed==2){ # 2

if(test>1){ws<- c(w_old,w)}else{ws<- w}
ws<- as.numeric(by(ws,ws,max))
nw<- length(ws)
namesw<- rep(0,nw)
for(k in 1:nw){if(k>1){namesw[k]<- paste("w=",ws[k])}else{namesw[k]<- paste("  w=",ws[k])}}

result<- data.frame(matrix(0,test+1,2*nw+7))

names<- rep(0,2*nw+7)
names[1]<- c(" ")

if(nw==1){names[2]<- c("#Events"); names[3]<- c("Relative Risk")}else{
if(nw==2){names[2]<- c("--") ; names[3]<- c("#Events --"); names[4]<- " --Relative"; names[5]<- "Risk--"}else{
if(nw==3){names[2]<- c("--");names[3]<- c("---"); names[4]<- c("#Events"); names[5]<- " --"; names[6]<- "Relative"; names[7]<- "Risk---"}else{
ac<- floor(nw/2+1);  names[2:(ac-1)]<- c("--"); names[ac]<- c("--"); names[ac+1]<- c("#Events "); names[(ac+2):(nw+1)]<- c("--")
if(nw==4){
names[(nw+2)]<- " ---"; names[nw+3]<- c("Relative"); names[nw+4]<- c("Risk"); names[nw+5]<- c("---") 
         }else{
names[(nw+2)]<- " ---"; names[(nw+3):(nw+ac)]<- "--" ;names[nw+ac+1]<- c("Relative"); names[nw+ac+2]<- c("Risk"); names[(nw+ac+3):(2*nw+1)]<- c("--"); names[2*nw+1]<- c("---")
              }          
                                                                                                                              }   
                                                                                                        }
                                                                      }
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "--Critical"; names[1+nw+nw+3]<- "Value--" ; names[1+nw+nw+4]<- "--alpha"; names[1+nw+nw+5]<- "-----" ; names[1+nw+nw+6]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "low"
result[1,1+nw+nw+3]<- "high"
result[1,1+nw+nw+4]<- "target"
result[1,1+nw+nw+5]<- "spent"
result[1,1+nw+nw+6]<- "Reject H0"

result[test+1,1]<- test

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)

for(k in 1:nw){
if(sum(ws[k]==w)>0){cum_events[k]<- sum(cases[ws[k]==w])+sum(controls[ws[k]==w]); cum_cases[k]<- sum(cases[ws[k]==w]); cum_controls[k]<- sum(controls[ws[k]==w])}
if(sum(w_old==ws[k])>0){
cum_events[k]<- cum_events[k]+sum(cases_old[w_old==ws[k]])+sum(controls_old[w_old==ws[k]])
cum_cases[k]<- cum_cases[k]+sum(cases_old[w_old==ws[k]])
cum_controls[k]<- cum_controls[k]+sum(controls_old[w_old==ws[k]]) 
                       }
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"}
 result[test+1,1+k]<- as.numeric(cum_events[k]) 
              }

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

if(CVl!="NA"){result[test+1,2*nw+3]<- specify_decimal(CVl/(sum(cum_events*ws)-CVl),2)}else{result[test+1,2*nw+3]<- "NA"}  # <= low critical value

if(CVu!="NA"){result[test+1,2*nw+4]<- specify_decimal(CVu/(sum(cum_events*ws)-CVu),2)}else{result[test+1,2*nw+4]<-"NA"}  # <= high critical value

result[test+1,2*nw+5]<- specify_decimal(current_alpha,4)
result[test+1,2*nw+6]<- specify_decimal(actualspent,4)
result[test+1,2*nw+7]<- paste("No")

if(test>1){
for(i in 1:(test-1)){
result[i+1,1]<- i

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)
for(k in 1:nw){
  if(sum(ws[k]==w_old)>0){
               cum_events[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])+sum(controls_old[ws[k]==w_old&test_old<=i])
               cum_cases[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])
               cum_controls[k]<- sum(controls_old[ws[k]==w_old&test_old<=i])
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(rre_old[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
               result[i+1,1+k]<- as.numeric(cum_events[k])
                         }else{result[i+1,nw+1+k]<- "NA"}
    
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVl_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVl_old[i]/(sum(cum_events*ws)-CVl_old[i]),2)}else{result[i+1,2*nw+3]<- "NA"}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+4]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+6]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+6]<- "NA"}
result[i+1,2*nw+7]<- paste("No")
                    }
          }


             } # close 2

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
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", zp= ", z_setup,", and M= ",M, ", H0: RR<=",R0, "."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)                                                         

                                              


                                                                                                   } # CLOSE




###############
### Situation 5: H0 rejected in the current test

if(reject==0&reject_new>0){# OPEN


if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,12))
result[2:(test+1),1]<- seq(1,test,1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases|H0]","RR estimate","LLR","target","actual","CV","Reject H0")

result[2:(test+1),1]<- seq(1,test)
result[test+1,2]<- sum(cases)
result[test+1,3]<- sum(controls)
result[test+1,4]<- sum(cases)+sum(cases_old)
result[test+1,5]<- sum(controls)+sum(controls_old)
result[test+1,6]<- specify_decimal(Exp(test),2)
result[test+1,7]<- specify_decimal(RelRisk(test),2)
result[test+1,8]<- specify_decimal(LLRcum(test),6)
result[test+1,9]<- specify_decimal(current_alpha,4)
result[test+1,10]<- specify_decimal(actualspent,4)
result[test+1,11]<-  CVu
result[test+1,12]<- paste("Yes")

if(test>1){

for(i in 1:(test-1)){

if(i==1){a1<- 1; a2<- sum(rows_current[1:i])}else{a1<- sum(rows_current[1:(i-1)])+1; a2<- sum(rows_current[1:i])}

result[i+1,2]<- sum(cases_current[a1:a2])
result[i+1,3]<- sum(controls_current[a1:a2])
result[i+1,4]<- sum(cases_current[1:sum(rows_current[1:i])])
result[i+1,5]<- sum(controls_current[1:sum(rows_current[1:i])])
result[i+1,6]<- specify_decimal(Exp(i),2)
result[i+1,7]<- specify_decimal(RelRisk(i),2)
result[i+1,8]<- specify_decimal(LLRcum(i),6)
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,9]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,9]<- NA}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,10]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,10]<- NA}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,11]<- CVu_old[i]}else{result[i+1,11]<- NA}
result[i+1,12]<- paste("No")

                    }

          }

###### Graphics

par(mfrow=c(2,2))

## GRAPH FOR CV, OBSERVED, AND EXPECTED CASES

limy<-  3*(sum(cases)+sum(cases_old))/2+10
plot(seq(1,test), rep(limy,test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Cumulative Cases"),main=title,ylim= c(0,limy) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)

points(seq(1,test), result[2:(test+1),4],col="blue",pch=20)
lines(seq(1,test), result[2:(test+1),4],col="blue",lty=1)
points(seq(1,test), result[2:(test+1),6],col="black",pch=2)
lines(seq(1,test), result[2:(test+1),6],col="black",lty=1)

ini<- 1
if(is.na(result[2,11])==TRUE){
while(is.na(result[ini+1,11])==TRUE){ini<- ini+1}; if(ini==test){cvs<- result[ini+1,11]}else{cvs<- rep(0,test-ini+1) ; cvs[1]<- result[ini+1,11]; g1<- 1; for(gg in (ini+1):test){g1<- g1+1; if(is.na(result[gg+1,11])==TRUE){cvs[g1]<-result[gg,11]}else{cvs[g1]<-result[gg+1,11] }}}
cvs<- as.numeric(cvs)
points(seq(ini,test),cvs ,col="red",pch=20)
lines(seq(ini,test), cvs,col="red",lty=1) 
legend("topleft",c("Needed to reject H0 (CV)","Observed","Expected = E[Cases|H0]"),col=c("red","blue","black"),pch=c(18,20,2),lty=c(2,1,1),bty="n")
                                            }else{
legend("topleft",c("Observed","Expected = E[Cases|H0]"),col=c("blue","black"),pch=c(20,2),lty=c(1,1),bty="n")
                                                 }

## GRAPH FOR ALPHA SPEND 

if(sum(is.na(result[2:(test+1),9]))<test){
plot(seq(1,test), rep(as.numeric(result[test+1,9])+0.02,test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Alpha spending"),main="Alpha Spending",ylim= c(0,as.numeric(result[test+1,9])+0.02 ) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(ini,test),result[(ini+1):(test+1),9],lty=2,col="red")
points(seq(ini,test),result[(ini+1):(test+1),9],pch=18,col="red")
lines(seq(ini,test),result[(ini+1):(test+1),10],lty=2,col="blue")
points(seq(ini,test),result[(ini+1):(test+1),10],pch=18,col="blue")
legend("topleft",c("Target","Actual spent"),col=c("red","blue"),pch=c(18,20),lty=c(2,1),bty="n")
                                           }else{
plot(seq(1,test,1),rep(alpha,test),pch=18,col="white",ylab="Alpha spending",xlab="Test",xaxt="n",cex.lab=1.3, 
sub=" ",font.sub=1,main="Alpha spending"); legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")
                                                } 



## GRAPH FOR RR ESTIMATES

if(sum(result[2:(test+1),7]==Inf)<test){
ini<- 1; while(result[ini+1,7]==Inf){ini<- ini+1}
plot(seq(1,test,1),rep(max(as.numeric(result[(ini+1):(test+1),7])),test),pch=18,col="white",ylab="Observed relative risk",xlab="Test",xaxt="n",cex.lab=1.3, 
ylim=c(0,max(as.numeric(result[(ini+1):(test+1),7]))),sub=" ",font.sub=1,main="Observed Relative Risk")
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(ini,test),result[(ini+1):(test+1),7],lty=2,col="blue")
points(seq(ini,test),result[(ini+1):(test+1),7],pch=18,col="blue")
                                        }else{
plot(seq(1,test,1),rep(1,test),pch=18,col="white",ylab="Observed relative risk",xlab="Test",xaxt="n",cex.lab=1.3, 
sub=" ",font.sub=1,main="Observed relative risk"); legend("top",c("Not Applicable"),col=c("red"),pch=c(18),lty=c(2),bty="n")
                                             }


## GRAPH FOR LLR

plot(seq(1,test), rep(max(as.numeric(result[2:(test+1),8])),test), col="white",pch=18,xlab="Test",xaxt="n",cex.main=1.3,cex.lab=1.3,cex.main=1.5,
ylab=c("Log-likelihood ratio"),main="Log-likelihood ratio",ylim= c(0,max(as.numeric(result[2:(test+1),8])) ) )
axis(1, at=seq(1,test), labels=seq(1,test),las=1,cex.axis=1.2)
lines(seq(1,test),result[2:(test+1),8],lty=2,col="blue")
points(seq(1,test),result[2:(test+1),8],pch=18,col="blue")



             } # close 1


####
if(Tailed==2){ # 2

if(test>1){ws<- c(w_old,w)}else{ws<- w}
ws<- as.numeric(by(ws,ws,max))
nw<- length(ws)
namesw<- rep(0,nw)
for(k in 1:nw){if(k>1){namesw[k]<- paste("w=",ws[k])}else{namesw[k]<- paste("  w=",ws[k])}}

result<- data.frame(matrix(0,test+1,2*nw+7))

names<- rep(0,2*nw+7)
names[1]<- c(" ")

if(nw==1){names[2]<- c("#Events"); names[3]<- c("Relative Risk")}else{
if(nw==2){names[2]<- c("--") ; names[3]<- c("#Events --"); names[4]<- " --Relative"; names[5]<- "Risk--"}else{
if(nw==3){names[2]<- c("--");names[3]<- c("---"); names[4]<- c("#Events"); names[5]<- " --"; names[6]<- "Relative"; names[7]<- "Risk---"}else{
ac<- floor(nw/2+1);  names[2:(ac-1)]<- c("--"); names[ac]<- c("--"); names[ac+1]<- c("#Events "); names[(ac+2):(nw+1)]<- c("--")
if(nw==4){
names[(nw+2)]<- " ---"; names[nw+3]<- c("Relative"); names[nw+4]<- c("Risk"); names[nw+5]<- c("---") 
         }else{
names[(nw+2)]<- " ---"; names[(nw+3):(nw+ac)]<- "--" ;names[nw+ac+1]<- c("Relative"); names[nw+ac+2]<- c("Risk"); names[(nw+ac+3):(2*nw+1)]<- c("--"); names[2*nw+1]<- c("---")
              }          
                                                                                                                              }   
                                                                                                        }
                                                                      }
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "--Critical"; names[1+nw+nw+3]<- "Value--" ; names[1+nw+nw+4]<- "--alpha"; names[1+nw+nw+5]<- "-----" ; names[1+nw+nw+6]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "low"
result[1,1+nw+nw+3]<- "high"
result[1,1+nw+nw+4]<- "target"
result[1,1+nw+nw+5]<- "spent"
result[1,1+nw+nw+6]<- "Reject H0"

result[test+1,1]<- test

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)

for(k in 1:nw){
if(sum(ws[k]==w)>0){cum_events[k]<- sum(cases[ws[k]==w])+sum(controls[ws[k]==w]); cum_cases[k]<- sum(cases[ws[k]==w]); cum_controls[k]<- sum(controls[ws[k]==w])}
if(sum(w_old==ws[k])>0){
cum_events[k]<- cum_events[k]+sum(cases_old[w_old==ws[k]])+sum(controls_old[w_old==ws[k]])
cum_cases[k]<- cum_cases[k]+sum(cases_old[w_old==ws[k]])
cum_controls[k]<- cum_controls[k]+sum(controls_old[w_old==ws[k]]) 
                       }
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"}
 result[test+1,1+k]<- as.numeric(cum_events[k]) 
              }

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

if(CVl!="NA"){result[test+1,2*nw+3]<- specify_decimal(CVl/(sum(cum_events*ws)-CVl),2)}else{result[test+1,2*nw+3]<- "NA"}  # <= low critical value

if(CVu!="NA"){result[test+1,2*nw+4]<- specify_decimal(CVu/(sum(cum_events*ws)-CVu),2)}else{result[test+1,2*nw+4]<-"NA"}  # <= high critical value

result[test+1,2*nw+5]<- specify_decimal(current_alpha,4)
result[test+1,2*nw+6]<- specify_decimal(actualspent,4)
result[test+1,2*nw+7]<- paste("Yes")

if(test>1){
for(i in 1:(test-1)){
result[i+1,1]<- i

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)
for(k in 1:nw){
  if(sum(ws[k]==w_old)>0){
               cum_events[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])+sum(controls_old[ws[k]==w_old&test_old<=i])
               cum_cases[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])
               cum_controls[k]<- sum(controls_old[ws[k]==w_old&test_old<=i])
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(rre_old[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
               result[i+1,1+k]<- as.numeric(cum_events[k])
                         }else{result[i+1,nw+1+k]<- "NA"}
    
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVl_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVl_old[i]/(sum(cum_events*ws)-CVl_old[i]),2)}else{result[i+1,2*nw+3]<- "NA"}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+4]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+6]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+6]<- "NA"}
result[i+1,2*nw+7]<- paste("No")
                    }
          }


             } # close 2

message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE) 
                                              message("=>    Reject H0. No further sequential analyses are needed.",domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", zp= ", z_setup,", and M= ",M, ", H0: RR<=",R0, "."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

            }# CLOSE




###############
### Situation 6: H0 rejected in previous tests

if(reject>0){# OPEN


if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,12))
result[2:(test+1),1]<- seq(1,test,1)
colnames(result)<- c(" "," "," ", "-----","Cumulative","--------"," "," ","--alpha","spent--"," "," ") 
result[1,]<- c("Test","Cases","Controls","Cases","Controls","E[Cases|H0]","RR estimate","LLR","target","actual","CV","Reject H0")

result[2:(test+1),1]<- seq(1,test)
result[test+1,2]<- sum(cases)
result[test+1,3]<- sum(controls)
result[test+1,4]<- sum(cases)+sum(cases_old)
result[test+1,5]<- sum(controls)+sum(controls_old)
result[test+1,6]<- specify_decimal(Exp(test),2)
result[test+1,7]<- specify_decimal(RelRisk(test),2)
result[test+1,8]<- specify_decimal(LLRcum(test),6)
result[test+1,9]<- NA
result[test+1,10]<- NA
result[test+1,11]<-  NA
result[test+1,12]<- paste("Yes")

if(test>1){

for(i in 1:(test-1)){

if(i==1){a1<- 1; a2<- sum(rows_current[1:i])}else{a1<- sum(rows_current[1:(i-1)])+1; a2<- sum(rows_current[1:i])}

result[i+1,2]<- sum(cases_current[a1:a2])
result[i+1,3]<- sum(controls_current[a1:a2])
result[i+1,4]<- sum(cases_current[1:sum(rows_current[1:i])])
result[i+1,5]<- sum(controls_current[1:sum(rows_current[1:i])])
result[i+1,6]<- specify_decimal(Exp(i),2)
result[i+1,7]<- specify_decimal(RelRisk(i),2)
result[i+1,8]<- specify_decimal(LLRcum(i),6)
if(i<=reject){
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,9]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,9]<- NA}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,10]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,10]<- NA}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,11]<- CVu_old[i]}else{result[i+1,11]<- NA}
if(i==reject){result[i+1,12]<- paste("Yes")}else{result[i+1,12]<- paste("No")}
            }else{result[i+1,9]<- NA; result[i+1,10]<- NA; result[i+1,11]<- NA ; result[i+1,12]<- paste("Yes")}

                    }

          }

             } # close 1


####
if(Tailed==2){ # 2

if(test>1){ws<- c(w_old,w)}else{ws<- w}
ws<- as.numeric(by(ws,ws,max))
nw<- length(ws)
namesw<- rep(0,nw)
for(k in 1:nw){if(k>1){namesw[k]<- paste("w=",ws[k])}else{namesw[k]<- paste("  w=",ws[k])}}

result<- data.frame(matrix(0,test+1,2*nw+7))

names<- rep(0,2*nw+7)
names[1]<- c(" ")

if(nw==1){names[2]<- c("#Events"); names[3]<- c("Relative Risk")}else{
if(nw==2){names[2]<- c("--") ; names[3]<- c("#Events --"); names[4]<- " --Relative"; names[5]<- "Risk--"}else{
if(nw==3){names[2]<- c("--");names[3]<- c("---"); names[4]<- c("#Events"); names[5]<- " --"; names[6]<- "Relative"; names[7]<- "Risk---"}else{
ac<- floor(nw/2+1);  names[2:(ac-1)]<- c("--"); names[ac]<- c("--"); names[ac+1]<- c("#Events "); names[(ac+2):(nw+1)]<- c("--")
if(nw==4){
names[(nw+2)]<- " ---"; names[nw+3]<- c("Relative"); names[nw+4]<- c("Risk"); names[nw+5]<- c("---") 
         }else{
names[(nw+2)]<- " ---"; names[(nw+3):(nw+ac)]<- "--" ;names[nw+ac+1]<- c("Relative"); names[nw+ac+2]<- c("Risk"); names[(nw+ac+3):(2*nw+1)]<- c("--"); names[2*nw+1]<- c("---")
              }          
                                                                                                                              }   
                                                                                                        }
                                                                      }
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "--Critical"; names[1+nw+nw+3]<- "Value--" ; names[1+nw+nw+4]<- "--alpha"; names[1+nw+nw+5]<- "-----" ; names[1+nw+nw+6]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "low"
result[1,1+nw+nw+3]<- "high"
result[1,1+nw+nw+4]<- "target"
result[1,1+nw+nw+5]<- "spent"
result[1,1+nw+nw+6]<- "Reject H0"

result[test+1,1]<- test

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)

for(k in 1:nw){
if(sum(ws[k]==w)>0){cum_events[k]<- sum(cases[ws[k]==w])+sum(controls[ws[k]==w]); cum_cases[k]<- sum(cases[ws[k]==w]); cum_controls[k]<- sum(controls[ws[k]==w])}
if(sum(w_old==ws[k])>0){
cum_events[k]<- cum_events[k]+sum(cases_old[w_old==ws[k]])+sum(controls_old[w_old==ws[k]])
cum_cases[k]<- cum_cases[k]+sum(cases_old[w_old==ws[k]])
cum_controls[k]<- cum_controls[k]+sum(controls_old[w_old==ws[k]]) 
                       }
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"}
 result[test+1,1+k]<- as.numeric(cum_events[k]) 
              }

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

result[test+1,2*nw+3]<- "NA"  # <= low critical value

result[test+1,2*nw+4]<- "NA"  # <= high critical value

result[test+1,2*nw+5]<- "NA"
result[test+1,2*nw+6]<- "NA"
result[test+1,2*nw+7]<- paste("Yes")

if(test>1){
for(i in 1:(test-1)){
result[i+1,1]<- i

cum_events<- rep(0,nw)
cum_cases<- rep(0,nw)
cum_controls<- rep(0,nw)
for(k in 1:nw){
  if(sum(ws[k]==w_old)>0){
               cum_events[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])+sum(controls_old[ws[k]==w_old&test_old<=i])
               cum_cases[k]<- sum(cases_old[ws[k]==w_old&test_old<=i])
               cum_controls[k]<- sum(controls_old[ws[k]==w_old&test_old<=i])
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(rre_old[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
               result[i+1,1+k]<- as.numeric(cum_events[k])
                         }else{result[i+1,nw+1+k]<- "NA"}
    
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVl_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVl_old[i]/(sum(cum_events*ws)-CVl_old[i]),2)}else{result[i+1,2*nw+3]<- "NA"}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+4]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+6]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+6]<- "NA"}
if(i<=reject){result[i+1,2*nw+7]<- paste("No");if(i==reject){result[i+1,2*nw+7]<- paste("Yes")}}else{result[i+1,2*nw+7]<- paste("Yes");result[i+1,2*nw+6]<- "NA";result[i+1,2*nw+5]<- "NA";result[i+1,2*nw+4]<- "NA";result[i+1,2*nw+3]<- "NA"}
                    }
          }


             } # close 2

message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)      
          message(paste(c("=>    H0 was rejected on test"," ",reject,". ","No further sequential analyses are needed.")),domain = NULL, appendLF = TRUE)  
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", zp= ", z_setup,", and M= ",M, ", H0: RR<=",R0, "."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

            }# CLOSE




############################################################
## UPDATING INFORMATION FOR FUTURE TESTES
############################################################

                                                                                                                                       
### For decision matters

if(test> ncol(inputSetUp) | sum(row_old)+length(z_w_levels)>ncol(inputSetUp)){inputSetUp<- cbind(inputSetUp,matrix(0,nrow(inputSetUp),length(z_w_levels)))}
inputSetUp[1,1]<- test
inputSetUp[2,1]<- start
if(start>0){inputSetUp[1,7]<- reject_new}
if(reject==0){inputSetUp[3,test]<- CVu}
if(reject==0){if(Tailed==2){inputSetUp[9,test]<- CVl}}
if(reject>0){inputSetUp[5,test]<- 0}else{inputSetUp[5,test]<- actualspent}
if(reject==0&start>0){inputSetUp[7,test]<- current_alpha}; if(start==0){inputSetUp[7,test]<- 0}
inputSetUp[10,test]<- length(z_w_levels)
inputSetUp[4,(sum(row_old)+1):(sum(row_old)+length(z_w_levels))]<- data_matrix[,3]
inputSetUp[8,(sum(row_old)+1):(sum(row_old)+length(z_w_levels))]<- data_matrix[,4]
inputSetUp[11,(sum(row_old)+1):(sum(row_old)+length(z_w_levels))]<- data_matrix[,1]
inputSetUp[12,1:length(alphaspend)]<- alphaspend
inputSetUp[13,(sum(row_old)+1):(sum(row_old)+length(z_w_levels))]<- data_matrix[,2]
inputSetUp[14,(sum(row_old)+1):(sum(row_old)+length(z_w_levels))]<- test

## inputSetUp matrix contains:
# line 1: (C11) the index for the order of the test (zero entry if we did not have a first test yet), (C12) SampleSize, (C13) alpha, (C14) M, (C15) base(the line of p where the looping will start in the next test), (C16) title, (C17) reject (the index indicating if and when H0 was rejected), (C18) rho (zero if Wald is used)
# line 2: says if the analysis has been started or not. 0 for not started. 1 for already started.
# line 3: CVu which is the upper critical values in the scale of Sn
# line 4: observed cases  
# line 5: actual alpha spent
# line 6: cumulative expected value of Sn E[sum(w*C)] under H0, test by test
# line 7: has the target alpha spending until the (test-1)th look.
# line 8: observed controls
# line 9: CVl which is the lower critical values in the scale of Sn
# line 10: the number of different weights per test
# line 11: matched case-control (old z)
# line 12: Target alpha spending defined with AnalyzeSetUp
# line 13: weight for each observation (old w)
# line 14: test associated to each weight



############################################################
## SAVING INFORMATION FOR FUTURE TESTS
############################################################

write.table(inputSetUp,name)

if(sum(cases_old)+sum(controls_old)+sum(cases)+sum(controls)>=M){
write.table(pold,paste(name1,"pold.txt",sep=""))
write.table(statesOld,paste(name1,"statesOld.txt",sep=""))
                                                                }

result2<- result[2:(test+1),]
colnames(result2)<- c("Test", "Cases", "Controls", "Cum. cases", "Cum. controls", "E[Cases|H0]" , "RR estimate", "LLR", "target alpha", "actual alpha", "CV", "Reject H0")
invisible(result2)

#####################################
}##### Close function Analyze.wBinomial
#####################################



#### Example. 
####
#    Note: cut off the symbol "#" before runing the lines below.

#AnalyzeSetUp.Binomial(name="test",N="n", alpha=0.05, pp=0.44, M=5, AlphaSpendType="optimal",address="C:/Users/ivair/POST-DOC/Material para construcao do pacote Sequential/Sequential_3.0/PASTA PARA TREINO")
#system.time(res<- Analyze.Binomial(name="test",test=1,p=0.44,cases=6,controls=2))








