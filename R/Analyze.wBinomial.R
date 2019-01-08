
Analyze.wBinomial<- function(name,test,z=1,w,ExposureA,ExposureB,AlphaSpend="n")
{


if( sum(is.numeric(z))==1){z<- max(z)}

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



cases<- ExposureA
controls<- ExposureB


if( sum(is.numeric(cases))!=1){stop("Symbols and texts are not applicable for 'ExposureA'. It must be a vector of counts.",call. =FALSE)}
if( sum(is.numeric(controls))!=1){stop("Symbols and texts are not applicable for 'ExposureB'. It must be a vector of counts.",call. =FALSE)}
if( sum(is.numeric(w))!=1){stop("Symbols and texts are not applicable for 'w'. It must be a vector of positive numbers.",call. =FALSE)}
if( sum(is.numeric(z))!=1){stop("Symbols and texts are not applicable for 'z'. It must be a vector of positive numbers.",call. =FALSE)}


if(sum(cases<0)>0){stop("ExposureA must contain only zeros and positive integers.",call. =FALSE)}
if(sum(controls<0)>0){stop("ExposureB must contain only zeros and positive integers.",call. =FALSE)}
if(sum(w<0)>0){stop("w must contain only positive numbers.",call. =FALSE)}
if(sum(z<0)>0){stop("z must contain only zeros and positive numbers.",call. =FALSE)}



if( sum(is.numeric(test))!=1|length(test)>1){stop("Symbols and texts are not applicable for 'test'. It must be an integer greater than zero.",call. =FALSE)}

if(test<0){stop("'test' must be an integer greater than zero.",call. =FALSE)}

if(sum(cases)+sum(controls)<=0){stop("It must be informed a number greater than zero for the total number of events.",call. =FALSE)}

if( length(names(table(c(length(cases),length(controls),length(w)))))>1 ){stop("ExposureA, ExposureB, z, and w must be of the same length.",call. =FALSE)}

w<- w[order(w)]; cases<- cases[order(w)]; controls<- controls[order(w)]

####
## Uploading information from previous tests
####

inputSetUp<- read.table(name)
if(inputSetUp[1,1]!=test-1){stop(c("The current test should be"," ",inputSetUp[1,1]+1,". ",
"If you do not have information about previous tests, see the user manual for more details."),call. =FALSE)}



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
# line 10: the number of rows (different weights/z ) per test
# line 11: matched case-control (old z)
# line 12: PENSAR
# line 13: weight for each observation (old w)
# line 14: test associated to each weight 

#### 

SampleSize<- inputSetUp[1,2]
alpha<- inputSetUp[1,3]
M<- inputSetUp[1,4]
start<- inputSetUp[2,1]
reject<- inputSetUp[1,7]
rho<- inputSetUp[1,8]
Tailed<- inputSetUp[1,9]
maxspent<- max(inputSetUp[5,])
actual_alpha_old<- inputSetUp[5,]
target_alpha_old<- inputSetUp[7,]
row_old<- inputSetUp[10,] 
CVu_old<- as.numeric(inputSetUp[3,]) 
cases_old<- as.numeric(inputSetUp[4,1:sum(row_old)])  # <= cases_old  
CVl_old<- as.numeric(inputSetUp[9,])  # <= colocar CVl aqui
controls_old<- as.numeric(inputSetUp[8,1:sum(row_old)]) # <= controls_old  
z_old<- as.numeric(inputSetUp[11,1:sum(row_old)])        # z_old   
#z2_old<- inputSetUp[12,]        # 
w_old<- as.numeric(inputSetUp[13,1:sum(row_old)])         # w_old   
test_old<- as.numeric(inputSetUp[14,1:sum(row_old)])      # vector with indexes for testing



#### More checks

if( sum(is.numeric(AlphaSpend))!=1&AlphaSpend!="n"){stop("Symbols and texts are not applicable for 'AlphaSpend'. If you want to use the default, use 'n'. Otherwise,  'AlphaSpend' must be a positive number smaller than or equal to 'alpha'.",call. =FALSE)}

if( sum(length(AlphaSpend))!=1){stop("'AlphaSpend' must be a single value, not a vector.",call. =FALSE)}

if(AlphaSpend<0&AlphaSpend!="n"){stop("'AlphaSpend' must be a positive number smaller than 'alpha'.",call.=FALSE)}

if(AlphaSpend>alpha&AlphaSpend!="n"){stop(c("'AlphaSpend' must be smaller than or equal to ",alpha,"."),call. =FALSE)}


if(start>0){
name2<- paste(name1,"pold.txt",sep="")
pold<- matrix(as.numeric(read.table(name2)[[1]]), ,1) # has the previous state probabilities 
name2<- paste(name1,"statesOld.txt",sep="")
statesOld<- as.numeric(read.table(name2)[[1]]) # has the possible state of Sn until (test-1)th analysis.
           }

##################################################################
## Setting up the target alpha spending.

kref<- sum(cases+controls)+sum(controls_old+cases_old); events<- kref

if(SampleSize<=kref){current_alpha<- alpha}else{
if(max(inputSetUp[5,])<alpha-0.00000001&inputSetUp[1,7]==0&kref>=M){

if(AlphaSpend=="n"){
  if(inputSetUp[1,8]==0){ # for Wald type of alpha spending
                         if(inputSetUp[1,2]<as.numeric(events+sum(inputSetUp[4,]))){
                         current_alpha<- alpha
                                                                                   }else{current_alpha<- as.numeric(target_alpha_old[events+sum(inputSetUp[4,])]) }
                        }else{
                         current_alpha<- alpha*(kref/inputSetUp[1,2])^rho
                             }
                   }else{
                         if(AlphaSpend<=max(inputSetUp[5,])&test>1){stop(c("For this test, 'AlphaSpend' must be selected in the (", specify_decimal(max(inputSetUp[5,]),6), ",", alpha,"] interval because it has already been spent up to ",specify_decimal(max(inputSetUp[5,]),6)," until the previous test."),call. =FALSE)}
                         current_alpha<- AlphaSpend
                        }

                                                                            }
                                                                         }


############################################################
###### INTERNAL AUXILIARY FUNCTIONS
############################################################


# Function to calculate the marginal distribution of Sn for an unique group of data
##########################################

sn_marginal<- function(z,w,cases,controls)
{
z_w_levels<- names(table(paste(z,w)))
data_matrix<- matrix(0,length(z_w_levels),5)
colnames(data_matrix)<- c("z","w","Cases","Controls","Total")

for(ii in 1:length(z_w_levels)){ # 1
sum_cases<- sum(cases[z_w_levels[ii]==paste(z,w)]) 
sum_controls<- sum(controls[z_w_levels[ii]==paste(z,w)])
sum_total<- sum_cases+sum_controls

data_matrix[ii,]<- c(max(z[z_w_levels[ii]==paste(z,w)]),max(w[z_w_levels[ii]==paste(z,w)]),sum_cases,sum_controls,sum_cases+sum_controls)

pii<- 1/(1+max(z[z_w_levels[ii]==paste(z,w)]))
 
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

probs_b<- probs

         } # close 2

states_b<- states                          
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

if(start==0){#1
     res<- sn_marginal(z,w,casesr,controlsr)
     statesNew<- res[[1]]
     pnew<- res[[2]]
            }#1 close if start==0



if(start>0){#2

resm<- sn_marginal(z,w,casesr,controlsr)
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

          }#2 close if start>0


perrorI<- current_alpha-maxspent 

if(Tailed==1){ # 3

jj<- length(pnew) ; i<- jj ; p<- 0 ; pp<-0
while(p<perrorI){p<- sum(pnew[length(pnew):jj]); if(p<perrorI){i<- jj;jj<- jj-1;pp<- p}}


alphahere<- pp              ## actual alpha spent in the current test 
if(pp>0){CVf<- statesNew[i]}else{CVf<- "NA"}
statesOld<- statesNew[1:i]
pold<- pnew[1:i]

return(list(CVf,alphahere,pold,statesOld))
            } # 3

if(Tailed==2){ # 4

jj<- length(pnew) ; iu<- jj ; p<- 0 ; pp1<-0 
while(p< perrorI/2){p<- sum(pnew[length(pnew):jj]); if(p< perrorI/2){iu<- jj;jj<- jj-1;pp1<- p}}

if(pp1>0){CVfu<- statesNew[iu]}else{CVfu<- "NA"}

jj<- 1 ; il<- 1 ; p<- 0 ; pp2<-0
while(p< perrorI/2){p<- sum(pnew[1:jj]); if(p< perrorI/2){il<- jj;jj<- jj+1;pp2<- p}}

if(pp2>0){CVfl<- statesNew[il]}else{CVfl<- "NA"}

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
zhh<- rep(z,length(w))
res<- critical_value(z=zhh,w,casesr=cases,controlsr=controls,current_alpha,maxspent)

if(Tailed==1){
CVu<- res[[1]]
actualspent<- res[[2]]+max(actual_alpha_old)
if(actualspent==0){CVu<- "NA"}
pold<- res[[3]]
statesOld<- res[[4]]
             }

if(Tailed==2){
CVu<- res[[1]]
CVl<- res[[2]]
actualspent<- res[[3]]+max(actual_alpha_old)
if(actualspent==0){CVu<- "NA"; CVl<- "NA"}
pold<- res[[4]]
statesOld<- res[[5]]
             }

# Surveillance started?
if(start>0){Sn<- sum(cases*w)+sum(cases_old*w_old)}else{start<- test; Sn<- sum(cases*w)}
             

# H0 rejected?
if(Tailed==1){
if(Sn>=CVu){reject_new<- test}else{reject_new<- 0}
             }else{
if( (Sn>=CVu&CVu!="NA")| (Sn<=CVl&CVl!="NA") ){reject_new<- test}else{reject_new<- 0}
                  }
                                                                     }else{reject_new<- max(0,reject)}

if(M>totalevents |reject==1| max(inputSetUp[5,])>=alpha-0.00000001 ){actualspent<- 0; CVu<- "NA";  if(Tailed==2){CVl<- "NA"} }



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



if(test>1){

w_levels<- names(table(c(w_old,w)))

dataCase<- c(cases_old,cases); dataControl<- c(controls_old,controls); dataW<- c(w_old,w); dataTest<- c(test_old,rep(test,length(cases)))

dataMat<- data.frame(matrix(0,length(dataCase),4)); dataMat[,1]<- dataTest; dataMat[,2]<- dataW; dataMat[,3]<- dataCase; dataMat[,4]<- dataControl

colnames(dataMat)<- c("dataTest","dataW","dataCase","dataControl")
#attach(dataMat,warn.conflicts=FALSE)

dataRRe<- matrix(0,test,length(w_levels))
for(i in 1:test){
  for(j in 1:length(w_levels)){
         if(sum(dataCase[which(dataTest<=i&dataW==w_levels[j])]+dataControl[which(dataTest<=i&dataW==w_levels[j])])>0){
                             casRRe<- sum(dataCase[which(dataTest<=i&dataW==w_levels[j])])
                             contRRe<- sum(dataControl[which(dataTest<=i&dataW==w_levels[j])])        
                             dataRRe[i,j]<- z*casRRe/contRRe 
                                                                                                      }
                              }
                }
  rre<- dataRRe[test,]
  rre_old<- dataRRe[1:(test-1),]
          }else{rre<- z*cases/controls} 




#rre<- z*cases/controls 
#if(test>1){rre_old<- z_old*cases_old/controls_old}


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

result<- data.frame(matrix(0,test+1,2*nw+6))

names<- rep(0,2*nw+6)
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
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "Critical"; names[1+nw+nw+3]<- "--alpha"; names[1+nw+nw+4]<- "-----" ; names[1+nw+nw+5]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "Value"
result[1,1+nw+nw+3]<- "target"
result[1,1+nw+nw+4]<- "spent"
result[1,1+nw+nw+5]<- "Reject H0"

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
result[test+1,1+k]<- as.numeric(cum_events[k])
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"} 
}

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

result[test+1,2*nw+3]<- "NA"  # <= critical value

result[test+1,2*nw+4]<- 0
result[test+1,2*nw+5]<- 0
result[test+1,2*nw+6]<- paste("No")


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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
                         }else{result[i+1,nw+1+k]<- "NA"}
    result[i+1,1+k]<- as.numeric(cum_events[k])
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+4]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
result[i+1,2*nw+6]<- paste("No")
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
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"}
 result[test+1,1+k]<- as.numeric(cum_events[k]) 
              }

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

result[test+1,2*nw+3]<- "NA"  # <= low critical value

result[test+1,2*nw+4]<- "NA"  # <= high critical value

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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
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
                                              message("=>    H0 cannot be rejected yet because the cumulative events is still smaller than M.",domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho, ", Tailed= ", Tailed,", and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)


           } # CLOSE




###############
### Situation 2: Surveillance started, H0 not rejected yet and sample size not achieved

if(reject==0&reject_new==0&start>0&totalevents<SampleSize){# OPEN

if(Tailed==1){ # 1

if(test>1){ws<- c(w_old,w)}else{ws<- w}
ws<- as.numeric(by(ws,ws,max))
nw<- length(ws)
namesw<- rep(0,nw)
for(k in 1:nw){if(k>1){namesw[k]<- paste("w=",ws[k])}else{namesw[k]<- paste("  w=",ws[k])}}

result<- data.frame(matrix(0,test+1,2*nw+6))

names<- rep(0,2*nw+6)
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
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "Critical"; names[1+nw+nw+3]<- "--alpha"; names[1+nw+nw+4]<- "-----" ; names[1+nw+nw+5]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "Value"
result[1,1+nw+nw+3]<- "target"
result[1,1+nw+nw+4]<- "spent"
result[1,1+nw+nw+5]<- "Reject H0"

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
result[test+1,1+k]<- as.numeric(cum_events[k])
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"} 
}

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls
  
if(CVu!="NA"){result[test+1,2*nw+3]<- specify_decimal(CVu/(sum(cum_events*ws)-CVu),2)}else{result[test+1,2*nw+3]<-"NA"} # <= critical value

result[test+1,2*nw+4]<- specify_decimal(current_alpha,4)
result[test+1,2*nw+5]<- specify_decimal(actualspent,4)
result[test+1,2*nw+6]<- paste("No")


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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
                         }else{result[i+1,nw+1+k]<- "NA"}
    result[i+1,1+k]<- as.numeric(cum_events[k])
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+4]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
result[i+1,2*nw+6]<- paste("No")
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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
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
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
message("=>   Do not reject H0. Proceed to a new test as soon as you have more data.", domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
print(result,right=TRUE,row.names=FALSE)

message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", Tailed= ", Tailed,", and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)


                                                                                }# CLOSE




###############
### Situation 3: SampleSize achieved with remaining alpha spending and "start>0"
if(reject>0){actualspent<- 0}
if(totalevents>=SampleSize&alpha-actualspent>0.00000001&start>0&reject==0&reject_new==0){ # OPEN


if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,2*nw+6))

names<- rep(0,2*nw+6)
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
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "Critical"; names[1+nw+nw+3]<- "--alpha"; names[1+nw+nw+4]<- "-----" ; names[1+nw+nw+5]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "Value"
result[1,1+nw+nw+3]<- "target"
result[1,1+nw+nw+4]<- "spent"
result[1,1+nw+nw+5]<- "Reject H0"

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
result[test+1,1+k]<- as.numeric(cum_events[k])
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"} 
}

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

if(CVu!="NA"){result[test+1,2*nw+3]<- specify_decimal(CVu/(sum(cum_events*ws)-CVu),2)}else{result[test+1,2*nw+3]<-"NA"} # <= critical value

result[test+1,2*nw+4]<- specify_decimal(current_alpha,4)
result[test+1,2*nw+5]<- specify_decimal(actualspent,4)
result[test+1,2*nw+6]<- paste("No")


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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
                         }else{result[i+1,nw+1+k]<- "NA"}
    result[i+1,1+k]<- as.numeric(cum_events[k])
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+4]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
result[i+1,2*nw+6]<- paste("No")
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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
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
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", Tailed= ", Tailed,", and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)                                                         
                                            
                                                                  } # CLOSE




###############
### Situation 4: SampleSize achieved without remaining alpha spending

if(totalevents>=SampleSize&alpha-actualspent<=0.00000001&start>0&reject==0&reject_new==0){ # OPEN

if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,2*nw+6))

names<- rep(0,2*nw+6)
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
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "Critical"; names[1+nw+nw+3]<- "--alpha"; names[1+nw+nw+4]<- "-----" ; names[1+nw+nw+5]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "Value"
result[1,1+nw+nw+3]<- "target"
result[1,1+nw+nw+4]<- "spent"
result[1,1+nw+nw+5]<- "Reject H0"

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
result[test+1,1+k]<- as.numeric(cum_events[k])
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"} 
}

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

if(CVu!="NA"){result[test+1,2*nw+3]<- specify_decimal(CVu/(sum(cum_events*ws)-CVu),2)}else{result[test+1,2*nw+3]<-"NA"} # <= critical value

result[test+1,2*nw+4]<- specify_decimal(current_alpha,4)
result[test+1,2*nw+5]<- specify_decimal(actualspent,4)
result[test+1,2*nw+6]<- paste("No")


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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
                         }else{result[i+1,nw+1+k]<- "NA"}
    result[i+1,1+k]<- as.numeric(cum_events[k])
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+4]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
result[i+1,2*nw+6]<- paste("No")
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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
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
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", Tailed= ", Tailed,", and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)                                                         

                                              


                                                                                                   } # CLOSE




###############
### Situation 5: H0 rejected in the current test

if(reject==0&reject_new>0){# OPEN


if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,2*nw+6))

names<- rep(0,2*nw+6)
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
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "Critical"; names[1+nw+nw+3]<- "--alpha"; names[1+nw+nw+4]<- "-----" ; names[1+nw+nw+5]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "Value"
result[1,1+nw+nw+3]<- "target"
result[1,1+nw+nw+4]<- "spent"
result[1,1+nw+nw+5]<- "Reject H0"

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
result[test+1,1+k]<- as.numeric(cum_events[k])
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"} 
}

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

if(CVu!="NA"){result[test+1,2*nw+3]<- specify_decimal(CVu/(sum(cum_events*ws)-CVu),2)}else{result[test+1,2*nw+3]<-"NA"} # <= critical value

result[test+1,2*nw+4]<- specify_decimal(current_alpha,4)
result[test+1,2*nw+5]<- specify_decimal(actualspent,4)
result[test+1,2*nw+6]<- paste("Yes")


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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
                         }else{result[i+1,nw+1+k]<- "NA"}
    result[i+1,1+k]<- as.numeric(cum_events[k])
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+4]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
result[i+1,2*nw+6]<- paste("No")
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
##################MODIFICAR PARA 2 no lugar de 8 abaixo
result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),8) # <= test statistic SA/SB, or, S_cases/S_controls

if(CVl!="NA"){result[test+1,2*nw+3]<- specify_decimal(CVl/(sum(cum_events*ws)-CVl),8)}else{result[test+1,2*nw+3]<- "NA"}  # <= low critical value

if(CVu!="NA"){result[test+1,2*nw+4]<- specify_decimal(CVu/(sum(cum_events*ws)-CVu),8)}else{result[test+1,2*nw+4]<-"NA"}  # <= high critical value

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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
               result[i+1,1+k]<- as.numeric(cum_events[k])
                         }else{result[i+1,nw+1+k]<- "NA"}
    
              }
##################MODIFICAR PARA 2 no lugar de 8 abaixo
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),8)
if(is.numeric(CVl_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVl_old[i]/(sum(cum_events*ws)-CVl_old[i]),8)}else{result[i+1,2*nw+3]<- "NA"}
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+4]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),8)}else{result[i+1,2*nw+4]<- "NA"}
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
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", Tailed= ", Tailed,", and M= ",M,"."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)

            }# CLOSE




###############
### Situation 6: H0 rejected in previous tests

if(reject>0){# OPEN


if(Tailed==1){ # 1

result<- data.frame(matrix(0,test+1,2*nw+6))

names<- rep(0,2*nw+6)
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
names[1+nw+nw+1]<- "Test";names[1+nw+nw+2]<- "Critical"; names[1+nw+nw+3]<- "--alpha"; names[1+nw+nw+4]<- "-----" ; names[1+nw+nw+5]<- " "
colnames(result)<- names


result[1,1]<- "Test"
result[1,2:(1+nw)]<- namesw
result[1,(2+nw):(1+nw+nw)]<- namesw
result[1,1+nw+nw+1]<- "Statistic"
result[1,1+nw+nw+2]<- "Value"
result[1,1+nw+nw+3]<- "target"
result[1,1+nw+nw+4]<- "spent"
result[1,1+nw+nw+5]<- "Reject H0"

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
result[test+1,1+k]<- as.numeric(cum_events[k])
if(sum(ws[k]==w)>0&max(rre[w==ws[k]])!="NaN"&max(rre[w==ws[k]])!="NA"&max(rre[w==ws[k]])!=Inf){result[test+1,nw+1+k]<- specify_decimal(mean(rre[w==ws[k]]),2)}else{result[test+1,nw+1+k]<- "NA"} 
}

result[test+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2) # <= test statistic SA/SB, or, S_cases/S_controls

result[test+1,2*nw+3]<- "NA"  # <= critical value

result[test+1,2*nw+4]<- "NA"
result[test+1,2*nw+5]<- "NA"
result[test+1,2*nw+6]<- paste("Yes")


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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
                         }else{result[i+1,nw+1+k]<- "NA"}
    result[i+1,1+k]<- as.numeric(cum_events[k])
              }
result[i+1,2*nw+2]<- specify_decimal(sum(cum_cases*ws)/sum(cum_controls*ws),2)
if(is.numeric(CVu_old[i])==TRUE){result[i+1,2*nw+3]<- specify_decimal(CVu_old[i]/(sum(cum_events*ws)-CVu_old[i]),2)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(target_alpha_old[[i]])==TRUE){result[i+1,2*nw+4]<- specify_decimal(target_alpha_old[i],4)}else{result[i+1,2*nw+4]<- "NA"}
if(is.numeric(actual_alpha_old[[i]])==TRUE){result[i+1,2*nw+5]<- specify_decimal(actual_alpha_old[i],4)}else{result[i+1,2*nw+5]<- "NA"}
if(i<reject){result[i+1,2*nw+6]<- paste("No")}else{result[i+1,2*nw+6]<- paste("Yes");result[i+1,2*nw+5]<- "NA";result[i+1,2*nw+4]<- "NA";result[i+1,2*nw+3]<- "NA"}
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
               if(sum(w_old==ws[k]&test_old==i)>0&max(rre_old[w_old==ws[k]&test_old==i])!="NaN"&max(rre_old[w_old==ws[k]&test_old==i])!="NA"&max(rre_old[w_old==ws[k]&test_old==i])!=Inf){result[i+1,nw+1+k]<- specify_decimal(mean(c(t(rre_old))[w_old==ws[k]&test_old==i]),2)}else{result[i+1,nw+1+k]<- "NA"}
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
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", Tailed= ", Tailed,", and M= ",M,"."),domain = NULL, appendLF = TRUE)
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
if(reject==0){inputSetUp[7,test]<- current_alpha}
inputSetUp[10,test]<- length(z_w_levels)
inputSetUp[4,(sum(row_old)+1):(sum(row_old)+length(z_w_levels))]<- data_matrix[,3]
inputSetUp[8,(sum(row_old)+1):(sum(row_old)+length(z_w_levels))]<- data_matrix[,4]
inputSetUp[11,(sum(row_old)+1):(sum(row_old)+length(z_w_levels))]<- data_matrix[,1]
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
# line 12: PENSAR
# line 13: weight for each observation (old w)
# line 14: test associated to each weight



############################################################
## SAVING INFORMATION FOR FUTURE TESTES
############################################################

write.table(inputSetUp,name)

if(start>0){
write.table(pold,paste(name1,"pold.txt",sep=""))
write.table(statesOld,paste(name1,"statesOld.txt",sep=""))
           }

invisible(result)

#####################################
}##### Close function Analyze.wBinomial
#####################################

