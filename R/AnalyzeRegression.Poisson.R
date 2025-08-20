

AnalyzeRegression.Poisson<- function(name,test,mu0="n",events,covariates,AlphaSpend="n")
{

lambdas<- 1 # baseline rate per event defined in Silva et al.(2025). It is part of mu0 in this code.  


### Important checks

if( sum(is.numeric(events))!=1){stop("Symbols and texts are not applicable for 'events'. It must be an integer number or zero.",call. =FALSE)}

if( sum(is.numeric(test))!=1|length(test)>1){stop("Symbols and texts are not applicable for 'test'. It must be an integer greater than zero.",call. =FALSE)}

if(test<0){stop("'test' must be an integer greater than zero.",call. =FALSE)}

if(sum(events<0)>0){stop("The count 'events' must be an integer greater than or equal to zero.",call. =FALSE)}
 
if(sum(is.numeric(mu0))!=1){stop("Symbols and texts are not applicable for 'mu0'. It must be a number greater than zero.",call. =FALSE)}

if(sum(mu0<=0)>0){stop("The entry of 'mu0' must be a number greater than zero.",call. =FALSE)}

if(length(events)!=length(mu0)){stop("'events' and 'mu0' are vectors that must have the same number of entries.",call. =FALSE)}


#### More checks

if( sum(is.numeric(AlphaSpend))!=1&AlphaSpend!="n"){stop("Symbols and texts are not applicable for 'AlphaSpend'. If you want to use the default, use 'n'. Otherwise,  'AlphaSpend' must be a positive number smaller than or equal to 'alpha'.",call. =FALSE)}

if( sum(length(AlphaSpend))!=1){stop("'AlphaSpend' must be a single value, not a vector.",call. =FALSE)}

if(AlphaSpend<0&AlphaSpend!="n"){stop("'AlphaSpend' must be a positive number smaller than or equal to 'alpha'.",call.=FALSE)}

if(AlphaSpend>alpha&AlphaSpend!="n"){stop(c("'AlphaSpend' must be smaller than or equal to ",alpha,"."),call. =FALSE)}






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



if( sum(is.numeric(test))!=1|length(test)>1){stop("Symbols and texts are not applicable for 'test'. It must be an integer greater than zero.",call. =FALSE)}
if(test<0){stop("'test' must be an integer greater than zero.",call. =FALSE)}



####
## Uploading information from previous tests
####

inputSetUp<- read.table(name)
if(inputSetUp[1,1]!=test-1){
aux_aux<- 1
message(c("The current test should be"," ",inputSetUp[1,1]+1,". ", "If you do not have information about previous tests, see the user manual for more details. If you are saving this run in an object, then the old data and analyses are uploaded on it."),domain = NULL, appendLF = TRUE)
                           }else{aux_aux<- 0} 


if(inputSetUp[1,1]>0&aux_aux==1){nameb<- paste(name1,"results.txt",sep=""); result2<- readRDS(nameb)}








#####
#####  OPEN IMPORTANT GLOBAL TEST FOR THE CASE WHEN THE WRONG test INPUT IS ENTERED. THUS, THIS FUNCTION ONLY RETURNS THE TABLE WITH INFORMATION FROM PREVIOUS TESTS
#####

if(aux_aux==0){


## inputSetUp matrix contains:
# line 1: (C11) the index for the order of the test (zero entry if we did not have a first test yet), (C12) SampleSize, (C13) alpha, (C14) mref (number of Monte Carlo simulations to calculate p-value), (C16) title, (C17) reject (the index indicating if and when H0 was rejected), (C18) pho
# line 2: says if the analysis has been started or not. 0 for not started. 1 for already started.
# line 3: critical values in the scale of the test statistic U
# line 4: events  
# line 5: actual alpha spent
# line 6: testv_old => identifier of look (test) per observed event
# line 7: has the target alpha spending until the (test-1)th look.
# line 8: 
# line 9: CVl
# line 10: mu_old => mus of each original entry per test 
# line 11: matched case-control <= ?
# line 12: Target alpha spending defined with AnalyzeSetUpRegression.Binomial
# line 13: MCpvalues per test
# line 14: number of entries (nent) per test
# line 15: the first column has R0, and 
# line 16: the first column has cases_fraction <= ?
# line 17: lower limit of the confidence interval per test
# line 18: upper limit of the confidence interval per test
# line 19--(number of covariates +1): matrix (y) with response variable (first collumn) and covariates 


SampleSize<- inputSetUp[1,2] ; N<- SampleSize
alpha<- inputSetUp[1,3]
rho<- inputSetUp[1,8]
r0<- inputSetUp[15,1] ; R0<- r0
nent<- length(events)
mref<- inputSetUp[1,4]
if(test==1){alphahs_old<- 0}else{alphahs_old<-  inputSetUp[7,1:(test-1)]}


if(test>1){   
   nent_old<-  inputSetUp[14,1:(test-1)]   
   testv_old<- inputSetUp[6,1:sum(inputSetUp[6,]>0)]
   mu0_old<- inputSetUp[10,1:sum(inputSetUp[10,]>0)]
   events_old<- inputSetUp[4,1:sum(inputSetUp[4,]>0)]
   yold<- inputSetUp[19:(19+ncol(covariates)), 1:sum(events_old)]
   covariates_old<- inputSetUp[(19+ncol(covariates)+1):((19+2*ncol(covariates))) , 1:sum(nent_old)]
   MCpvalues_old<- inputSetUp[13,1:(test-1)] 
          }





# Setting the names of the covariates
if(sum(colnames(covariates)==NULL)==0){colnames(covariates)<- paste(rep("X",ncol(covariates)),seq(1,ncol(covariates)),sep=""); covarnames<- paste(rep("X",ncol(covariates)),seq(1,ncol(covariates)),sep="")}






# Organizing the information per event according to Section 3.2 of Silva et al(2025)
ns<- events 

for(i in 1:length(events)){
 
  if(i==1){
           y<- matrix(0,events[i],ncol(covariates)+1); for(j in 1:events[i]){y[j,]<- c(mu0[i]/events[i],covariates[i,])}
          }

  if(i>1){
          yaux<- matrix(0,events[i],ncol(covariates)+1); for(j in 1:events[i]){yaux[j,]<- c(mu0[i]/events[i],covariates[i,])}; y<- rbind(y,yaux)
         }
                           } 












###########################################################
## HERE STARTS THE CODE FOR THE REGRESSION MODEL WITH A FIXED CHUNK OF POISSON DATA
###########################################################      


RegPoisson<- function(yt,alpha)
  
{

## alpha: to be used in the construction of the confidence intervals for the regression coefficients
     
if(ncol(yt)>1){k<- ncol(yt)-1; if(k>1){x<- yt[,2:ncol(yt)]}else{x<- matrix(yt[,2:ncol(yt)],ncol=k)}} ## k is the number of covariates
  
  mus<- yt[,1]
    







 ####### function for relative risk (rr) given theta= c(delta0,delta1,...,deltak) and data    
      
      
      rr<- function(theta){
        
        rrs<- rep(0,length(mus))
        
        for(tt in 1:length(mus)){
          
          if(ncol(yt)>1){rrs[tt]<- exp(sum(c(1,as.numeric(x[tt,1:k]))*theta[1:(1+k)]))}else{rrs[tt]<- exp(theta[1])}
          
        }
        
        return(rrs)
        
      }








      




    
    ################ FUNCTION FOR ESTIMATING THE PARAMETERS ###############
    
   
      
      
            ################ AUXILIARY FUNCTIONS ################################
      
      ######################################################################################
      
         
      
      ####### log-likelihood given rrs, theta and data
      
      
      
      
      LLR<- function(rrs){  
        
      return( sum( log(rrs)-rrs*mus ) )
        
                         }
      
      
      
      
   
      
      
      ####### First derivative of the log-likelihood with respect to delta_j,  j=0, 1, ..., k:
      
      
      
      
      d_deltaj<- function(j,rrs,theta){
        
        if(j==0){             
          return( sum( (1-rrs*mus) ) )
                }else{
                      return( sum( x[,j]*(1-rrs*mus) ) )
                     }
        
                                      }
      
     
      
     
      



      
      ####### Second derivative of the log-likelihood with respect to delta_j and delta_u
      
      
      
      d_deltaj_deltau<- function(j,u,rrs,theta){
               
                if(j==0){xxj<- rep(1,nrow(yt))}else{xxj<- x[,j]}
        
                if(u==0){xxu<- rep(1,nrow(yt))}else{xxu<- x[,u]}

                return( -sum( xxj*xxu*rrs*mus ) ) 
        
                                               }
      
            
      
      
      



      
      ######## Vector with first derivatives for given theta and rrs
      
       F<- function(rrs,theta){
        
        FF<- matrix(0,nrow(theta),1)
        
        FF[1,1]<- d_deltaj(0,rrs,theta)
        
        if(k>0){for(l in 2:(k+1)){FF[l,1]<- d_deltaj(l-1,rrs,theta)}}
        
        return(FF) 
        
      }
      
      
      
      
      


     ######## Matrice with second derivatives for given theta
      
      J<- function(rrs,theta){
        
        JJ<- matrix(0,nrow(theta),nrow(theta))
        
        JJ[1,1]<- d_deltaj_deltau(0,0,rrs,theta)   
        
        if(k>0){
          
          for(l in 2:(k+1)){JJ[l,l]<- d_deltaj_deltau(l-1,l-1,rrs,theta)}
          
          
          
          for(l in 1:k){
            
            for(cc in (l+1):(k+1)){
              
              secondDeriv<- d_deltaj_deltau(cc-1,l-1,rrs,theta)
              
              JJ[l,cc]<- secondDeriv; JJ[cc,l]<- secondDeriv 
              
            }
            
          }
          
        }   
        
        
        
        return(JJ)
        
        
        
      }## close J
      
            
      




      
      ####### ESTIMATING THE PARAMETERS THROUGH Newton-Raphson
      
      
      
       
      thetaN<- matrix(rep(0,k+1),ncol=1); thetaN[1]<- 1
      
      
      
      theta_old<- matrix(0,nrow(thetaN),1)
      
      
      
      
      count<- 0
      Fh<- 1
      
      while(max(abs(Fh))>0.0001&count<100){
        
        count<- count +1
        
        rrs<- rr(thetaN)
        
        Fh<- F(rrs,thetaN)
        
        Jh<- J(rrs,thetaN)
        
        theta_old<- thetaN
        
        thetaN<- thetaN - solve(Jh)%*%Fh 
        
      }
      
      
      rrs<- rr(thetaN)
      llr<- LLR(rrs)
      rrs0<- rrs; rrs0[rrs0>r0]<- r0 
      llr_stat<- LLR(rrs) - LLR(rrs0)  # likelihood ratio test statistic
      

      #Jh is the matrix needed to get the observed information matrix
      Jh<- J(rrs,thetaN)
      
          

delta_matrix= Jh


#observed information matrix
    information_obs=-delta_matrix
#variance-covariance matrix 
    sigma_matrix=solve(information_obs)
#lower and upper bounds of the intervals
    lowerCI=thetaN-sqrt(diag(sigma_matrix))*qnorm(1-alpha/2)
    upperCI=thetaN+sqrt(diag(sigma_matrix))*qnorm(1-alpha/2)
    Confidence_Interval<- matrix(c(lowerCI,upperCI),ncol=2)
        
     return(list(thetaN,llr,Fh,Jh,rrs,Confidence_Interval,count,llr_stat,sigma_matrix))

} ## CLOSE RegPoisson


###########################################################
## HERE CLOSES THE CODE FOR THE REGRESSION MODEL WITH A FIXED CHUNK OF DATA
###########################################################










# Joint of old and current data

if(test==1){
testv<- rep(test,nrow(y))
           }else{
testv<- c(as.numeric(testv_old), rep(test,nrow(y)))
y<- rbind(t(yold),y)
covariates<- rbind(t(covariates_old),covariates); colnames(covariates)<-  covarnames
events<- c(as.numeric(events_old),events)
mu0<- c(as.numeric(mu0_old),mu0)
                }










###################################
################################### FUNCTION FOR THE TEST STATISTIC FOR EACH CHUNK OF DATA (test)
###################################


test_stat<- function(y)
{
for(i in 1:max(testv) ){
   if(ncol(y)>1){ yh<- y[testv<=i,]}else{yh<- matrix(y[testv<=i,],ncol=1)}
    res<- RegPoisson(yh,alpha)
    if(i==1){u<- max(res[[8]]); intervals<- matrix(res[[6]],ncol(y),2); coefs<- matrix(res[[1]],1,ncol(y))}else{u<- c(u,max(res[[8]])); intervals<- rbind(intervals,res[[6]]); coefs<- rbind(coefs,matrix(res[[1]],1,ncol(y)))}
                      }
    sigma_matrix<- res[[9]]
return(list(u,intervals,coefs,sigma_matrix))
}






###################################
################################### VECTOR WITH OBSERVED TEST STATISTICS FOR EACH CHUNK OF DATA (test)
###################################

res<- test_stat(y)
u0<- res[[1]]
ConfInter<- res[[2]]
coefs<- res[[3]]
sigma_matrix<- res[[4]]








###################################
################################### CALCULATING THE MONTE CARLO P-VALUES
###################################



# calculating the alphas and minimum Monte Carlo replications per test

for(i in 1:max(testv)){

# Setting the test specific alpha spending 
if(i==1){n<- sum(testv<=i); alphahs<- alpha*(n/N)^rho; m<- max(2*ceiling(1/alphahs)-1,mref)}else{
         n1<- sum(testv<=i-1); n2<- sum(testv<=i); alphahs<- c(alphahs,alpha*(n2/N)^rho-alpha*(n1/N)^rho)
         m<- max(m,2*ceiling(1/alphahs)-1)
                                                                                                }
         
                      }

if(is.numeric(AlphaSpend)==T){ alphahs[test]<- max(0,AlphaSpend - sum(alphahs_old) ) } 
if(test>1){alphahs[1:(test-1)]<- as.numeric(alphahs_old)}








# Monte Carlo replications

gera_mus<- function(mus){return(rexp(1,mus*r0))}

G<- rep(0,max(test))
i<- 1
while(i<=m){     
     w<- apply(matrix(y[,1],ncol=1),1,gera_mus)
     yy<- y; yy[,1]<- w  
     
testehh<- tryCatch({
res<- test_stat(yy) 
umc<- res[[1]]  # vector with the Monte Carlo test statistics
     G<- G+ (umc>=u0)*1
     i<- i+1
}, error=function(e){})

           }





MCpvalues<- (G+1)/(m+1) ;  if(test>1){MCpvalues<- c(as.numeric(MCpvalues_old),MCpvalues[test])}


names(MCpvalues)<- paste(rep("test",test),rep("_",test),seq(1,test),sep="")
test_statistic_per_test<- u0; names(test_statistic_per_test)<- paste(rep("test",test),rep("_",test),seq(1,test),sep="")
AlphaSpending_per_test<- alphahs; names(AlphaSpending_per_test)<- paste(rep("test",test),rep("_",test),seq(1,test),sep="")
Coefficients_per_test<- coefs; rownames(Coefficients_per_test)<- paste(rep("test",test),rep("_",test),seq(1,test),sep="")
colnames(Coefficients_per_test)[2:ncol(y)]<- colnames(covariates)
colnames(Coefficients_per_test)[1]<- "Intercept"  
Confidence_Intervals_per_test<- ConfInter; colnames(Confidence_Intervals_per_test)<- c("Lower_limit","Upper_limit")
rownames(Confidence_Intervals_per_test)<- paste(rep("test",test*(ncol(covariates)+1)),rep("_",test*(ncol(covariates)+1)),rep(seq(1,test),rep(ncol(covariates)+1,test)), rep("_",test*(ncol(covariates)+1)),rep(colnames(Coefficients_per_test),test),sep="")

cum_alpha<- matrix(as.numeric(alphahs),1,)%*% upper.tri(matrix(0,test,test),diag=T)

Results<- data.frame(matrix(0,test,5))
Results[,5]<- rep("No",length(AlphaSpending_per_test))
if(sum(MCpvalues<=AlphaSpending_per_test)>0){Results[min(seq(1,length(MCpvalues))[MCpvalues<=AlphaSpending_per_test]):length(MCpvalues),5]<- "Yes"}
Results[,1]<- round(test_statistic_per_test,6) ; Results[,4]<- round(AlphaSpending_per_test,6)
Results[,3]<- round(MCpvalues,6) ; Results[,2]<- Results[,2]<- round(as.numeric(cum_alpha),6)
rownames(Results)<- paste(rep("test",test),rep("_",test),seq(1,test),sep="")
colnames(Results)<- c("Test statistic","Alpha spending","p-value","Test specific alpha spending","H0 rejected?")


if(test==1){nent_glob<- nent}else{nent_glob<- c(nent_old,nent)} 
for(i in 1:length(nent_glob)){
if(i==1){nomesRR<- paste(rep("test",nent_glob[i]),rep("_",nent_glob[i]),rep(i,nent_glob[i]),sep="")}else{
  nomesRR<- c(nomesRR,paste(rep("test",nent_glob[i]),rep("_",nent_glob[i]),rep(i,nent_glob[i]),sep=""))
                                                                                                        }
                             }



Risco_relativo<- round( exp(cbind( matrix( rep(1,nrow(covariates)), ncol=1) ,  covariates  ) %*% Coefficients_per_test[test,] ), 4)
A<- cbind( matrix( rep(1,nrow(covariates)), ncol=1) ,  covariates  )
B<- A%*%sigma_matrix%*%t(A)
lower<- round( Risco_relativo - sqrt(diag(B))*qnorm(1-alpha/2) , 4); for(i in 1:length(lower)){lower[i]<- max(0,lower[i])}
upper<- round( Risco_relativo + sqrt(diag(B))*qnorm(1-alpha/2), 4)
Relative_risk_per_data_entry<- cbind( matrix(events,ncol=1), t(t(matrix(events,ncol=1)) %*% upper.tri(matrix(0,length(events),length(events)),diag=T)), covariates, matrix(Risco_relativo,ncol=1), matrix(lower,ncol=1), matrix(upper,ncol=1) )
rownames(Relative_risk_per_data_entry)<- nomesRR
colnames(Relative_risk_per_data_entry)[1]<- "Events"
colnames(Relative_risk_per_data_entry)[2]<- "Cum. events"
colnames(Relative_risk_per_data_entry)[ncol(Relative_risk_per_data_entry)-2]<- c("Relative risk estimates")
colnames(Relative_risk_per_data_entry)[ncol(Relative_risk_per_data_entry)-1]<- c("Relative risk lower limit")
colnames(Relative_risk_per_data_entry)[ncol(Relative_risk_per_data_entry)]<- c("Relative risk upper limit")

 



message(c("                                ",title),domain = NULL, appendLF = TRUE)
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
          message(paste(c("=>    H0 is rejected after (p-value) <= (Test specific alpha spending) for the first time.")),domain = NULL, appendLF = TRUE)  
message("-------------------------------------------------------------------------------------------",domain = NULL, appendLF = TRUE)
options("width"=300)
result2<- list(Results,Relative_risk_per_data_entry,Coefficients_per_test,Confidence_Intervals_per_test)
names(result2)<- c("Decision_table","Relative_risk_estimates","Coefficients","Confidence_Intervals_for_coefficients")
print(result2,right=TRUE,row.names=FALSE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)
message(c("Parameter settings: N= ",SampleSize,", alpha= ",alpha,", rho= ",rho,", and ", "H0: RR<=",R0, "."),domain = NULL, appendLF = TRUE)
message(c("Analysis performed on ",date(),"."),domain = NULL, appendLF = TRUE)
message("===========================================================================================",domain = NULL, appendLF = TRUE)



#plot(seq(1,length(MCpvalues)),MCpvalues,type="l",col="blue",xlab="Test",ylab="P-value",ylim=c(0,0.05))
#lines(seq(1,length(MCpvalues)) , AlphaSpending_per_test , col="red")




############################################################
## SAVING INFORMATION FOR FUTURE TESTS
############################################################

inputSetUp[1,1]<- test

inputSetUp[6,1:nrow(y)]<- testv

inputSetUp[7,test]<- alphahs[test] 

inputSetUp[10,1:length(mu0)]<- mu0

if(nrow(inputSetUp)<19+2*ncol(covariates)){inputSetUp<- rbind(inputSetUp, matrix(0,19+2*ncol(covariates)-nrow(inputSetUp)+1,ncol(inputSetUp)))}
if(ncol(inputSetUp)<nrow(y)){inputSetUp<- cbind(inputSetUp, matrix(0,nrow(inputSetUp),nrow(y)-ncol(inputSetUp)+1))}
inputSetUp[19:(19+ncol(covariates)), 1:nrow(y)] <- t(y)
inputSetUp[(19+ncol(covariates)+1):(19+2*ncol(covariates)), 1:nrow(covariates)] <- t(covariates)

inputSetUp[13,test]<- MCpvalues[test] 

inputSetUp[14,test]<- nent

if(test>1){
inputSetUp[4,1:(sum(nent_old)+nent)]<- events
          }else{
inputSetUp[4,1:nent]<- events
               }




write.table(inputSetUp,name)

saveRDS(result2,file=paste(name1,"results.txt",sep=""))


#####
#####  CLOSES IMPORTANT GLOBAL TEST
#####
                   } 


invisible(result2)

###############################---------------------------
}## CLOSE THE AnalyzeRegression.Poisson FUNCTION 
##########################################################






