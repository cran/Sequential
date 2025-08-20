
# -------------------------------------------------------------------------
# Function to setup parameters for the Analyze.Multinomial function
# to perform the unpredictable multinomial marginal MaxSPRT surveillance - Version edited at May-25-2025
# -------------------------------------------------------------------------

AnalyzeSetUp.Multinomial<- function(name,N=200,alpha=0.05,AlphaSpendType=1,k,R0=1,R1=2,rho=1,pmin=0.05, pmax=0.95, target_power=0.8,Rmin=1,Rmax=10,gamma=0.9, m=100000,title="n",ExposuresNames="n",address="n")
{

# name: name to be used in each analysis to read the information saved from previus test.
# N: maximum length of surveillance.
# alpha: overall significance level.
# k: length of the multinomial vector.
# R0: test margin under the null hypothesis. Must be a positive number. Default is R0=1.  
# R1: relative risk under the alternative hypothesis given the target_power.
# rho: parameter to setup the shape of the power-type alpha spending.
# m: Monte Carlo replications of the multinomial for critical values calculations
# title: Optional. title of the table with results of analysis
# ExposuresNames: Optional. This is to inform the name of the exposures related to each entry of the multinomial vector. For example, it can be c("A","B","AB") indicating that the count entries in the "cases" vector are related to populations exposed to vaccines A, B, and AB, respectively.
# AlphaSpendType: the possible values are 1 and 2 according to Silva and Maro 2025.
# Rmin: minimum value for the relative risks of adjacent exposures when constructing the confidence interval for RR of a given exposure
# Rmax: maximum value for the relative risks of adjacent exposures when constructing the confidence interval for RR of a given exposure
# gamma: confidence coefficient of the confidence intervals for RR's 


#### Setting the address to save analyzes information

# Example of addresses: "C:/Users/Ivair/Documents"

if(address=="n"){stop(c("Please, provide a valid directory address to save the setup information."),call. =FALSE)}

safedir<- getwd()
if(title== "n"){title<- 0}
name1<- name

address1<- tempdir()
y<- substr(address1,1:4,4) ; i<- 2
while(i<=nchar(address1)-3&y!="Temp"&y!="TEMP"){y<- substr(address1,i:(i+3),i+3);i<- i+1}
address1<- substr(address1,1:(i+3),i+3)

address2<- data.frame(c(0))
address2[1,1]<- address
setwd(address1)
write.table(address2,paste(name,"address.txt",sep=""),sep=";")
setwd(address)
name<- paste(name,".","txt",sep="")
if(file.exists(name)==TRUE){
stop(c("There already exists a file called"," ",name1,".
","You may want check if some test has been performed for this monitoring before. 
If you really want to overwrite the existent file, please, delete that file by using the 
following commands: ", "setwd(","'",getwd(),"'",")","; ", "file.remove(","'",name,"'","), 
then try 'AnalyzeSetUp.Multinomial' again."),call. =FALSE)
                           }






#### Verifying the consistency of the input parameters
if(rho<=0){stop("rho must be greater than zero or equal to the default (rho=1).",call. =FALSE)}
if( sum(is.numeric(alpha))!=1){stop("Symbols and texts are not applicable for 'alpha'. It must be a number in the (0,0.5) interval.",call. =FALSE)}
if(is.numeric(N)==FALSE){stop("Symbols and texts are not applicable for 'N'. It must be an integer greater than zero or equal to the default 'n'.",call. =FALSE)}
if(N<=0){stop("'N' must be an integer greater than zero or equal to the default 'n'.",call. =FALSE)}
if(alpha<=0|alpha>0.5||is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}






#### Setting the power-type alpha spending function
x<- seq(1/N,by=1/N,1)
sum_sa<- alpha*(x^rho)






#### Information to be used in the analysis:


active<- rep(1,k) # entries of the multinomial that are still alive for future tests
if(length(ExposuresNames)!=k){ExposuresNames<- seq(1,k)}





#############################################################################################################
##   HERE WE SAVE THE KEY CONTENT TO SETUP THE SURVEILLANCE. THE CONTENT IS SAVED IN THE MATRIX  CALLED inputSetUp 
#############################################################################################################


inputSetUp<- as.data.frame(matrix(0, 19+6*k+6, max(N,10,k) ))


## inputSetUp matrix contains:
# line 1: (C11) the index for the order of the test (zero entry if we did not have a first test yet), (C12) Maximum SampleSize, (C13) alpha, (C14) k, (C15) m, (C16) title, (C17) reject (the index indicating if and when H0 was rejected), (C18) rho, (C19) vazio, (C1,10) vazio
# line 2: says if the analysis has been started or not. 0 for not started. 1 for already started.
# line 3: critical values in the scale of MaxSPRT
# line 4: observed cases  
# line 5: AlphaOld, actual alpha spent per test
# line 6: cumulative expected value of Sn E[sum(w*C)] under H0, test by test
# line 7: has the target alpha spending until the (test-1)th look.
# line 8: observed controls
# line 9: reject time
# line 10: the number of rows (different weights/z ) per test
# line 11: matched case-control
# line 12: Target alpha spending defined with AnalyzeSetUp
# line 13: weight for each observation
# line 14: test per weight
# line 15: test margin R0.
# line 16: AlphaSpendType
# line 17: active, vector informing which entries are still elegible for alpha spending in the current test. Entries must be settled equal to 1 for active(elegible), or equal to zero if the null hypothesis related to it has been rejected in some previous test.
# line 18: Ns, vector with test-specific number of events
# line 19: N_old, total cummulative number of events in previous tests


inputSetUp[1,]<- 0
inputSetUp[1,1:10]<- c(0,N,alpha,k,m,0,0,rho,0,0) 
inputSetUp[2,]<- 0
inputSetUp[2,1]<- 0
inputSetUp[3,]<- 0
inputSetUp[4,]<- 0
inputSetUp[5,]<- 0
inputSetUp[6,]<- 0
inputSetUp[7,]<- 0
inputSetUp[8,]<- 0
inputSetUp[9,]<- 0
inputSetUp[10,]<- 0
inputSetUp[11,]<- 0
inputSetUp[12,]<- 0
inputSetUp[12,1:length(sum_sa)]<- sum_sa
inputSetUp[13,]<- 0
inputSetUp[14,]<- 0
inputSetUp[15,1]<- R0
inputSetUp[16,1:k]<- 0; inputSetUp[16,1]<- AlphaSpendType; inputSetUp[16,2]<- pmin; inputSetUp[16,3]<- pmax; inputSetUp[16,4]<- target_power; inputSetUp[16,5]<- R1 ; inputSetUp[16,6]<- Rmin; inputSetUp[16,7]<- Rmax; inputSetUp[16,8]<- gamma  
inputSetUp[17,1:k]<- active
inputSetUp[18,]<- 0 
inputSetUp[19,1]<- 0




write.table(inputSetUp,name)
write.table(ExposuresNames,"ExposuresNames.txt")
titlecheck<- data.frame(matrix(0,1,1))
if(title!=0){titlecheck[1,1]<- title}
 
message(c("The parameters were successfully set at '",address,"'."),domain = NULL, appendLF = TRUE)
message(c("The temporary directory of your computer has the address of the directory where the settings information of this sequential analysis is saved.
Thus, do not clean the temporary directory before finishing this sequential analysis."),domain = NULL, appendLF = TRUE)

write.table(titlecheck,paste(name1,"title.txt",sep=""))

setwd(safedir)

} ## end function AnalyzeSetUp.Multinomial





