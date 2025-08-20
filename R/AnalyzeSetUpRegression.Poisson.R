
# -------------------------------------------------------------------------
# Function to perform the unpredictable Poisson MaxSPRT regression surveillance - Version edited at July-29-2025
# -------------------------------------------------------------------------

AnalyzeSetUpRegression.Poisson<- function(name,N="n",alpha=0.05,R0=1,rho=1,mref=999,title="n",address="n")
{


# Example of address: "C:/Users/Ivair/Documents"

if(address=="n"){stop(c("Please, provide a valid directory address to save the setup information."),call. =FALSE)}





pho<- rho


phoref<- rho

if(pho<=0){stop("rho must be greater than zero or equal to the default (rho=1).",call. =FALSE)}





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
then try 'AnalyzeSetUpRegression.Binomial' again."),call. =FALSE)
                           }






if( sum(is.numeric(alpha))!=1){stop("Symbols and texts are not applicable for 'alpha'. It must be a number in the (0,0.5) interval.",call. =FALSE)}
if(is.numeric(N)==FALSE&N!="n"){stop("Symbols and texts are not applicable for 'N'. It must be an integer greater than zero or equal to the default 'n'.",call. =FALSE)}


if(N<=0){stop("'N' must be an integer greater than zero or equal to the default 'n'.",call. =FALSE)}

if(alpha<=0|alpha>0.5||is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}






alpha1<- alpha
#posi<- length(total_cases)+1
posi<- 2

rejt<- 0








r0<- R0





#############################################################################################################
##   HERE WE SAVE THE KEY CONTENT TO SETUP THE SURVEILLANCE. THE CONTENT IS SAVED IN THE MATRIX  CALLED inputSetUp 
#############################################################################################################
## inputSetUp matrix contains:
# line 1: (C11) the index for the order of the test (zero entry if we did not have a first test yet), (C12) SampleSize, (C13) alpha, (C14) mref (number of Monte Carlo simulations to calculate p-value), (C16) title, (C17) reject (the index indicating if and when H0 was rejected), (C18) pho
# line 2: says if the analysis has been started or not. 0 for not started. 1 for already started.
# line 3: critical values in the scale of the test statistic U
# line 4:   
# line 5: actual alpha spent
# line 6: testv_old => identifier of look (test) per observed event
# line 7: has the target alpha spending until the (test-1)th look.
# line 8: none
# line 9: CVl
# line 10: mu_old => mus's per event
# line 11: matched case-control <= ?
# line 12: Target alpha spending defined with AnalyzeSetUpRegression.Binomial
# line 13: MCpvalues per test
# line 14: 
# line 15: the first column has R0, and 
# line 16: the first column has cases_fraction <= ?
# line 17: lower limit of the confidence interval per test
# line 18: upper limit of the confidence interval per test
# line 19--(number of covariates +1): matrix (y) with response variable (first collumn) and covariates 


inputSetUp<- as.data.frame(matrix(0,30,max(N,10)))

inputSetUp[1,]<- 0
inputSetUp[1,1:10]<- c(0,N,alpha,mref,1,0,0,rho,0,0) 
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
inputSetUp[13,]<- 0
inputSetUp[14,]<- 0
inputSetUp[15,1]<- R0
 




write.table(inputSetUp,name)
titlecheck<- data.frame(matrix(0,1,1))
if(title!=0){titlecheck[1,1]<- title}
 
message(c("The parameters were successfully set at '",address,"'."),domain = NULL, appendLF = TRUE)
message(c("The temporary directory of your computer has the address of the directory where the settings information of this sequential analysis is saved.
Thus, do not clean the temporary directory before finishing this sequential analysis."),domain = NULL, appendLF = TRUE)


write.table(titlecheck,paste(name1,"title.txt",sep=""))

setwd(safedir)



}
########### CLOSES THE AnalyzeSetUpRegression.Poisson
######################################################





