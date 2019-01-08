
AnalyzeSetUp.wBinomial<- function(name,N,alpha=0.05,M=1,rho=0.5,title="n",address="n",Tailed=1)
{

# Example of address: "C:/Users/Ivair/Documents"

if(address=="n"){stop(c("Please, provide a valid directory address to save the setup information."),call. =FALSE)}

if(is.numeric(rho)!=TRUE){stop("Symbols and texts are not applicable for 'rho'. It must be a positive number.",call. =FALSE)}
if(rho<=0){stop("rho must be greater than zero or equal to the default (rho='n')",call. =FALSE)}

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
then try 'AnalyzeSetUp.wBinomial' again."),call. =FALSE)
                           }
MinCases<- M

if( sum(is.numeric(alpha))!=1){stop("Symbols and texts are not applicable for 'alpha'. It must be a number in the (0,0.5) interval.",call. =FALSE)}
if( sum(is.numeric(MinCases))!=1){stop("Symbols and texts are not applicable for 'M'. It must be an integer greater than zero.",call. =FALSE)}
if(is.numeric(N)==FALSE){stop("Symbols and texts are not applicable for 'N'. It must be an integer greater than zero.",call. =FALSE)}


if(N<=0){stop("'N' must be an integer greater than zero.",call. =FALSE)}

if(alpha<=0|alpha>0.5||is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(MinCases>N||is.numeric(MinCases)==FALSE){stop("'M' must be an integer smaller than or equal to 'N'.",call. =FALSE)}
if(MinCases<1){stop("'M' must be an integer greater than zero.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("'M' must be an integer.",call. =FALSE)}

#############################################################################################################
##   HERE WE SAVE THE KEY CONTENT TO SETUP THE SURVEILLANCE. THE CONTENT IS SAVED IN THE MATRIX  CALLED inputSetUp 
#############################################################################################################
## inputSetUp matrix contains:
# line 1: (C11) the index for the order of the test (zero entry if we did not have a first test yet), (C12) SampleSize, (C13) alpha, (C14) M, (C15) base(the line of p where the looping will start in the next test), (C16) title, (C17) reject (the index indicating if and when H0 was rejected), (C18) pho (zero if Wald is used), (C19) Tailed (1 for one-sided test and 2 for two-sided)
# line 2: says if the analysis has been started or not. 0 for not started. 1 for already started.
# line 3: critical values in the scale of Sn
# line 4: observed cases  
# line 5: actual alpha spent
# line 6: cumulative expected value of Sn E[sum(w*C)] under H0, test by test
# line 7: has the target alpha spending until the (test-1)th look.
# line 8: observed controls
# line 9: CVl
# line 10: the number of rows (different weights/z ) per test
# line 11: matched case-control
# line 12: PENSAR
# line 13: weight for each observation
# line 14: test per weight

inputSetUp<- as.data.frame(matrix(0,14,9))

inputSetUp[1,]<- 0
inputSetUp[1,1:9]<- c(0,N,alpha,M,1,0,0,rho,Tailed)
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

write.table(inputSetUp,name)
titlecheck<- data.frame(matrix(0,1,1))
if(title!=0){titlecheck[1,1]<- title}
 
message(c("The parameters were successfully set at '",address,"'."),domain = NULL, appendLF = TRUE)
message(c("The temporary directory of your computer has the address of the directory where the settings information of this sequential analysis is saved.
Thus, do not clean the temporary directory before finishing this sequential analysis."),domain = NULL, appendLF = TRUE)

write.table(titlecheck,paste(name1,"title.txt",sep=""))

setwd(safedir)

} ## end function AnalyzeSetUp.Binomial








