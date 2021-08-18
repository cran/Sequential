
# -------------------------------------------------------------------------
# Function to perform the unpredictable binomial MaxSPRT surveillance - Version edited at August-2021
# -------------------------------------------------------------------------

AnalyzeSetUp.Poisson<- function(name,SampleSize,alpha=0.05,D=0,M=1,AlphaSpendType="Wald",rho="n",title="n",address="n",Tailed="upper")
{

M1<- M
if(Tailed!="upper"){stop("For this version of the Sequential package, CV.Binomial works only for 'Tailed=upper'.",call. =FALSE)}

# Example of address: "C:/Users/Ivair/Documents"

if(address=="n"){stop(c("Please, provide a valid directory address to save the setup information."),call. =FALSE)}
                                    

pho<- rho

phoref<- rho
T<- SampleSize

if(AlphaSpendType!="Wald"&AlphaSpendType!="power-type"){stop("Set AlphaSpendType= 'Wald' or AlphaSpendType= 'power-type'.",call. =FALSE)}
if(AlphaSpendType=="power-type"&is.numeric(pho)!=TRUE){stop("Symbols and texts are not applicable for 'rho'. It must be a positive number.",call. =FALSE)}
if(pho<=0&AlphaSpendType=="power-type"){stop("rho must be greater than zero or equal to the default (rho='n')",call. =FALSE)}

if(pho=="n"){pho<- 0}

if(AlphaSpendType=="Wald"){pho<- 0}

safedir<- getwd()
if(title== "n"){title<- 0}
name1<- name

# address<- choose.dir(default = "", caption = paste("Select the folder where the file '",name,"' is going to be saved.")) ## Old version of choose dir



address1<- tempdir()
y<- substr(address1,1:4,4) ; i<- 2
while(i<=nchar(address1)-3&y!="Temp"&y!="TEMP"){y<- substr(address1,i:(i+3),i+3);i<- i+1}
address1<- substr(address1,1:(i+3),i+3)

if(paste(address)=="temp"){address<- address1}
address2<- data.frame(c(0))
address2[1,1]<- address
setwd(address1)
write.table(address2,paste(name,"address.txt",sep=""),sep=";")
setwd(address)
name<- paste(name,".","txt",sep="")
if(file.exists(name)==TRUE){
stop(c("There already exists a file called"," ",name1,".
","You may want check if some test has been performed for this monitoring before. 
If you really want to overwrite the existent file, please, go to ",address," 
to delete the file '",name,"'. Alternatively, you can delete that file by using the 
following commands: ", "setwd(","'",getwd(),"'",")","; ", "file.remove(","'",name,"'","), 
then try 'AnalyzeSetUp.Binomial' again."),call. =FALSE)
                           }
MinCases<- M

if( sum(is.numeric(alpha))!=1){stop("Symbols and texts are not applicable for 'alpha'. It must be a number in the (0,0.5) interval.",call. =FALSE)}
if( sum(is.numeric(MinCases))!=1){stop("Symbols and texts are not applicable for 'M'. It must be an integer greater than zero.",call. =FALSE)}
if( sum(is.numeric(D))!=1){stop("Symbols and texts are not applicable for 'D'. It must be a number greater than or equal to zero.",call. =FALSE)}

if(is.numeric(SampleSize)==FALSE){stop("Symbols and texts are not applicable for 'SampleSize'. It must be an integer greater than zero.",call. =FALSE)}

if(SampleSize<=0){stop("'SampleSize' must be an integer greater than zero.",call. =FALSE)}

if(alpha<=0|alpha>0.5||is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(MinCases>SampleSize||is.numeric(MinCases)==FALSE){stop("'M' must be an integer smaller than or equal to 'SampleSize'.",call. =FALSE)}
if(MinCases<1){stop("'M' must be an integer greater than zero.",call. =FALSE)}
if(MinCases!=round(MinCases)){stop("'M' must be an integer.",call. =FALSE)}

if(D>SampleSize||is.numeric(D)==FALSE){stop("'D' must be a number smaller than or equal to 'SampleSize'.",call. =FALSE)}
if(D<0){stop("'D' must be a number greater than or equal to zero.",call. =FALSE)}



alpha1<- alpha

posi<- 2

rejt<- 0

# -------------------------------------------------------------------------
# Function produces alpha spending associated to flat critical values - continuous Poisson MaxSPRT
# -------------------------------------------------------------------------

SalphafLAtcv <- function(SampleSize,alpha,MinCases) {

#---------------------------------------------------------------------
# Function that calculates the product log through a recursive formula
#---------------------------------------------------------------------
ProdLog <- function(z){
	x = z-z^2+1.5*z^3-(8/3)*z^4+(125/24)*z^5-(54/5)*z^6+(16807/720)*z^7
	for(i in 1:10) x = x-(x*exp(x)-z)/(exp(x)+x*exp(x))
	x
	                } # end ProdLog function 

#----------------------------------------------------------------------------------------------
# Function that calculates the probability of type I error for a given set of IMPUT parameters
#----------------------------------------------------------------------------------------------
Perror_I<- function(cv){
 
z = -exp(-1-cv/c)
mu = -c * ProdLog(z) 		#The expected counts under H0 that is needed to reject the null with i number of adverse events
mtemp = c(0,mu)
mmu = diff(mtemp) 		#The marginal difference of the mu[] vector
imin=MinCases


while (mu[imin] < D) imin=imin+1
if(imin>MinCases) { 
	mu[imin-1]=D
	mmu[imin]=mu[imin]-D
	} # end if 


imax=1
while (mu[imax] < T) imax=imax+1    		# imax is the maximum number of cases that will generate a signal.            

# NOTE: If imax=1, this code will not work


if(imin<imax){

# Defining the p[][] matrix
# -------------------------

p = seq(length=(imax-1)*imax, from=0, by=0)				# p[i,j] is the probability of having j-1 cases at time mu[i]
dim(p) = c(imax-1,imax)								# i in 1:imax-1 is the rows and j in 1:imax is the column

# Calculating the first row p[imin][] in the matrix for which there is a chance to reject H0
# When MinCases=1, there is no skipping, and it is the first row in the matrix (p[1][]).
# ------------------------------------------------------------------------------------------

if(imin==MinCases) {
	for(s in 1:imin) p[imin,s] = dpois(s-1,mu[imin])			# Probability of having s-1 cases at time mu[imin], not rejectinh H0
	p[imin,imin+1] = 1-ppois(imin-1,mu[imin])				# Probability of having s+ cases at time mu[imin], rejectinh H0
	} # end if

if(imin>MinCases) {
	for(s in 1:imin) p[imin-1,s]=dpois(s-1,mu[imin-1])		# Probability of having s-1 cases at time mu[imin-1], not rejecting H0
	p[imin-1,imin+1] = 1-ppois(imin-1,mu[imin-1])				# Probability of having s+ cases at time mu[imin-1], rejecting H0
	for(s in 1:imin) 								# Probability of having s-1 cases at time mu[imin], not rejectinh H0
		for(k in 1:s) 
			p[imin,s]=p[imin,s]+p[imin-1,k]*dpois(s-k,mmu[imin])	
	for(k in 1:imin) 
		p[imin,imin+1] = p[imin,imin+1] + p[imin-1,k]*(1-ppois(imin-k,mmu[imin]))
} # end if 

funcaux1<- function(ii){j<- matrix(seq(1,(ii-1)),,1); ptes<- apply(j,1,funcaux2,ii); return(ptes)}
funcaux2<- function(jj,ii){k<- seq(1,jj); return(sum(p[ii-1,k]*dpois(jj-k,mmu[ii])) ) }
funcaux3<- function(ii){k<- seq(1,ii-1); return(sum(p[ii-1,k]*dpois(ii-k,mmu[ii])) ) }
funcaux4<- function(ii){k<- seq(1,ii-1); return(sum(p[ii-1,k]*(1-ppois(ii-k,mmu[ii])) ) ) }

# Calculating the remaining rows in the p[][] matix
# -------------------------------------------------

if(MinCases+1<=imax-1&((imin+1)<=(imax-1)))
probaux3<- 0
i<- (imin+1)
while(i <=(imax-1)&probaux3<=alpha+PRECISION) {

p[i,1:((i-1))]<- funcaux1(i)
p[i,i]<- funcaux3(i)
p[i,i+1]<- funcaux4(i)

probaux3<- probaux3 + p[i,i+1]
i<- i+1
} # end for i	


pp=0
if(imax>imin)
for(k in 1:(imax-1)) pp=pp+p[imax-1,k]*(1-ppois(imax-k,T-mu[imax-1])) #Calculates the last probability to signal before time SampleSize


# Sums up the probabilities of absorbing states when a signal occurs, to get the alpha level
# ------------------------------------------------------------------------------------------

alpha_I=0 ; Salpha<- rep(0,imax-1) 

if(imin>MinCases){ alpha_I=p[imin-1,imin+1] ; Salpha[imin-1]<- p[imin-1,imin+1]}
for(i in imin:(imax-1)){ alpha_I=alpha_I+p[i,i+1] ; Salpha[i]<- p[i,i+1] }					
alpha_I=alpha_I+pp ; Salpha[imax-1]<- Salpha[imax-1]+pp


}else{alpha_I<- 1-ppois(imax-1,mu[imax]); Salpha<- rep(0,imax-1) ; Salpha[imax-1]<-  1-ppois(imax-1,mu[imax])} # end if(imin<imax)

return(list(alpha_I,Salpha,mu[1:(imax-1)]))

                      } # end Perror_I
#############################################################################################

PRECISION<- 0.00000001
CV1<- 0 ; CV2<- 10 ; trunc<- ceiling(log(1/PRECISION)/log(2))

cc<- 0
mum<- 0
while(mum<T){
cc<- cc+1
zm = -exp(-1-CV2/cc)
mum = -cc * ProdLog(zm)
            }

c = 1:cc


cont<- 1
alphaobs<- 0

while(abs(alphaobs-alpha)>PRECISION&cont<=trunc){
        CVm<- (CV1+CV2)/2
        resE<- Perror_I(cv=CVm)
        alphaobs<- resE[[1]]
        if(alphaobs>alpha){CV1<- CVm}else{CV2<- CVm}
        cont<- cont+1         
                                                }
Salpha<- resE[[2]]
mut<- resE[[3]]
return(list(Salpha,mut))

} # CLOSES FUNCTION THAT OBTAINS ALPHA SPENDING FOR THE POISSON MAXSPRT
 




##############################################################################################################
## HERE THE TARGET ALPHA SPENDING IS DEFINED WHEN AlphaSpendTyp=Wald
##############################################################################################################

if(AlphaSpendType=="Wald"){
resE<- SalphafLAtcv(SampleSize,alpha=alpha1,MinCases)
sa<- resE[[1]]
mut<- resE[[2]]
j<- length(mut)
if(sum(sa)==0){stop("Choose larger SampleSize. It is not possible to find a solution for the desired alpha with the current SampleSize choice.",call. =FALSE)}
sum_sa<- sa%*%(upper.tri(matrix(0,length(sa),length(sa)),diag=T))
                          }else{j<- 0}

#############################################################################################################
##   HERE WE SAVE THE KEY CONTENT TO SETUP THE SURVEILLANCE. THE CONTENT IS SAVED IN THE MATRIX  CALLED inputSetUp 
#############################################################################################################
## inputSetUp matrix contains:
# line 1: (C11) the index for the order of the test (zero entry if we did not have a first test yet), (C12) SampleSize, (C13) alpha, (C14) M, (C15) base(the line of p where the looping will start in the next test), (C16) title, (C17) reject (the index indicating if and when H0 was rejected), (C18) pho (zero if Wald is used), (C19) j (the sample size in the scale of the events if rho=0, and j=0 otherwise), (C1,10) has D, (C1,11) has the M entered by the user 
# line 2: says if the analysis has been started or not. 0 for not started. 1 for already started.
# line 3: critical values in the scale of the events for each test
# line 4: observed events
# line 5: actual alpha spent
# line 6: expected number of events under H0, mu0, test by test
# line 7: has the target alpha spending actually used until the (test-1)th look.


if(AlphaSpendType=="Wald"){k<- length(sa)}
inputSetUp<- as.data.frame(matrix(0,7,11))

if(AlphaSpendType=="Wald"){alphaspend<- sum_sa; M<- sum(alphaspend==0)+1} #Target alpha spending to be spent event by event.

inputSetUp[1,]<- 0
inputSetUp[1,1:11]<- c(0,SampleSize,alpha,M,1,0,0,pho,j,D,M1) 
inputSetUp[2,]<- 0
inputSetUp[2,1]<- 0 # says if the surveillance was started or not.
inputSetUp[3,]<- 0
inputSetUp[4,]<- 0
inputSetUp[5,]<- 0
inputSetUp[6,]<- 0
inputSetUp[7,]<- 0


write.table(inputSetUp,name)
titlecheck<- data.frame(matrix(0,1,1))
if(title!=0){titlecheck[1,1]<- title}
 
message(c("The parameters were successfully set at '",address,"'."),domain = NULL, appendLF = TRUE)
message(c("The temporary directory of your computer has the address of the directory where the settings information of this sequential analysis is saved.
Thus, do not clean the temporary directory before finishing this sequential analysis."),domain = NULL, appendLF = TRUE)

if(AlphaSpendType=="Wald"&phoref!="n"){message(c("The value of 'rho' is ignored, as it is not used when AlphaSpendType='Wald'."),domain = NULL, appendLF = TRUE)}

write.table(titlecheck,paste(name1,"title.txt",sep=""))

if(AlphaSpendType=="Wald"){write.table(alphaspend, paste(name1,"alphaspend.txt",sep=""))}

setwd(safedir)

} ## end function AnalyzeSetUp.Poisson


# AnalyzeSetUp.Poisson(name="teste",SampleSize=6,alpha=0.05,D=2,M=1,AlphaSpendType="Wald",rho="n",title="teste",address="C:/Users/ivair/POST-DOC/Material para construcao do pacote Sequential/PASTA PARA TREINO",Tailed="upper")
