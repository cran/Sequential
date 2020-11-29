
###### This code finds optimal alpha spending to minimize expected lenght of surveillance restricted to expcted time to signal with binomial data

Optimal.Binomial<- function(Objective="ETimeToSignal",N="n",z="n",p="n",alpha,power,RR,GroupSizes="n",Tailed= "upper",ConfIntWidth=3,gama=0.9,R0=1)
{


# Objective: statistical performance measure to minimize. Options are "ETimeToSignal" and "ESampleSize". Default is "ETimeToSignal"
# N: Maximum length of surveillance
# alpha = desired alpha level
# z = matching ratio between exposed and unexposed cases  
# p = probability of having a case under the null hypothesis 
# RR: target relative risk. It will be a vector of length 2 in case of "Tailed=two", where RR[1] is a value in (0,1), and RR[2]>1.
# power: target power. It will be a vector of length 2 in case of "Tailed=two", where power[1] is the target power under RR[1], and power[2] under RR[2].
# GroupSizes:  Vector with the number of adverse events (exposed+unexposed) between two looks at the data, i.e, irregular group sizes. Important: Must sums up N. The default is continuous sequential testing.
# Tailed= "lower", "upper", "two".  Default is "upper".
# ConfIntWidth: a positive value indicating that the actual relative risk R will belong to the rang [Rhat-ConfIntWidth/2, Rhat+ConfIntWidth/2], where Rhat is the maximum likelihood estimator of R. Default is 3. 
# gama: confidence coefficient. Default is 0.9.
# R0: one number (>0) for the relative risk under H0, or a two-dimensional vector for H0: R0[1]<= R <= R0[2].

if(p=="n"&z=="n"){stop("Please, at least z or p must be provided.",call. =FALSE)}

if( z!="n"){if(sum(is.numeric(z))!=1){stop("Symbols and texts are not applicable for 'z'. It must be a number greater than zero.",call. =FALSE)}}
if(z<=0){stop("'z' must be a number greater than zero.",call. =FALSE)}

if(p!="n"){
if(is.numeric(p)!=TRUE){stop("Symbols and texts are not applicable for 'p'. It must be a probability measure.",call. =FALSE)}
if(z!="n"&p!="n"){if(p!= 1/(1+z)){stop("Both z and p are specified, but the required relationship that p=1/(1+z) does not hold. Please remove either the definition of z or the definition of p. Only one of them is needed.",call. =FALSE)}}
if(p<=0|p>=1){stop("p must be a number greater than zero and smaller than 1.",call. =FALSE)}
           }
if(p!="n"){z<- 1/p-1}

if(Tailed=="lower"){z<- 1/z}

if(is.numeric(RR)==FALSE){stop("'RR' must have numeric positive entries.",call. =FALSE)}

if(Tailed=="two"){if(length(RR)!=2){stop("'RR' must be a number greater than 1.",call. =FALSE)}}

if(sum(RR<=0)>0){stop("'RR' must be numeric and positive.",call. =FALSE)}
if(Tailed=="two"){
if(length(RR)!=2){stop("'RR' must be a two-dimensional vector with entries greater than zero.",call. =FALSE)}
if(length(power)!=2){stop("'power' must be a two-dimensional vector with entries in [alpha, 1).",call. =FALSE)}
                 }



if(alpha<=0|alpha>0.5|is.numeric(alpha)==FALSE){stop("'alpha' must be a number greater than zero and smaller than 0.5.",call. =FALSE)}
if(is.numeric(power)==FALSE){stop("'power' must be a number greater than 'alpha' and smaller than 1.",call. =FALSE)}
if(length(power)>2){stop("'power' must be a number, or a two-dimensional vector, with entries in [alpha, 1).",call. =FALSE)}
if(min(power)<alpha){stop("'power' must have numbers greater than 'alpha' and smaller than 1.",call. =FALSE)}
if(length(RR)>2){stop("RR must be a two-dimensional vector, with entries smaller than 1.",call. =FALSE)}

#### Finding minimum N in case of having it as one of the results of the function.

if(Tailed!="two"){if(length(RR)>1|length(power)>1){stop("For 'Tailed=lower' or 'Tailed=upper', the inputs 'power' and 'RR' must be single numbers, not vectors.",call. =FALSE)}}

if(min(power)<=alpha){stop("Please, choose a 'power' value greater than 'alpha'.",call. =FALSE)}

if(Tailed=="upper"){if(RR<1){stop("Please, for 'Tailed=upper' choose 'RR>1'.",call. =FALSE)}}
if(Tailed=="lower"){if(RR>1){stop("Please, for 'Tailed=lower' choose 'RR<1'.",call. =FALSE)}}

if(Tailed=="two"){
if(RR[1]>=1){stop("'RR[1]' must be smaller than 1.",call. =FALSE)}
if(RR[2]<=1){stop("'RR[2]' must be greater than 1.",call. =FALSE)}
                 }


if(Tailed=="two"){power<- power[2]; power_inf<- power[1]}

if(Tailed=="two"){
Nauxpow<- 1; poweaux<- 0; while(poweaux<power){Nauxpow<- Nauxpow+1; cst<- qbinom(1-alpha,Nauxpow,1/(1+z)) ; if(1-pbinom(cst-1,Nauxpow,1/(1+z))>alpha){cst<- cst+1}; poweaux<- 1-pbinom(cst-1,Nauxpow,1/(1+z/RR[2]))}
Nauxpow2<- 1; poweaux2<- 0; while(poweaux2<power_inf){Nauxpow2<- Nauxpow2+1; cst<- qbinom(alpha,Nauxpow2,1/(1+z)) ; if(pbinom(cst,Nauxpow2,1/(1+z))>alpha){cst<- cst-1}; poweaux2<- pbinom(cst,Nauxpow2,1/(1+z/RR[1]))}
if(N=="n"){N<- 20+ max(Nauxpow,Nauxpow2)}
                 }else{
Nauxpow<- 1; poweaux<- 0; while(poweaux<power){Nauxpow<- Nauxpow+1; cst<- qbinom(1-alpha,Nauxpow,1/(1+z)) ; if(1-pbinom(cst-1,Nauxpow,1/(1+z))>alpha){cst<- cst+1}; poweaux<- 1-pbinom(cst-1,Nauxpow,1/(1+z/RR))}
if(N=="n"){N<- 20+Nauxpow}
                      }

if(length(N)>1){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}
if(length(GroupSizes)>1&N=="n"){stop(" 'N' must be equal to the sum of the entries in 'GroupSizes'.",call. =FALSE)}


if(is.numeric(N)==FALSE){stop("The maximum length of surveillance, 'N', must be a positive integer or held equal to the default 'n'.",call. =FALSE)}
if(N!=round(N)|N<=0){stop("The maximum length of surveillance, 'N', must be a positive integer.",call. =FALSE)}



Groups<- GroupSizes
if(length(Groups)==1){
if(GroupSizes!="n"){
if(is.numeric(Groups)==FALSE){stop("'GroupSizes' must be an integer smaller than or equal to 'N'.",call. =FALSE)}
if(sum(Groups<=0)>0){stop("'GroupSizes' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'GroupSizes' must be a positive integer smaller than or equal to 'N'.",call. =FALSE)}

if(Groups==0){stop("'GroupSizes' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(N/Groups!=round(N/Groups)){stop("The maximum length of surveillance, 'N', must be a multiple of 'GroupSizes'.",call. =FALSE)}
if(Groups>N){stop("The maximum length of surveillance, 'N', must be a multiple of 'GroupSizes'.",call.=FALSE)}
                   }
}




if(length(Groups)>1){
if(sum(is.numeric(Groups))==0){stop("'GroupSizes' must be a vector of positive integers.",call. =FALSE)}else{
if(sum(Groups<=0)>0){stop("'GroupSizes' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(sum(Groups==round(Groups))!=length(Groups)){stop("'GroupSizes' must be a positive integer or a vector with positive integers.",call. =FALSE)}
if(sum(Groups)!=N){stop("'GroupSizes' must sum up equal to 'N'.",call. =FALSE)}
}
}





if(length(GroupSizes)<=20){Nref<- 300}else{Nref<- 240}

if(N>Nref){stop(c("Please, choose a N value smaller than or equal to ",Nref,"."),call. =FALSE)}

#### Checking if the desired power can be reached with the alpha and N provided

if(Tailed!="two"){## open 2
cst<- qbinom(1-alpha,N,1/(1+z)) ; if(1-pbinom(cst-1,N,1/(1+z))>alpha){cst<- cst+1}; powcheck<- 1-pbinom(cst-1,N,1/(1+z/max(RR)))
if(powcheck<power&N<Nauxpow){
if(Nauxpow<=Nref){   ###########<<<<<<<<<<<<<<<<<<<<<<<<<<<<<+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
stop(c("For the desired alpha, power, RR, and the specified z, the maximum length of surveillance, 'N', must be greater than or equal to ",Nauxpow),call. =FALSE)
                }else{
           stop(c("For the desired alpha, power, RR, and the specified z, the required 'N' is greater than ",Nauxpow, " which can lead to high computational effort. Please, choose a smaller power, or a larger alpha, or a larger RR."),call. =FALSE)           
                     }  
                            }
                 }## close 2


if(Tailed=="two"){ ## open 3
cst<- qbinom(1-alpha,N,1/(1+z)) ; if(1-pbinom(cst-1,N,1/(1+z))>alpha){cst<- cst+1}; powcheck<- 1-pbinom(cst-1,N,1/(1+z/min(RR)))
if(powcheck<power_inf&N<Nauxpow2){
if(Nauxpow2<=240){
stop(c("For the desired alpha, power[1], RR[1], and the specified z, the maximum length of surveillance, 'N', must be greater than or equal to ",Nauxpow2),call. =FALSE)
                }else{
           stop(c("For the desired alpha, power[1], RR[1], and the specified z, the required 'N' is greater than ",Nauxpow2, " which can lead to high computational effort. Please, choose a smaller power[1], or a larger alpha, or a smaller RR[1]."),call. =FALSE)           
                     }  
                  }

cst<- qbinom(1-alpha,N,1/(1+z)) ; if(1-pbinom(cst-1,N,1/(1+z))>alpha){cst<- cst+1}; powcheck<- 1-pbinom(cst-1,N,1/(1+z/max(RR)))
if(powcheck<power&N<Nauxpow){
if(Nauxpow<=Nref){
stop(c("For the desired alpha, power[2], RR[2], and the specified z, the maximum length of surveillance, 'N', must be greater than or equal to ",Nauxpow),call. =FALSE)
                }else{
           stop(c("For the desired alpha, power[2], RR[2], and the specified z, the required 'N' is greater than ",Nauxpow, " which can lead to high computational effort. Please, choose a smaller power[2], or a larger alpha, or a larger RR[2]."),call. =FALSE)           
                     }  
                  }

                 }## close 3



if(N>130&N<=170&length(GroupSizes)>20){message("For this N value the computation can take two or three hours.")}
if(N>170&length(GroupSizes)>20){message("For this N value the computation can take a day.")}


if(length(GroupSizes)==1){
   if(GroupSizes=="n"){GroupSizes<- seq(1,N)}else{GroupSizes<- seq(GroupSizes,N,GroupSizes)}
                         }else{G<- length(GroupSizes) ; GroupSizes<- GroupSizes%*%(upper.tri(matrix(0,G,G),diag=T)*1);GroupSizes<- GroupSizes[1,1:G]}
   

# auxiliary function

posi<- function(ii,jj){return( ii*(ii-1)/2+jj )}


#################
################# Function for relative risk estimate  <<<<IMPORTANT: NEEDS NOTATION ADJUSTMENTS
#################

### MAXIMUM LILKELIHOOD ESTIMATE FOR R
MLE_R<- function(cases,SampleSizes,z)
{
recand<- matrix(seq(0.01,10,0.01),,1)
lr<- function(rr){
return(prod( choose(SampleSizes,cases)*((1/(1+z/rr))^(cases))*((1-1/(1+z/rr))^(SampleSizes-cases)) ))
}
veccand<- apply(recand,1,lr)
Rhat<- seq(0.01,10,0.01)[veccand==max(veccand)]
return(Rhat)
} 





#################
################# Function that calculates parameters order
#################


### Auxiliar function for calculation of parameters order
countj<- function(ii,jj){if(jj<=ks_inf(ii)| (jj>ks_inf(ii)&j<ks(ii)) ){return(jj+1)}else{return(jj-ks(ii)+1+ks_inf(ii)+1)}}

### Function that calculates paramters order 
posi2<- function(ii,jj)
{

if(Tailed!="two"){
if(ii==imin){return(jj-ks(ii)+1)}else{
return( (imin+(ii-1))*(ii-imin)/2 - sum(apply(matrix(seq(imin,ii-1),,1),1,ks))+(ii-imin) + jj-ks(ii)+1)
                                     }
                 }else{ 
                       if(ii==imin){ return( countj(ii,jj) ) }else{
 return( (imin+(ii-1))*((ii-1)-imin+1)/2 - sum(apply(matrix(seq(imin,ii-1),,1),1,ks)) + (ii-1)-imin+1   +   sum(apply(matrix(seq(imin,ii-1),,1),1,ks_inf))+(ii-1)-imin+1  + countj(ii,jj) )
                                                                   }
                      }

}



################
################ <<< HERE STARTS THE SEARCH FOR THE SOLUTION DEPENDING ON Objective="ETimeToSignal" or Objective="ESampleSize"
################


if(Objective=="ETimeToSignal"){ NN<- N}else{if(Tailed!="two"){NN<- Nauxpow}else{NN<- max(Nauxpow,Nauxpow2)}}

Nin<- N 
obj_old<- Nin
N<- NN-1
auxnn<- 0

while(auxnn==0){ ## 99
N<- N+1
# total number of parameters

M<- (1+N)*N/2

if(length(R0)==1){p0<- 1/(1+z/R0)}else{p0<- 1/(1+z/R0[1]); p02<- 1/(1+z/R0[2])}
if(Tailed=="two"){pA<- 1/(1+z/RR[2]) ; pA_inf<- 1/(1+z/RR[1])}else{pA<- 1/(1+z/RR)}


imin<- 1
pini<- 1-pbinom(0,1,p0)
while(pini>alpha){imin<- imin+1;pini<- 1-pbinom(imin-1,imin,p0)}

imin_inf<- 1
pini<- pbinom(0,1,p0)
while(pini>alpha){imin_inf<- imin_inf+1;pini<- pbinom(0,imin_inf,p0)}

imin<- min(imin,imin_inf)

# auxiliary function

ks<- function(ii){return(qbinom(1-alpha,ii,p0))}
ks_inf<- function(ii){if(Tailed=="two"){ axx<- qbinom(alpha,ii,p0); if(pbinom(axx-1,ii,p0)>alpha|pbinom(axx-1,ii,p0)==0){return(0)}else{return(axx+1)} }else{return(-2)}}

# auxiliary to define the objective function to be optimized (the expected time to signal)
#auxETS<- rep(seq(1,N),seq(1,N,1)) 



## Restriction of the type <= with objects Ar1 and a1

if(Tailed!="two"){M2<- (imin+N)*(N-imin+1)/2 - sum(apply(matrix(seq(imin,N),,1),1,ks)) + N-imin+1}else{
                  M2<- (imin+N)*(N-imin+1)/2 - sum(apply(matrix(seq(imin,N),,1),1,ks)) + N-imin+1   +   sum(apply(matrix(seq(imin,N),,1),1,ks_inf))+N-imin+1
                                                                                                      }

Ar1s<- matrix(0,M,M2)

Ar1<- matrix(0,M2+2,M2)
a1<- rep(1,M2+2)

prob0<- rep(0,M2) # prob0[posi(i,j)] probability of having, at time i, exactly j cases.
probA<- rep(0,M2); if(Tailed=="two"){probA_inf<- rep(0,M2)}
probINT<- rep(0,M2) 
auxETS<- rep(0,M2); if(Tailed=="two"){auxETS2<- rep(0,M2)}

if(Tailed=="two"){jin<- 0}else{jin<- 1}

for(i in 1:N){
   for(j in jin:i){


            if( i==imin & (j>=ks(i)|j<=ks_inf(i)+1) &sum(i==GroupSizes)>0){
      
                if(j<i){Ru<- z*j/(i-j)+ConfIntWidth/2 ; Rl<- z*j/(i-j)-ConfIntWidth/2}else{Ru<- z+ConfIntWidth/2; Rl<- z-ConfIntWidth/2}
                    pRu<- 1/(1+z/Ru) ; pRl<- 1/(1+z/Rl)

                    auxETS[posi2(i,j)]<- i  ; if(Tailed=="two"){auxETS2[posi2(i,j)]<- j}
                    if(length(R0)==1){prob0[posi2(i,j)]<-  dbinom(j,i,p0)}else{prob0[posi2(i,j)]<-  max(dbinom(j,i,p0),dbinom(j,i,p02))}
                    probA[posi2(i,j)]<-  dbinom(j,i,pA) ; if(Tailed=="two"){probA_inf[posi2(i,j)]<-  dbinom(j,i,pA_inf)}
                    if(pRl>0&Rl>0){probINT[posi2(i,j)]<-  dbinom(j,i,pRu)+dbinom(j,i,pRl)}else{probINT[posi2(i,j)]<-  dbinom(j,i,pRu)}
                    Ar1s[posi(i,j),posi2(i,j)]<- 1 
                        Ar1[posi2(i,j),]<- Ar1s[posi(i,j),]
                        a1[posi2(i,j)]<- 1
 
                                                                                     }


            if(i>imin){
         if(j<i){Ru<- z*j/(i-j)+ConfIntWidth/2 ; Rl<- z*j/(i-j)-ConfIntWidth/2}else{Ru<- z+ConfIntWidth/2; Rl<- z-ConfIntWidth/2}
                    pRu<- 1/(1+z/Ru) ; pRl<- 1/(1+z/Rl)

     Ar1s[posi(i,j),]<- (choose(i-1,j-1)*Ar1s[posi(i-1,j-1), ]+choose(i-1,j)*Ar1s[posi(i-1,j), ])/choose(i,j)

              if( (j>=ks(i)|j<=ks_inf(i)+1) &sum(i==GroupSizes)>0){

                    auxETS[posi2(i,j)]<- i ; if(Tailed=="two"){auxETS2[posi2(i,j)]<- j}
                    if(length(R0)==1){prob0[posi2(i,j)]<-  dbinom(j,i,p0)}else{prob0[posi2(i,j)]<-  max(dbinom(j,i,p0),dbinom(j,i,p02))}
                    probA[posi2(i,j)]<-  dbinom(j,i,pA) ; if(Tailed=="two"){probA_inf[posi2(i,j)]<-  dbinom(j,i,pA_inf)}
                     if(pRl>0&Rl>0){probINT[posi2(i,j)]<-  dbinom(j,i,pRu)+dbinom(j,i,pRl)}else{probINT[posi2(i,j)]<-  dbinom(j,i,pRu)}                  
                        Ar1s[posi(i,j),posi2(i,j)]<- 1 
                        Ar1[posi2(i,j),]<- Ar1s[posi(i,j),]
                        a1[posi2(i,j)]<- 1

                                                                }
           
                     }

           
           }
               }


## Restriction of the type <=

Ar1[M2+1,]<- prob0
Ar1[M2+2,]<- probINT
a1[M2+1]<- alpha
a1[M2+2]<- 1-gama

## Restriction of the type ==

Ar3<- probA; a3<- max(power); if(Tailed=="two"){Ar3<- rbind(probA,probA_inf); a3<- c(a3,power_inf)}

## Coefficients of the objective function

if(Tailed=="two"){a<- auxETS*(probA+probA_inf)/2}else{a<- auxETS*probA}

## Finding the optimal result
#library(boot)

#memory.limit(size=6212)

rm(Ar1s)


result_new<- simplex(a, A1=Ar1, b1=a1,A2=NULL, b2=NULL , A3=Ar3 , b3=a3,eps = 0.01)



op_res<- result_new$soln # pis
powerres<- sum(probA*op_res)

if(Tailed!="two"){
if(N==NN){result_old<- result_new; probA_old<- probA; prob0_old<- prob0; auxETS_old<- auxETS}
if(Objective=="ETimeToSignal"){obj<- sum(auxETS*probA*op_res)/powerres}else{obj<- sum(auxETS*probA*op_res)+N*(1-powerres)}
if((obj>obj_old&result_new$solved==1)|N==Nin){auxnn<- 1; result<- result_old; probA<- probA_old; prob0<- prob0_old; auxETS<- auxETS_old}else{result_old<- result_new; obj_old<- obj; probA_old<- probA; prob0_old<- prob0; auxETS_old<- auxETS}
                 }else{
if(N==NN){result_old<- result_new; probA_old<- probA; probA_inf_old<- probA_inf; prob0_old<- prob0; auxETS_old<- auxETS; auxETS2_old<- auxETS2}
powerres2<- sum(probA_inf*op_res) 
if(Objective=="ETimeToSignal"){obj<- sum(auxETS*probA*op_res)/(2*powerres)+sum(auxETS*probA_inf*op_res)/(2*powerres2)}else{obj<- (sum(auxETS*probA*op_res)+N*(1-powerres))/2+(sum(auxETS*probA_inf*op_res)+N*(1-powerres2))/2}
if((obj>obj_old&result_new$solved==1)|N==Nin){auxnn<- 1; result<- result_old; probA<- probA_old; probA_inf<- probA_inf_old; prob0<- prob0_old; auxETS<- auxETS_old; auxETS2<- auxETS2_old; N<- N-1}else{result_old<- result_new; obj_old<- obj; probA_old<- probA; probA_inf_old<- probA_inf ; prob0_old<- prob0; auxETS_old<- auxETS; auxETS2_old<- auxETS2}
                      }  

} ## CLOSES 99


op_res<- result$soln # pis
powerres<- sum(probA*op_res) ; if(Tailed=="two"){powerres2<- sum(probA_inf*op_res) }
if(Tailed!="two"){ETimeToSignal<- sum(auxETS*probA*op_res)/powerres ; ELS<- sum(auxETS*probA*op_res)+N*(1-powerres)}else{
ETimeToSignal<- sum(auxETS*probA*op_res)/(2*powerres)+sum(auxETS*probA_inf*op_res)/(2*powerres2)
ELS<- (sum(auxETS*probA*op_res)+N*(1-powerres))/2+(sum(auxETS*probA_inf*op_res)+N*(1-powerres2))/2
                                                                                                                        }

if(Tailed!="two"){ #1                           
alpha_spending<- rep(0,N); cases.threshold<- rep(0,N)
for(i in imin:N){alpha_spending[i]<- sum(op_res[auxETS==i]*prob0[auxETS==i])}
alpha_spending<- as.numeric(alpha_spending%*%(upper.tri(matrix(rep(0,length(alpha_spending)*length(alpha_spending)),length(alpha_spending),),diag=T)*1))

## Reporting only the alpha spending and threshold at testing times
alpha_spending<- alpha_spending[GroupSizes]

## Calculating the optimal sample size in the case of continuous sequential analysis
if(length(Groups)==1){if(Groups=="n"|Groups==1){Nmax<- 1; while(alpha-alpha_spending[Nmax]>0.0001&Nmax<N){Nmax<- Nmax+1} }}else{Nmax<- N}

if(N<Nin&length(Groups)==1){alpha_spending[(N+1):(Nin)]<- alpha}

## Preparing objects to deliver results

res<- list(alpha_spending,ETimeToSignal,ELS,powerres,result$solved,Nmax-1)

if(Objective=="ETimeToSignal"){
names(res)<- c("optimal_alpha_spending","minETimeToSignal","ESampleSize","Power","solved","Optimal_N")
                              }else{ names(res)<- c("optimal_alpha_spending","ETimeToSignal","minESampleSize","Power","solved","Optimal_N")  }

plot(GroupSizes,alpha_spending,type="b",xlab="n",ylab="Alpha spending")

              }# close #1



if(Tailed=="two"){ #2
alpha_spending_lower<- rep(0,N); alpha_spending_upper<- rep(0,N)
   for(i in imin:N){
  for(j in 0:i){
if(j<=ks_inf(i)){alpha_spending_lower[i]<- alpha_spending_lower[i]+op_res[auxETS==i&auxETS2==j]*prob0[auxETS==i&auxETS2==j]}
if(j>=ks(i)){alpha_spending_upper[i]<- alpha_spending_upper[i]+op_res[auxETS==i&auxETS2==j]*prob0[auxETS==i&auxETS2==j]}
               }
                  }
alpha_spending_lower<- as.numeric(alpha_spending_lower%*%(upper.tri(matrix(rep(0,length(alpha_spending_lower)*length(alpha_spending_lower)),length(alpha_spending_lower),),diag=T)*1))          
alpha_spending_upper<- as.numeric(alpha_spending_upper%*%(upper.tri(matrix(rep(0,length(alpha_spending_upper)*length(alpha_spending_lower)),length(alpha_spending_upper),),diag=T)*1))                

## Reporting only the alpha spending and threshold at testing times
alpha_spending_lower<- alpha_spending_lower[GroupSizes]
alpha_spending_upper<- alpha_spending_upper[GroupSizes]

## Calculating the optimal sample size in the case of continuous sequential analysis
if(Groups=="N"|Groups==1){Nmax<- 1; while(alpha-alpha_spending_lower[Nmax]-alpha_spending_upper[Nmax]>0.0001&Nmax<N){Nmax<- Nmax+1} }else{Nmax<- N}

## Preparing objects to deliver results

res<- list(alpha_spending_lower,alpha_spending_upper,ETimeToSignal,ELS,powerres,result$solved,Nmax-1)

if(Objective=="ETimeToSignal"){
names(res)<- c("optimal_alpha_spending_lower","optimal_alpha_spending_upper","minETimeToSignal","ESampleSize","Power","solved","Optimal_N")
                              }else{ names(res)<- c("optimal_alpha_spending_lower","optimal_alpha_spending_upper","ETimeToSignal","minESampleSize","Power","solved","Optimal_N")  }

par(mfrow=c(2,1))
plot(GroupSizes,alpha_spending_upper,type="b",xlab="n",ylab="Alpha spending",main="Alpha spending for the upper threshold")
plot(GroupSizes,alpha_spending_lower,type="b",xlab="n",ylab="Alpha spending",main="Alpha spending for the lower threshold")



               }# close #2 

return(res)

}# CLOSE FUNCTION FOR OPTIMAL ALPHA SPENDING








