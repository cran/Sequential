\name{Performance.CondPoisson}
\alias{Performance.CondPoisson}
\title{Statistical power, expected time to signal and expected sample size for the continuous sequential CMaxSPRT for Poisson data with limited information from historical cohort.}
\description{
The function \code{Performance.CondPoisson} calculates several performance metrics for the continuous CMaxSPRT for selected computation approach,
the type of upper limit definition and its value, critical value, number of historical data events, criteria of delayed start, the Type I error rate,
and a relative risk. The metrics calculated are the statistical power, the expected time to signal when the null hypothesis
is rejected, and the expected sample size at the end of the analysis whether the null hypothesis was rejected or not.}
\usage{
Performance.CondPoisson(Inference="exact",cv,StopType="Cases",
T="n",K="n",cc,D=0,M=1,RR=2)
}

\arguments{
\item{Inference}{Inference='liberal', 'exact', or 'conservative' for the computation approach. Inference='liberal' for the liberal approach with possibly underestimated critival values and higher-than-nominal Type I error rate, Inference='exact' for the exact approach with exact critival values and nominal Type I error rates, Inference='conservative' for the conservative approach with possibly overestimated critival values and lower-than-nominal Type I error rates. The default is Inference='exact'.}
\item{cv}{Critical value, defining the upper rejection boundary for the log-likelihood ratio test statistic. The null hypothesis is rejected when the log-likelihood value is greater than 'cv'. The 'cv' parameter is
usually obtained by first running \code{CV.CondPoisson}. It must be a positive number,and there is no default.}
\item{StopType}{StopType='Tal' or 'Cases' for the type of surveillance length definition. The default value is 'StopType=Cases'. See details.}
\item{T}{The upper limit on length of surveillance expressed in terms of the ratio of the cumulative person-time in the surveillance population divided by the total cumulative person-time in historical data, i.e., \eqn{P_k/V \leq T}. This argument \eqn{T} is used if and only if StopType='Tal'. There is no default value.}
\item{K}{The upper limit on length of surveillance expressed in terms of the observed number of events in the surveillance population, i.e., \eqn{k \leq K}. This argument \eqn{K} is used if and only if StopType='Cases'. There is no default value.}
\item{cc}{The total number of observed adverse events in the historical data. There is no default value.}
\item{D}{The minium number for the ratio \eqn{P_k/V} before the null hypothesis can be rejected. This argument is used together with StopType='Tal'. The default value is \eqn{D=0}.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. This argument is used together with StopType='Cases'. The default value is M=1.}
\item{RR}{The target relative risk for which statistical power, expected signal time and expected length of surveillance are calculated. The default is 'RR=2'.}
}
\details{
For continuous sequential analysis with Poisson data with limited historical information, the \code{Performance.CondPoisson} function calculates the
statistical power, the expected time to signal when the null hypothesis is rejected and the expected sample size until the analysis ends whether the
null is rejected or not. The sample size (i.e., the upper limit) can be expressed either in terms of the ratio "T" of the cumulative person-time in the surveillance population
divided by the total cumulative person-time in historical data (StopType="Tal"), i.e., \eqn{P_k/V \leq T}, or in terms of the observed number of events "K" in the surveillance
population (StopType="Cases"), i.e., \eqn{k \leq K}. Large values of the SampleSize, greater than say 1000, may leads to long computing times. When the statistical power
is close to 1, then the expected time to signal will be very close to the expected sample size.

For the parameter of delayed start, "D", a delayed start with \eqn{D>0} is recommended to avoid signaling very early on such that very little
information would be available to judge whether the signal is more likely to be a true signal or chance finding. Similarly, if the delayed start is defined in
terms of the number of events, "M", a setting such that \eqn{M>1} is recommended to avoid signaling very early on such that very little information would
be available to judge whether the signal is more likely to be a true signal or chance finding.
}

\value{
\item{Power}{The statistical power.}
\item{ESignalTime}{The expected time to signal given that the null hypothesis is rejected.}
\item{ESampleSize}{The expected sample size when the sequential analysis ends (length of surveillance) whether the null hypothesis was rejected or not.}
}
\author{Ivair Ramos Silva, Lingling Li}
\references{
Li L, Kulldorff M. (2010). A conditional maximized sequential probability ratio test for pharmacovigilance. Statistics in Medicine, 30;29(2):284--95.\cr

Silva IR, Li L, Kulldorff M. (2019). Exact conditional maximized sequential probability ratio test adjusting for covariates. Sequential Analysis, in press.

Silva IR., Lopes LM., Dias P., Yih WK. (2019). Alpha Spending for Historical Versus Surveillance Poisson Data With CMaxSPRT. Statistics in Medicine, DOI: 10.1002/sim.8097, 1--13. 

}
\section{Acknowledgements}{
Development of the \code{Performance.CondPoisson} function was funded by:\cr
- National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0.1, v2.0.2).
- Foundation for Research Support of Minas Gerais State (FAPEMIG), MG, Brazil, through the grant Demanda Universal.
}

\keyword{Continuous CmaxSPRT performance}

\section{See also}{
\code{\link[Sequential]{CV.CondPoisson}}: calculating the critical value for continuous CMaxSPRT.  \cr
\code{\link[Sequential]{SampleSize.CondPoisson}}: calculating the sample size for continuous CMaxSPRT.\cr
}

\examples{

## First calculate the critical value with upper limit defined in terms of
## the number of observed events in surveillance population (K=50), with 50
## events in historical data, no delayed start, and alpha=0.05: 
# res<-CV.CondPoisson(Inference="exact",StopType="Cases",K=20,cc=50,M=1,
# alpha=0.05)

# cvt<- res[[2]] 

## calculate the performance using the critical value 'cvt' from the previous
## step, under RR=1.5:
#Performance.CondPoisson(Inference="exact",cv=cvt,StopType="Cases",K=20,cc=50,
# M=1,RR=1.5) 

}

