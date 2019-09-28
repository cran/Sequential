\name{CV.CondPoisson}
\alias{CV.CondPoisson}
\title{Critical values for continuous sequential CMaxSPRT for Poisson data with limited information from historical cohort.}
\description{
The function \code{CV.CondPoisson} calculates critical values for the continuous sequential CMaxSPRT, using a Wald-type upper boundary,
which is flat with respect to the likelihood ratio function, and a pre-specified upper limit on surveillance length. 
}
\usage{
CV.CondPoisson(Inference="exact", StopType="Cases",T="n",K="n",cc,
D=0,M=1,alpha=0.05,Tailed="upper")
}

\arguments{
\item{Inference}{Inference='liberal', 'exact', or 'conservative' for the computation approach. Inference='liberal' for the liberal approach with possibly underestimated critival values and higher-than-nominal Type I error rate, Inference='exact' for the exact approach with exact critival values and nominal Type I error rates, Inference='conservative' for the conservative approach with possibly overestimated critival values and lower-than-nominal Type I error rates. The default is Inference="exact".}
\item{StopType}{StopType='Tal' or 'Cases' for the type of surveillance length definition. With StopType='Tal', the maximum surveillance length (i.e., the upper limit) is defined in terms of the ratio of the cumulative person-time in the surveillance population divided by the total cumulative person-time in historical data, i.e., \eqn{P_k/V \leq T}; with StopType='Cases', the maximum surveillance length is defined interms of the observed number of events in the surveillance population, i.e., \eqn{k \leq K}. The default is StopType="Cases"}
\item{cc}{The total number of observed adverse events in the historical data. There is no default value.}
\item{K}{The upper limit on length of surveillance expressed in terms of the observed number of events in the surveillance population, i.e., \eqn{k \leq K}. This argument \eqn{K} is used if and only if StopType='Cases'. There is no default value.}
\item{T}{The upper limit on length of surveillance expressed in terms of the ratio of the cumulative person-time in the surveillance population divided by the total cumulative person-time in historical data, i.e., \eqn{P_k/V \le T}. This argument \eqn{T} is used if and only if StopType='Tal'. There is no default value.}
\item{D}{The minium number for the ratio \eqn{P_k/V} before the null hypothesis can be rejected. This argument is used together with \eqn{T}. The default value is \eqn{D=0}. A delayed start with \eqn{D>0} is recommended to avoid signaling very early on such that very little information would be available to judge whether the signal is more likely to be a true signal or chance finding.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. This argument is used together with \eqn{K}. A delayed start with \eqn{M>1} is recommended to avoid signaling very early on such that very little information would be available to judge whether the signal is more likely to be a true signal or chance finding. The default value is M=1.}
\item{alpha}{The significance level, or the type 1 error probability, which is the probability of rejecting the null hypothesis when it is true. The alpha level  must be in the range (0,0.5]. The default value is alpha=0.05.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
}
\details{
For continuous sequential analysis with CMaxSPRT by Li and Kulldorff (2010) for Poisson data and limited historical data, \code{CV.CondPoisson} calculates the critical
value that constitutes the upper boundary used to determine if the null hypothesis should be rejected. This is done for pre-specified values
of the statistical significance level (alpha) and an upper limit which can be defined based on either the observed number of events, "K", or the ratio "T"
between the cumulative person-times in the surveillance population versus the historical data, as well as other parameter settings.

The test is one-sided, so that the null hypothesis is only rejected when there are more events than expected. 

Following the results of Silva et al. (2016), the function offers three computation approaches which calculate liberal, exact, and conservative critical values respectively.
When the upper limit is medium (e.g., \eqn{K=50}) or large, the computational requirements for the exact approach can be high.
The recommendation is to use the exact approach when the upper limit is small (e.g., \eqn{K=10}), use the conservative approach
when the upper limit is medium (\eqn{K=50}) or large but \eqn{cc} is small, and use the liberal approach when \eqn{cc} is medium (e.g., 50)
or large. Exact numerical results show that the three approaches yield very similar results when \eqn{K} and \eqn{cc} are reasonably large. 

}

\value{
\item{Type_I_Error}{The actual Type I error, for the exact approach. It equals the nominal level specified by the argument "alpha".
For the liberal approach, the actual Type I error rate may be higher than the specified nominal level. For the conservative approach,
the actual Type I error rate may be lower than the specified nominal level.}
\item{cv}{The critical value for a significance level. For the exact approach, it is the exact critical value, for the liberal
approach, it is the smallest liberal vaue, for the conservative approach, it is the largest conservative value.}
}
\author{Ivair Ramos Silva, Lingling Li}
\references{
Li L, Kulldorff M. (2010). A conditional maximized sequential probability ratio test for Pharmacovigilance. Statistics in Medicine, 29 (2), 284--295. 

Silva IR, Li L, Kulldorff M. (2019). Exact Conditional Sequential Testing for Poisson Data. Sequential Analysis, in press.

Silva IR., Lopes LM., Dias P., Yih WK. (2019). Alpha Spending for Historical Versus Surveillance Poisson Data With CMaxSPRT. Statistics in Medicine, DOI: 10.1002/sim.8097, 1--13. 

}
\section{Acknowledgements}{
Development of the \code{CV.CondPoisson} function was funded by:\cr
- National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0.1, v2.0.2).
}

\keyword{Continuous CmaxSPRT CV}

\section{See also}{
\code{\link[Sequential]{SampleSize.CondPoisson}}: calculating the upper limit with given Alpha, RR, and a desired power level for continuous CMaxSPRT.\cr
\code{\link[Sequential]{Performance.CondPoisson}}: calculating the statistical power, expected time to signal and expected
time of analysis for continuous CMaxSPRT.\cr
}

\examples{

# Calculates the exact critical value with upper limit of
# K=20 and a delayed start of a minimum of 2 cases, historical
# data has 20 events, and for a statistical significance level
# of 0.05. 
# res<- CV.CondPoisson(Inference="exact", StopType="Cases",K=20,cc=20,
# M=2,alpha=0.05)

# which gives the results:
# res
# $Type_I_Error
# [1] 0.05
# $cv
# [1] 3.149115

# Calculates the liberal critical value with a uppe limit of
# T=0.5 and a delayed start of D=0.1, i.e., the cumulative
# person-time in the surveillance population is at least
# one-tenth of the total cumulative person-time in historical
# data, historical data has 20 events, and for a statistical
# significance level of 0.05. 
# res2<- CV.CondPoisson(Inference="liberal",StopType="Tal",T=0.5,cc=20,
# D=0.1,alpha=0.05)

# which gives the results:
# res2
# $Type_I_Error
# [1] 0.05
# $cv
# [1] 2.874993


}

