\name{SampleSize.CondPoisson}
\alias{SampleSize.CondPoisson}
\title{Sample size calculation for the continuous sequential CMaxSPRT for Poisson data with limited information from historical cohort.}
\description{
The function \code{SampleSize.CondPoisson} obtains the required sample size (length of surveillance) needed
to guarantee a desired statistical power for a pre-specified relative risk, when doing continuous sequential CMaxSPRT, using a Wald-type upper boundary,
which is flat with respect to the likelihood ratio function.}
\usage{
SampleSize.CondPoisson(cc,D=0,M=1,alpha=0.05,power=0.9,RR=2)
}
\arguments{
\item{cc}{The total number of observed adverse events in the historical data. There is no default.}
\item{D}{The minium number for the ratio \eqn{P_k/V} before the null hypothesis can be rejected. The default value is \eqn{D=0}. A delayed start with \eqn{D>0} is recommended to avoid signaling very early on such that very little information would be available to judge whether the signal is more likely to be a true signal or chance finding.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. The default value is M=1. A delayed start with \eqn{M>1} is recommended to avoid signaling very early on such that very little information would be available to judge whether the signal is more likely to be a true signal or chance finding.}
\item{alpha}{The significance level. It must be in the range (0,0.5]. The default value is alpha=0.05.}
\item{power}{The target overall statistical power to detect an increased relative risk (RR). The default value is power=0.90.}
\item{RR}{The target relative risk to be detected with the requested statistical power. The default value is RR=2.}
}
\details{
When using the CMaxSPRT and the \code{CV.CondPoisson} function to conduct continuous sequential analysis
for Poisson data and limited historical data, the null hypothesis is rejected when the log likelihood ratio
exceeds the predetermined critical value calculated by \code{CV.CondPoisson}. The sequential analysis ends
without rejecting the null hypothesis when a predetermined upper limit on the sample size is reached, expressed 
either in terms of the ratio of the cumulative person-time in the surveillance population divided by the total
cumulative person-time in historical data (StopType="Tal"), or in terms of the observed number of events in the
surveillance population (StopType="Cases"). For example, the sequential analysis may end as soon as the sample size
is such that the cumulative person-time in the surveillance population is twice the cumulative person-time in historical data,
or there are 50 observed events in the surveillance population. The function \code{SampleSize.CondPoisson} calculates what
the upper limit on the sample size (length of surveillance) that is required for the continuous CMaxSPRT to achieve the desired
statistical power for a pre-specified relative risk RR. It outputs the upper limit on sample size for both definitions of the
surveillance length, one expressed in terms of the ratio of the cumulative person-time in the surveillance population divided by
the total cumulative person-time in historical data (\eqn{T}), and the other one expressed in terms of the observed number of
events in the surveillance population (\eqn{K}). To save computing time, the liberal computation approach proposed by Silva et al. (2016)
is used in \code{SampleSize.CondPoisson} to find the solution.
}

\value{
\item{T}{The minimum required sample size (length of surveillance), expressed in terms of the ratio of the cumulative person-time in the surveillance population divided by the total cumulative person-time in historical data, needed to achieve the desired statistical power for the pre-specified relative risk (RR).}
\item{K}{The minimum required sample size (length of surveillance), expressed in terms of the observed number of events in the surveillance population, needed to achieve the desired statistical power for the pre-specified relative risk (RR).}
\item{cv}{the corresponding critical value}
}

\author{Ivair Ramos Silva, Lingling Li}
\references{
Li L, Kulldorff M. (2010). A conditional maximized sequential probability ratio test for pharmacovigilance. Statistics in Medicine, 30;29(2):284-95.\cr
Silva IR, Li L, Kulldorff M. (2016). Exact conditional maximized sequential probability ratio test adjusting for covariates. Working paper.
}
\section{Acknowledgements}{
Development of the \code{SampleSize.CondPoisson} function was funded by:\cr
- National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0.1, v2.0.2).
}

\keyword{Continuous CmaxSPRT sample size}

\section{See also}{
\code{\link[Sequential]{CV.CondPoisson}}: calculating the critical value for continuous CMaxSPRT.\cr
\code{\link[Sequential]{Performance.CondPoisson}}: calculating the statistical power, expected time to signal and expected time of analysis for continuous CMaxSPRT.\cr
}

\examples{

# Sample size required to obtain a power of 90\%, for a relative risk of 1.5,
# no delay for starting the surveilla# nce (D=0), under an alpha level of 5\%,
# with 5 events in the historical data.

# res<- SampleSize.CondPoisson(cc=50,D=0,M=1,alpha=0.05,power=0.9,RR=2)

# If we write:
# res
# then we get:

# $K
# [1] 66

# $Tal
# [1] 0.7890625

# $cv
# [1] 3.626436

}

