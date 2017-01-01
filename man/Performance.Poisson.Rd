\name{Performance.Poisson}
\alias{Performance.Poisson}
\title{Calculates Statistical Power, Expected Time to Signal and Expected Sample Size for Continuous Sequential Analysis with Poisson Data.}
\description{The \code{Performance.Poisson} function calculates three different performance metrics for the continuous sequential analysis with Poisson data:
the statistical power, the expected time to signal when the null hypothesis is rejected and the expected sample size at the end of the analysis whether the
null hypothesis was rejected or not. The user specifies the relative risk under the alternative hypothesis (RR), as well as the sequential analysis parameters. 
To calculate the statistical significance level alpha, RR=1, in which case the power output value is the alpha level.}
\usage{
Performance.Poisson(SampleSize, D = 0, M = 1, cv, RR = 2)
}
\arguments{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) expressed in terms of the expected number of events under the null hypothesis.
The SampleSize must be greater than 0. There is no default value.}
\item{D}{The expected number of events under the null hypothesis at the first look at the data. The default is D=0, which is
also the best choice. This means that there is no delay in the start of the sequential analysis. It is required that D\eqn{\leq}SampleSize.}
\item{M}{The minimum number of events needed to be observed before the null hypothesis can be rejected. The default is M=1, 
which means that even a single event can reject the null hypothesis if it occurs sufficiently early. A good rule of thumb is to set M=4
(Kulldorff and Silva, 2015).}
\item{cv}{The critical value constituting the upper rejection boundary. This can be calculated using the \code{\link[Sequential]{CV.Poisson}} function.}
\item{RR}{The relative risk under the alternative hypothesis. It is required that RR\eqn{\geq 1}. The default value is RR=2.}
}
\details{
For continuous sequential analysis with Poisson data, the \code{Performance.Poisson} function calculates the statistical power, the expected time to signal when
the null hypothesis is rejected and the expected sample size until the analysis ends whether the null is rejected or not. The sample size is expressed in terms of
the expected number of events under the null hypothesis. Large values of the SampleSize, greater than say 1000, may leads to long computing times.
When the statistical power is close to 1, then the expected time to signal will be very close to the expected sample size, since both are measured in information time,
using the expected events under the null hypothesis as the unit.

To avoid very large computation times, we suggest not using values greater than 1000. Typically, this is not a major restriction. For example,
for RR=1.1 and alpha=0.01, the statistical power is approximately 1 for a maximum sample size greater than 500.
}
\value{
\item{Power}{The statistical power.}
\item{ESignalTime}{The expected time to signal given that the null hypothesis is rejected.}
\item{ESampleSize}{The expected sample size when the sequential analysis ends (length of surveillance) whether the null hypothesis was rejected or not.}
}
\author{Ivair Ramos Silva and Martin Kulldorff
}
\references{
Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Vaccine Safety Surveillance. Sequential Analysis, 30: 58--78.
Kulldorff M, Silva IR. (2015). Continuous post-market sequential safety surveillance with minimum events to signal. arxiv:1503.01978 [stat.ap].
}
\section{Acknowledgements}{
Development of the Performance.Poisson function was funded by:\cr
-	Food and Drug Administration, Center for Biologics Evaluation and Research, through the Mini-Sentinel Post-Rapid Immunization Safety Monitoring (PRISM) program (v1.0).\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil (v1.0).\cr
-	Bank for Development of the Minas Gerais State (BDMG), Brazil (v1.0).\cr
-     National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0.1,2.0.2).  
}

\section{See also}{
\code{\link[Sequential]{CV.Poisson}}: Calculates critical values for continuous sequential analysis with Poisson data.\cr 
\code{\link[Sequential]{SampleSize.Poisson}}: Sample size calculations for continuous sequential analysis with Poisson data.\cr
\code{\link[Sequential]{Performance.G.Poisson}}: Calculates the statistical power, expected time to signal and expected sample size for group sequential analysis with Poisson data.
}

\examples{
## Suppose we want to find the statistical power to detect a relative risk
## of 2 when doing up to at most 20 months of surveillance, as well as the
## expected time to signal when the null hypothesis is rejected. During 
## each month, we expected to see 0.5 events if the null hypothesis is true. 
## This means that the upper limit on the sample size is 20*0.5=10 expected
## events under the null hypothesis. We will then first calculate the critical
## value for an upper limit on the sample size equal to 10 and a significance
## level of alpha=0.05:

# cvt<- CV.Poisson(SampleSize=10,alpha=0.05)
# cvt
# [1] 3.467952

## After that, we use the Performance.Poisson function to calculate the
## power and the expected time to signal when the null hypothesis is
## rejected for the alternative hypothesis with a relative risk equal to 2:

## Power, expected signal time and expected sample size for a relative risk
## equal to 2:
# Performance.Poisson(SampleSize=10,cv=cvt,RR=2)
#         Power    ESignalTime    ESampleSize 
# [1,] 0.6850634     4.130985      5.979353)

## From the results, we see that the statistical power is 68.5%. When the null
## is rejected, the expected time to signal is 4.13 in the unit of events
## expected under the null. If data is collected uniformly over time at the
## rate of 0.5 expected counts per month,the expected time to signal is
## 4.13/0.5= 8.26 months.

## The above calculations can also be accomplished using one single command line:

#  Performance.Poisson(SampleSize=10,cv=CV.Poisson(SampleSize=10,alpha=0.05),RR=2)
#         Power    ESignalTime    ESampleSize
#  [1,] 0.6850634     4.130985      5.979353
}