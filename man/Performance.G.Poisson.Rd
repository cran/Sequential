\name{Performance.G.Poisson}
\alias{Performance.G.Poisson}
\title{Power, Expected Signal Time and Expected Length of Serveillance for the Poisson Group MaxSPRT.}
\description{The \code{Performance.G.Poisson} function calculates several performance 
metrics for the continuous Poisson based MaxSPRT. The user specifies the upper limit on the sample 
size (SampleSize), the upper boundary critical value (cv), the maximum number of tests to be perfomed (Looks), the minimum number of events 
required before rejecting the null hypothesis (M) and a relative risk (RR). 
The metrics calculated are the statistical power, the expected time to signal when the null hypothesis is 
rejected and the expected sample size at the end of the analysis whether the null hypothesis was rejected or not. If the user wants to calculate
the statistical significance level alpha, that can be done by setting RR=1, in which case the power output value is the alpha level.}
\usage{
Performance.G.Poisson(SampleSize,cv,Looks=1,M=1,RR=2)
}
\arguments{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) expressed in terms of the expected number of events under the null hypothesis. 
"SampleSize" must be greater than 0. To avoid very large computation times, we suggest not using values greater than 1000. Typically,
this is not a major restriction. For example, for "RR=1.1" and "alpha=0.01", the statistical power is approximately 1 for a maximum 
sample size greater than 500. There is no default value.}
\item{cv}{the critical value for a group sequential MaxSPRT test (see \code{\link[Sequential]{CV.G.Poisson}}).}
\item{Looks}{The maximum number of tests to be performed, which is, the number of equally spaced looks at the data if the null hypothesis is never rejected. The value must be an integer greater or equal to 1. The default value is 'Looks=1'.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. According to Kulldorff and Silva(2015), a reasonable choice is 'M=4', but the default value is 'M=1', and it must be a positive integer.}
\item{RR}{the relative risk (>=1) to calculate the power and the expected time to have a signal. The default value is 'RR=2'.}
}
\details{
The function \code{Performance.G.Poisson} calculates the power, the expected time to have a signal and the expected
time of surveillance by using the group sequential MaxSPRT for poisson data, given a flat boundary, an upper limit on the number "SampleSize"
of expected events under the null, and a delayed start in the vigilance given by "M".
It is important to point out that large values of "SampleSize", as "SampleSize>1000", for
example, leads to a high computational effort. 

The times when each test will be performed is defined by 
using regular time intervals that are referenced by constant expected number of adverse events between 
looks at the accumulated data. In another words, the MaxSPRT statistic is compared against cv whenever the time is a multiple of \eqn{SampleSize/Looks}. 
}
\value{
\item{Power}{the statistical power for a given flat upper boundary (CV).}
\item{ESignalTime}{the expected time to have a signal.}
\item{ESampleSize}{the expected sample size at the end of the sequential analysis.}
}
\author{Ivair Ramos Silva and Martin Kulldorff
}
\references{
Silva IR, Kulldorff M. (2015). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, in press.
}
\section{Acknowledgements}{
Development of the Performance.G.Poisson function was funded by:\cr
-	Food and Drug Administration, Center for Biologics Evaluation and Research, through the Mini-Sentinel Post-Rapid Immunization Safety Monitoring (PRISM) program (v1.0);\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil (v1.0);\cr
-	Bank for Development of the Minas Gerais State (BDMG), Brazil (v1.0).  
}

\section{See also}{
\code{\link[Sequential]{Performance.G.Binomial}}: for calculating the statistical power, expected time to signal and expected time of analysis in group sequential analysis with binomial data.\cr
\code{\link[Sequential]{SampleSize.Binomial}}: for calculating the minimum sample size given a target power in continuous sequential analysis with binomial data.\cr 
\code{\link[Sequential]{CV.Binomial}}: for calculating critical values in continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{Analyze.Binomial}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion.\cr  
\code{\link[Sequential]{CV.Poisson}}: for calculating critical values for continuous sequential analysis with Poisson data.\cr 
\code{\link[Sequential]{CV.G.Poisson}}: for calculating critical values for group sequential analysis with Poisson data.
}
\examples{
## Critical value for a significance level of 0.05:
cvt<- CV.G.Poisson(SampleSize=30,alpha=0.05,Looks=5,M=1)
# Power, expected signal time and expected surveillance time for a relative risk equal
# to 2:
(Performance.G.Poisson(SampleSize=30,cv=cvt,Looks=5,M=1,RR=1.5))
}

