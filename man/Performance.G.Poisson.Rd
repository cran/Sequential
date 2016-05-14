\name{Performance.G.Poisson}
\alias{Performance.G.Poisson}
\title{Calculates Statistical Power, Expected Time to Signal and Expected Sample Size for Group Sequential Analysis with Poisson Data.}
\description{The \code{Performance.G.Poisson} function calculates statistical power, expected time to signal and the expected sample size
for group sequential analysis with Poisson data. It is required that the number of looks at the data is pre-specified, and equally spaced in terms of the sample size.
This function can also be used to  calculate the statistical significance level alpha, by setting RR=1, in which case the power output value is the alpha level.}
\usage{
Performance.G.Poisson(SampleSize,cv,GroupSizes,M=1,RR=2)
}
\arguments{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) expressed in terms of the expected number of events under the null hypothesis. 
SampleSize must be greater than 0. There is no default value.}
\item{cv}{The critical value for a group sequential analysis with Poisson data. This is calculated using the \code{\link[Sequential]{CV.G.Poisson}}.
There is no default value.}
\item{GroupSizes}{Vector containing the expected number of events under H0 for each test. The values must be positive numbers. The dimension of this vector must be equal to the maximum number of sequential tests. Thus, the sum of the entries in GroupSizes has to be equal to SampleSize. There is no default value.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. M must be a positive integer, and the default value is M=1. According to Kulldorff and Silva (2015), a reasonable choice is M=4.}
\item{RR}{The relative risk under the alternative hypothesis. Can be any number greater or equal to 1. The default value is RR=2.}
}
\details{
The function \code{Performance.G.Poisson} calculates the statistical power, the expected time to signal when the null hypothesis is rejected, and
the expected sample size at the end of the sequential analysis, whether the null hypothesis was rejected or not. The user specifies the upper limit
on the sample size (SampleSize), the upper boundary critical value (cv), a vector with the expected number of cases under H0 at each test (GroupSizes), the minimum number of
events required before rejecting the null hypothesis (M) and the relative risk under the alternative hypothesis (RR).

The calculations are done using an upper rejection boundary that is flat with respect to the log likelihood ratio, as expressed using the critical value.
The different looks at the data, or the different tests, are assumed to be equally spaced with respect to the accrual of the sample size.

To avoid very large computation times, we suggest not using SampleSize values greater than 1000. Typically, this is not a major restriction.
For example, for RR=1.1 and alpha=0.01, the statistical power is approximately 1 for a maximum sample size greater than 500.
}
\value{
\item{Power}{The statistical power.}
\item{ESignalTime}{The expected time to signal given that the null hypothesis is rejected.}
\item{ESampleSize}{The expected sample size when the sequential analysis ends (length of surveillance), whether the null hypothesis is rejected or not.}
}
\author{Ivair Ramos Silva and Martin Kulldorff
}
\references{
Silva IR, Kulldorff M. (2015). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, in press.
}
\section{Acknowledgements}{
Development of the Performance.G.Poisson function was funded by:\cr
-	Food and Drug Administration, Center for Biologics Evaluation and Research, through the Mini-Sentinel Post-Rapid Immunization Safety Monitoring (PRISM) program (v1.0).\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil (v1.0).\cr
-	Bank for Development of the Minas Gerais State (BDMG), Brazil (v1.0).\cr
-     National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0.1,2.0.2).  
}

\section{See also}{
\code{\link[Sequential]{CV.G.Poisson}}: Calculates critical values for group sequential analysis with Poisson data.\cr 
\code{\link[Sequential]{Performance.Poisson}}: Calculates the statistical power, expected time to signal and expected sample size for continuous sequential analysis with Poisson data.\cr
\code{\link[Sequential]{Performance.G.Binomial}}: Calculates the statistical power, expected time to signal and expected sample size for group sequential analysis with binomial data.
}
\examples{
## First use the CV.G.Poisson function to calculate the critical value for 5 sequential
## looks  at the data, spaced six units apart, and with a statistical significance level
## of 0.05:

cvt<- CV.G.Poisson(SampleSize=30,alpha=0.05,GroupSizes=c(6,6,6,6,6))

## For an alternative hypothesis of a relative risk of RR=1.5, calculates the statistical 
## power, the expected time to signal, and the expected sample size at the 
## end of the sequential analysis.

(Performance.G.Poisson(SampleSize=30,cv=cvt,GroupSizes=c(6,6,6,6,6),RR=1.5))
}

