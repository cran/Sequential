\name{Performance.Poisson}
\alias{Performance.Poisson}
\title{Calculates Statistical Power, Expected Signal Time and Sample Size for Continuous Sequential Analysis with Poisson Data.}
\description{The \code{Performance.Poisson} function calculates several performance 
metrics for the continuous Poisson based MaxSPRT. The user specifies the upper limit on the sample 
size (SampleSize), the delayed start of analysis (D), the minimujm number of events required before rejecting the null hypothesis (M), the upper boundary critical value (cv) and a relative risk (RR).
The metrics calculated are the statistical power, the expected time to signal when the null hypothesis is 
rejected and the expected sample size at the end of the analysis whether the null hypothesis was rejected or not. If the user wants to calculate
the statistical significance level alpha, that can be done by setting RR=1, in which case the power output value is the alpha level.}
\usage{
Performance.Poisson(SampleSize, D = 0, M = 1, cv, RR = 2)
}
\arguments{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) expressed in terms of the epected number of events under the null hypothesis. 
"SampleSize" must be greater than 0. To avoid very large computation times, we suggest not using values greater than 1000. Typically,
this is not a major restriction. For example, for "RR=1.1" and "alpha=0.01", the statistical power is approximately 1 for a maximum 
sample size greater than 500.}
\item{D}{The expected number of events under the null hypothesis before the first look at the data. The default is "D=0", which is
also the best choice. This means that there is no delay in the start of the sequential analyses. It is required that 'D<=SampleSize'. }
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. The default is 'M=1', 
which means that even a single event can reject the null hypothesis if it occurs sufficiently early. A good rule of thumb is to set 'M=4'
(Kulldorff and Silva, 2015).}
\item{cv}{The critical value constituting the upper rejection boundary. This can be calculated using the \code{\link[Sequential]{CV.Poisson}}) function.}
\item{RR}{The relative risk (>=1) for which statistical power, expected signal time and expected length of surveillance is calculated. The default is RR=2.}
}
\details{
For the continuous Poisson based MaxSPRT, the \code{Performance.Poisson} function calculates the statistical power, the expected time to signal 
when the null hypothesis is rejected and the expected sample size until the analysis ends whether the null is rejected or not. 
Large values of "SampleSize", such as "SampleSize>1000", leads to a long comcputing times.
}
\value{
\item{Power}{The statistical power.}
\item{ESignalTime}{The expected time to signal given that the null hypothesis is rejected.}
\item{ESampleSize}{The expected sample size when the analysis ends (length of surveillance) whether the null hypothesis was rejected or not.}
}
\author{Ivair Ramos Silva and Martin Kulldorff
}
\references{
Kulldorff, M. and Silva, I. R. (2015). Continuous post-market sequential safety surveillance with minimum events to signal. arxiv:1503.01978 [stat.ap].
}
\section{Acknowledgements}{
Development of the Performance.Poisson function was funded by:\cr
-	Food and Drug Administration, Center for Biologics Evaluation and Research, through the Mini-Sentinel Post-Rapid Immunization Safety Monitoring (PRISM) program (v1.0);\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil (v1.0);\cr
-	Bank for Development of the Minas Gerais State (BDMG), Brazil (v1.0).  
}

\section{See also}{
\code{\link[Sequential]{CV.Poisson}}: for calculating critical values for continuous sequential analysis with Poisson data.\cr 
\code{\link[Sequential]{CV.G.Poisson}}: for calculating critical values for group sequential analysis with Poisson data.\cr
\code{\link[Sequential]{Performance.G.Poisson}}: for calculating the statistical power, expected time to signal and expected time of analysis in group sequential analysis with Poisson data.\cr
\code{\link[Sequential]{CV.Binomial}}: for calculating critical values in continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{Analyze.Binomial}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion.  
}

\examples{
# Critical value for an upper limit on the sample size equal to 10 and a significance 
# level equal to 0.05:
#cvt<- CV.Poisson(SampleSize=10,D=0,M=1,alpha=0.05)
## Power, expected signal time and expected sample size for a relative risk equal to 2:
#(Performance.Poisson(SampleSize=10,D=0,M=3,cv=cvt,RR=2)) 
}

