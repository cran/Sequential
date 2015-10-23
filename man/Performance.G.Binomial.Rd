\name{Performance.G.Binomial}
\alias{Performance.G.Binomial}
\title{Statistical Power, Expected Signal Time and Expected Length of Surveillance for Group Sequential Analysis for Binomial Data.}
\description{The function \code{Performance.G.Binomial} calculates several performance 
metrics for the group binomial MaxSPRT for fixed upper limit on the sample 
size ("N"), critical value ("cv"), maximum number of tests to be perfomed, minimum number of events 
required before rejecting the null hypothesis (M), and true relative risk (RR). 
The metrics calculated are the statistical power, the expected time to signal when the null hypothesis is 
rejected and the expected sample size at the end of the analysis whether the null hypothesis was rejected or not.}
\usage{
Performance.G.Binomial(N,M=1,cv,z=1,RR=2,GroupSizes)
}
\arguments{
\item{N}{The upper limit on the sample size (length of surveillance) expressed in terms of the total number of events. 
"N" must be a positive integer. To avoid very large computation times, we suggest not using values greater than 1000. Typically,
this is not a major restriction. For example, for "RR=1.1" and "alpha=0.01" and "z=1", the statistical power is approximately 1 for "N>500". There is no default value.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. "M" must be a positive integer, and the default value is M=1.}
\item{cv}{Critical value, defining the upper rejection boundary. The null hypothesis is rejected when the log-likelihood value is greater than "cv". The "cv" parameter is usually obtained by first running \code{CV.Binomial}.It must be a positive number, and there is no default.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case. For example, if there are 3 controls matched to each case, z=3.  In a self-control analysis, z is the ratio of the control interval to the risk interval. For example, if the risk interval is 2 days long and the control interval is 7 days long, z=7/2. In terms of p, the binomial probability under the null hypothesis, p=1/(1+z), or equivalently, z=1/p-1. The parameter z must be a positive number. The default value is z=1.}
\item{RR}{The relative risk under the alternative hypothesis. If p and q are binomial probabilities under the null and the alternative hypotheses, respectively, then "RR=(1/p-1)/(1/q-1)".}
\item{GroupSizes}{Vector with the number of events (cases+controls) between two consecutive looks (tests) at the data, i.e, the group sizes. The
length of this vector is equal to the maximum number of looks. The entries do not have the same, but they sum up to N.
If the group sizes is an integer instead of a vector, then that integer is the group size for all looks at the data, and the number of looks is "N/GroupSizes".}
}
\details{
The function \code{Performance.G.Binomial} calculates the statistical power, the expected time to signal when the null hypothesis is rejected and the expected
sample size by using the group sequential MaxSPRT for binomial data given a fized relative risk (RR). When the null hypothesis is true, the probability of
having a case, instead of a control, is \eqn{p=1/(1+z)}. But, if the null hypothesis is false, and the true relative risk is a value RR>1, then the
probability of having a case is \eqn{p=RR/(RR+z)}.

If the user wants to calculate the actual overall size (or the maximum Type I error probability) of the test, that can be done by setting RR=1, in which case the power output value is the test size.
This is made for a given flat critical value, an upper limit on the sample size (N),
and a delayed start in the vigilance given by M.
 }
\value{
\item{Power}{Statistical power.}
\item{ESignalTime}{The expected time to signal given that the null hypothesis was rejected.}
\item{ESampleSize}{The expected sample size when the null was rejected or not.}
}
\author{Ivair Ramos Silva and Martin Kulldorff
}
\references{
Silva IR, Kulldorff M. (2015). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, in press.
}
\section{Acknowledgements}{
Development of the Performance.G.Binomial function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999. 
}

\section{See also}{
\code{\link[Sequential]{CV.G.Binomial}}: for calculating critical values in group sequential analysis with binomial data.\cr
\code{\link[Sequential]{Performance.Binomial}}: for calculating the statistical power, expected time to signal and expected time of analysis for continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{Analyze.Binomial}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion.\cr  
\code{\link[Sequential]{CV.Poisson}}: for calculating critical values for continuous sequential analysis with Poisson data.\cr 
\code{\link[Sequential]{Performance.G.Poisson}}: for calculating statistical power, expected time to signal and expected length of surveillancefor group sequential analysis with Poisson data.\cr
}
\examples{

result<- Performance.G.Binomial(N=40,M=1,cv=2.5,z=1,RR=2,GroupSizes=c(2))
# if you type:
result
# then you will get the following output:
# $Power
# [1] 0.6594118

# $`ESignalTime`
# [1] 17.18626

# $`ESampleSize`
# [1] 24.95635
}
