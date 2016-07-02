\name{Performance.Binomial}
\alias{Performance.Binomial}
\title{Statistical Power, Expected Signal Time and Sample Size for Continuous Sequential Analysis with Binomial Data.}
\description{The function \code{Performance.Binomial} calculates several performance 
metrics for the continuous binomial MaxSPRT for fixed upper limit on the sample 
size ("N"), minimum number of events required before rejecting the null hypothesis ("M"), critical value ("cv") and a relative risk ("RR").
The metrics calculated are the statistical power, the expected time to signal when the null hypothesis is 
rejected, and the expected sample size at the end of the analysis whether the null hypothesis was rejected or not.}
\usage{
Performance.Binomial(N,M=1,cv,z="n",p="n",RR=2)
}
\arguments{
\item{N}{The upper limit on the sample size (length of surveillance) expressed in terms of the total number of events. 
"N" must be greater than 0. To avoid very large computation times, we suggest not using values greater than 1000. Typically,
this is not a major restriction. For example, for "RR=1.1" and "alpha=0.01" and "z=1", the statistical power is approximately 1 for "N>500". There is no default value.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. The default value is 'M=1', 
which means that even a single event can reject the null hypothesis if it occurs sufficiently early. A good rule of thumb is to set 'M=4'
(Kulldorff and Silva, 2015). It must be a positive integer.}
\item{cv}{Critical value, defining the upper rejection boundary. The null hypothesis is rejected when the log-likelihood value is greater than "cv". The "cv" parameter is usually obtained by first running \code{CV.G.Binomial}.It must be a positive number, and there is no default.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case under the null hypothesis. There is no default value.}
\item{p}{The probability of having a case under the null hypothesis. There is no default value.}
\item{RR}{The relative risk (>=1) for which statistical power, expected signal time and expected length of surveillance are calculated. The default is RR=2.}
}
\details{
For the continuous Binomial MaxSPRT, the function \code{Performance.Binomial} calculates the statistical power, the expected time to signal 
when the null hypothesis is rejected, and the expected sample size until the analysis ends whether the null is rejected or not. When the null hypothesis is true, the probability of
having a case, instead of a control, is \eqn{p=1/(1+z)}. But, if the null hypothesis is false, and the true relative risk is a value 'RR>1', then the
probability of having a case is \eqn{p=RR/(RR+z)}.
If the user wants to calculate the exact Type I error probability for a given "cv", that can be done by setting "RR=1", in which case the power output value is the exact size of the test.

The input z represents the number of controls matched to each case. For example, if there are 3 controls matched to each case, "z=3". 
In a self-control analysis, z is the ratio of the control interval to the risk interval.
For example, if the risk interval is 2 days long and the control interval is 7 days long, z=7/2.
In terms of p, the binomial probability under the null hypothesis, p=1/(1+z), or equivalently, z=1/p-1. The parameter z must be a positive number.

Alternatively, instead of z the user can specify p directly.
Note that only one of these inputs, z or p, has to be specified, but if both are entered the code will only work if z and p are such that p=1/(1+z).
Otherwise, an error message will appear to remind that such condition must be complied.

}
\value{
\item{Power}{The statistical power.}
\item{ESignaltime}{The expected time to signal given that the null hypothesis is rejected.}
\item{EsampleSize}{The expected sample size when the analysis ends (length of surveillance) whether the null hypothesis was rejected or not.}
}
\author{Ivair Ramos Silva and Martin Kulldorff.
} 
\references{
Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

Kulldorff M, Silva IR. (2015). Continuous post-market sequential safety surveillance with minimum events to signal. arxiv:1503.01978 [stat.ap].
}
\section{Acknowledgements}{
Development of the Performance.Binomial function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999.
}

\section{See also}{
\code{\link[Sequential]{CV.Binomial}}: for calculating critical values in continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{SampleSize.Binomial}}: for calculating the minimum sample size given a target power in continuous sequential analysis with binomial data.\cr 
\code{\link[Sequential]{Performance.G.Binomial}}: for calculating the statistical power, expected time to signal and expected sample size for group sequential analysis with binomial data.\cr
}

\examples{

result<- Performance.Binomial(N=30,M=1,cv=2,z=1,RR=2)
# if you type:
result
# then you will get the following output:
# $power
# [1] 0.658732

# $signaltime
# [1] 10.7893

# $surveillancetime
# [1] 17.3453

}

