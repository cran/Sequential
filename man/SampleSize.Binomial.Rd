\name{SampleSize.Binomial}
\alias{SampleSize.Binomial}
\title{Sample Size Calculation for Continuous Sequential Analysis with Binomial Data.}
\description{The function \code{SampleSize.Binomial} obtains the sample size needed to guarantee a desired statistical power, for a fixed true relative risk, when doing continuous sequential analysis for binomial data with a Wald-type upper boundary, which is flat with respect to the log-likelihood ratio. It can also be used to approximate the sample size needed when doing group sequential analysis for binomial data.}
\usage{
SampleSize.Binomial(RR,alpha=0.05,power=0.9,M=1,z=1)
      }
\arguments{
\item{RR}{A target relative risk to be detected with the requested statistical power.}
\item{alpha}{The significance level. The default value is "alpha=0.05". Must be in the range (0, 0.5].}
\item{power}{The target overall statistical power to detect an increased risk of the relative risk (RR). The default value is "power=0.90".}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. It must be a positive integer. The default value is "M=1".}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case. For example, if there are 3 controls matched to each case, z=3.  In a self-control analysis, z is the ratio of the control interval to the risk interval. For example, if the risk interval is 2 days long and the control interval is 7 days long, z=7/2. In terms of p, the binomial probability under the null hypothesis, p=1/(1+z), or equivalently, z=1/p-1. The parameter z must be a positive number. The default value is z=1.}
}
\details{
The function \code{SampleSize.Binomial} calculates the sample size "N", to be used for the continuous binomial MaxSPRT in order to provide the desired statistical power for a user-specified relative risk "RR". The required sample size ("Required_N") is expressed in terms of the total number of observations, and it is the number of observations by which the sequential analysis will end without rejected the null hypothesis. The solution is  exact using iterative numerical calculations. 

The required sample size, N, increases for increasing values of "power", while N decreases for increasing values of "alpha", the relative risk "RR" and the minimum number of events needed to signal "M". For increasing values of "z", the required sample size N can either decrease or increase.

While this function calculates the required sample size for continuous sequential analysis, it can also be used as an approximation for group sequential analyses. With the same "Required_N", and all other parameters being the same, a group sequential analysis will always give higher statistical power than a continuous sequential analysis, so \code{SampleSize.Binomial} can be use to ensure the required statistical power for group sequential analyses. 
}
\value{
\item{Required_N}{The minimum sample size, N, that is required to achieve the requested statistical "power" for the relative risk "RR" provided. The sample size N is provided in terms of the total number of observations. In the case-control setting, this is equal to the total number of cases and controls. In the self-control setting, it is equal to the total number of events in either the risk or the control interval.}
\item{cv}{Critical value for the sequential binomial surveillance for the overall significance level "alpha".}
\item{Type_I_Error}{Exact size of the test by fixing "alpha", "N=Required_N" and "cv=Associated_cv". Always less than or equal to the requested "alpha", it is typically smaller due to the discrete nature of the binomial probability function.}
\item{Actual_power}{The exact actual power provided with a sample size of "Required_N". Always greater than or equal to the requested "power", it is typically greater due to the discrete nature of the binomial probability function.}
}
\author{ Ivair Ramos Silva, Martin Kulldorff.
}
\references{
Kulldorff M, Silva IR. (2015). Continuous post-market sequential safety surveillance with minimum events to signal. arxiv:1503.01978 [stat.ap].
}
\section{Acknowledgements}{
Development of the SampleSize.Binomial function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999.
}

\keyword{Continuous MaxSPRT analysis}

\section{See also}{
\code{\link[Sequential]{CV.Binomial}}: for calculating critical values in continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{Performance.Binomial}}: for calculating the statistical power, expected time to signal and expected time of analysis for continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{SampleSize.Poisson}}: sample size calculation for continuous sequential analysis with Poisson data.  
}


\examples{

result<- SampleSize.Binomial(RR=5,alpha=0.01,power=0.88,M=1,z=2)
# if you type:
result
# then you will get the following output:
# $Required_N
# [1] 25

# $cv
# [1] 4.59581

# $Type_I_Error
# [1] 0.009755004

# $Actual_power
# [1] 0.8855869

}


