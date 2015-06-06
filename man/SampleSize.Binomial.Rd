\name{SampleSize.Binomial}
\alias{SampleSize.Binomial}
\title{Sample Size Calculation for Continuous Sequential Analysis with Binomial Data.}
\description{The function \code{SampleSize.Binomial} obtains the sample size ("N") needed to guarantee a desired statistical power, for a fixed true relative risk,
when applying the MaxSPRT test with a Wald-type upper boundary, which is flat with respect to the log-likelihood ratio.}
\usage{
SampleSize.Binomial(RR,alpha=0.05,power=0.9,M=1,z=1)
      }
\arguments{
\item{RR}{A target relative risk to be detected with the requested statistical power.}
\item{alpha}{The significance level. The default value is 'alpha=0.05'. To avoid numerical complications when alpha is close to one, 'alpha' must be in the range (0,0.5].}
\item{power}{The target overall statistical power to detect an increased risk of the relative risk (RR). The default value is "power=0.90".}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. It must be a positive integer. The default value is 'M=1.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case. For example, if there are 3 controls matched to each case, z=3.  In a self-control analysis, z is the ratio of the control interval to the risk interval. For example, if the risk interval is 2 days long and the control interval is 7 days long, z=7/2. In terms of p, the binomial probability under the null hypothesis, p=1/(1+z), or equivalently, z=1/p-1. The parameter z must be a positive number. The default value is z=1.}
}
\details{
The function \code{SampleSize.Binomial} calculates the sample size ("N") to be used for the
continuous binomial MaxSPRT in order to provide the desired statistical power given a true "RR". 
The solution is exact using iterative numerical calculations. For details on how to abtain critical values under the binomial MaxSPRT (Kulldorff et al., 2011).}
\value{
\item{Required_N}{The lower limit on the sample size ("N") needed to achive the desired statistical "power" given the relative risk "RR".}
\item{cv}{Critical value for the sequential binomial serveillance for the overall significance level "alpha".}
\item{Type_I_Error}{Exact size of the test by fixing "alpha", "N=Required_N" and "cv=Associated_cv".}
\item{Actual_power}{The exact power promoted by the "Required_N".}
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


