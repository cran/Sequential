\name{SampleSize.Poisson}
\alias{SampleSize.Poisson}
\title{Sample size calculation for continuous sequential analysis with Poisson data.}
\description{The function \code{SampleSize.Poisson} obtains the required sample size (length of surveillance) needed to guarantee a desired statistical
power for a pre-specified relative risk, when doing continuous sequential analysis for Poisson data with a flat upper boundary in the scale of the Wald
type MaxSPRT (log-likelihood ratio scale), Pocock, OBrien-Fleming, or Wang-Tsiatis scales. Alternatively, \code{SampleSize.Poisson} calculates sample
sizes for non-flat signaling thresholds for user-defined alpha spending functions. It can also be used to approximate the sample size needed when doing
group sequential analysis for Poisson data.
}
\usage{
SampleSize.Poisson(alpha=0.05,power=0.9,M=1,D=0,RR=2,
precision=0.001,alphaSpend="n",rho="n",gamma="n",
Statistic=c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis"),
Delta="n",Tailed="upper")
      }
\arguments{
\item{alpha}{The significance level. The default value is alpha=0.05. Must be in the range (0,0.5].}
\item{power}{The target vector of overall statistical powers to detect an increased relative risk (RR). The default value is power=0.90.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. It must be a positive integer.
A good rule of thumb is to set M=4 (Kulldorff and Silva, 2015). The default value is M=1, which means that even a single event
can reject the null hypothesis if it occurs sufficiently early.}
\item{D}{The expected number of events under the null hypothesis at the first look at the data.
This is used when there is an initial large chunk of data arriving, followed by continuous sequential analysis. The default value is D=0, which is
also the best choice. This means that there is no delay in the start of the sequential analyses. If D is very large, the maximum 
sample size will be set equal to D if a non-sequential analysis provides the desired power.}
\item{RR}{The target vector of relative risks to be detected with the requested statistical vector of powers. The default value is RR=2.}
\item{precision}{The tolerance for the difference between the requested and actual statistical power. Should be very small. The default value is precision=0.001.}
\item{alphaSpend}{An integer between 1 to 4. Default is "n". See Details.}
\item{rho}{Positive number used for the power-type alpha spending function (\code{AlphaSpend=1}) only. The default value is "rho=0.5". See Details.}
\item{gamma}{Positive number used for the gamma-type alpha spending function (\code{AlphaSpend=4}) only. There is no default value. See Details.}
\item{Statistic}{The test statistic scale to deliver the signaling threshold. See Details.}
\item{Delta}{Parameter needed for calculation of Wang-Tsiatis test statistic if this is the option selected in "Statistic". Must be a number in the (0, 0.5] interval. There is no default value.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
}
\details{
When using the MaxSPRT and the \code{CV.Poisson} function to conduct continuous sequential analysis for Poisson data, the null
hypothesis is rejected when the log likelihood ratio exceeds the pre-determined critical value calculated by \code{CV.Poisson}.
The sequential analysis ends without rejecting the null hypothesis when a predetermined upper limit on the sample size is
reached, expressed in terms of the expected number of events under the null hypothesis. For example, the sequential analysis
may end as soon as the sample size is such that there are 50 expected events under the null.

The default in the function \code{SampleSize.Poisson} is for calculating the upper limit on the sample size (length of surveillance) required
for the continuous Poisson based MaxSPRT (alphaSpend="n") to achieve the desired statistical power for a pre-specified relative risk RR.
The solution is exact using iterative numerical calculations (Kulldorff et al., (2011).

While designed for continuous sequential analysis with flat threshold in the scale of the MaxSPRT statistic, the \code{SampleSize.Poisson} function can also be used to approximate the
required upper limit on the sample size that is needed when doing group sequential analysis for Poisson data, using the \code{CV.G.Poisson function}.

There is also the possibility to calculate the sample size for an user-defined alpha spending plan. This is possible with the input parameter \code{alphaSpend}.
The user can select among one of the four classical alpha spending shapes bellow:\cr
\eqn{F_{1}(t)=\alpha t^{\rho}}, where \eqn{\rho>0},\cr 
\eqn{F_{2}(t)=2-2\Phi(x_{\alpha}\sqrt{t^{-1}})}, where \eqn{x_{\alpha}=\Phi^{-1}(1-\alpha/2)},\cr 
\eqn{F_{3}(t)= \alpha \times log(1+[exp{1}-1]\times t) },\cr
\eqn{F_{4}(t)=\alpha[1-exp(-t\gamma)]/[1-exp(-\gamma)]} with \eqn{\gamma \in \Re},\cr
and \eqn{t} represents a fraction of the maximum length of surveillance. 

To select one of the four alpha spending types above, and using an integer \eqn{i} to indicate the type among
\eqn{i=} 1, 2, 3, and 4, for \eqn{F_{1}(t)}, \eqn{F_{2}(t)}, \eqn{F_{3}(t)} and \eqn{F_{4}(t)}, respectively,
one needs to set \code{alphaSpend=i}. Specifically for \code{alphaSpend=1}, it is necessary to choose a \code{rho} value,
or a \code{gamma} value if \code{alphaSpend=4} is used.

For more details on these alpha spending choices,
see the paper by Silva et al. (2021), Section 2.7.

When one sets \code{alphaSpend=i}, the threshold impplied by the correspondent alpha spending is calculated.
The function delivers the threshold in the scale of a test statistic selected by the user with the input
\code{Statistic} among the classic methods:
MaxSPRT (Kulldorf et al., 2011), Pocock (Pocock, 1977), OBrien-Fleming (O'Brien and Fleming, 1979), or Wang-Tsiatis (Jennison and Turnbull, 2000). 
For \code{Statistic="Wang-Tsiatis"}, the user has to choose a number in the (0, 0.5] interval for \code{Delta}.

}

\value{
\item{SampleSize_by_RR_Power}{A table containing the main performance measures associated to the required samples sizes, expressed in terms of the expected number of events under the null hypothesis, for each combination of RR and power.}
}

\author{ Ivair Ramos Silva, Martin Kulldorff.}
\references{
Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. Sequential Analysis, 30: 58--78.
Kulldorff M, Silva IR. (2015). Continuous Post-market Sequential Safety Surveillance with Minimum Events to Signal. REVSTAT Statistical Journal, 15(3): 373--394.
Silva IR, Maro J, Kulldorff M. (2021). Exact sequential test for clinical trials and post-market drug and vaccine safety surveillance with Poisson and binary data. Statistics in Medicine, DOI: 10.1002/sim.9094.
}
\section{Acknowledgements}{
Development of the \code{SampleSize.Poisson} function was funded by:\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil (v1.0).\cr
-	Bank for Development of the Minas Gerais State (BDMG), Brazil (v1.0).\cr
-     National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0.1,2.0.2). 
}

\keyword{Continuous MaxSPRT analysis}

\section{See also}{ 
\code{\link[Sequential]{CV.Poisson}}: for calculating critical values for continuous sequential analysis with Poisson data.\cr
\code{\link[Sequential]{Performance.Poisson}}: for calculating the statistical power, expected time to signal and expected sample size for continuous sequential analysis with Poisson data\cr
\code{\link[Sequential]{SampleSize.Binomial}}: for calculating the minimum sample size given a target power in continuous sequential analysis with binomial data. 
}


\examples{

### Example 1:   
##  Sample size required to obtain a power of 80%, for a relati-
##  ve risk of 3, no delay for starting the surveillance (D=0),
##  and when the null hypothesis can be rejected with one event
##  (M=1) under an alpha level of 5%.

# result1<- SampleSize.Poisson(alpha=0.05,power=0.8,M=1,D=0,RR=3)
# result1


## Example 2:
##  Sample size required to obtain a power of 90%, for a relati-
##  ve risk of 2, no delay for starting the surveillance (D=0),
##  and when the null hypothesis can be rejected only after 2
##  events (M=2) under an alpha level of 10%.
##
# result2<- SampleSize.Poisson(alpha=0.1,power=0.9,M=2,D=0,RR=2)
# result2

## Example 3:
##  Sample size calculated for the non-flat threshold using the
##  power-type alpha spending (alphaSpend=1), with rho=1,
##  to obtain a power of 80% for a relative risk of 2.5, delay to
##  start the surveillance equal to 2 (D=1), and the null hypo-
##  thesis can be rejected with 3 events (M=3) under an alpha
##  level of 5%. The critical values will be shown in the scale
##  of the MaxSPRT statistic.

# result3<- SampleSize.Poisson(alpha=0.05,power=0.8,M=3,D=1,RR=2.5,
# alphaSpend=1,rho=1,Statistic="MaxSPRT")
# result3

}


