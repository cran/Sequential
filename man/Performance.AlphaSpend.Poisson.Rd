\name{Performance.AlphaSpend.Poisson}
\alias{Performance.AlphaSpend.Poisson}
\title{Calculates performance and signaling thresholds for user-defined alpha spending for sequential analysis with Poisson data.}
\description{The function \code{Performance.AlphaSpend.Poisson} calculates power, expected time to signal, expected sample size and the signaling threshold (critical values) associated to any user-specified alpha spending for continuous sequential analysis with Poisson data,
all for a pre-specified upper limit on the sample size.}
\usage{
Performance.AlphaSpend.Poisson(SampleSize,alpha=0.05,D=0,M=1,RR,alphaSpend=1,rho=0.5,
gamma="n",Statistic=c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis"),
Delta="n",Tailed="upper")
      }
\arguments{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) expressed in terms of the expected number of events under the null hypothesis. There is no default value.}
\item{alpha}{The overall significance level.}
\item{D}{The expected number of events under the null hypothesis at the first look at the data.
This is used when there is an initial large chunk of data arriving, followed by continuous sequential analysis. The default value is D=0, which is
also the best choice. This means that there is no delay in the start of the sequential analyses. If D is very large, the maximum 
sample size will be set equal to D if a non-sequential analysis provides the desired power.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. It must be a positive integer.
A good rule of thumb is to set M=4 (Kulldorff and Silva, 2015). The default value is M=1, which means that even a single event
can reject the null hypothesis if it occurs sufficiently early.}
\item{RR}{Vector of relative risks for performance calculation. There is no default value.}
\item{alphaSpend}{A vector with the cummulative Type I error probabiliy to be spent up to each test. Alternatively, one can use an integer between 1 to 4. Default is 1. See Details.}
\item{rho}{Positive number used for the power-type alpha spending function (\code{AlphaSpend=1}) only. The default value is "rho=0.5". See Details.}
\item{gamma}{Positive number used for the gamma-type alpha spending function (\code{AlphaSpend=4}) only. There is no default value. See Details.}
\item{Statistic}{The test statistic scale to deliver the signaling threshold. See Details.}
\item{Delta}{Parameter needed for calculation of Wang-Tsiatis test statistic if this is the option selected in "Statistic". Must be a number in the (0, 0.5] interval. There is no default value.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
}
\details{
For continuous and group sequential analysis based on monitoring Poisson data, the threshold impplied by user-specified alpha spending is calculated with
\code{Performance.AlphaSpend.Poisson}. The function delivers the threshold in the scale of a test statistic selected by the user with the input
\code{Statistic} among the classic methods:
MaxSPRT (Kulldorf et al., 2011), Pocock (Pocock, 1977), OBrien-Fleming (O'Brien and Fleming, 1979), or Wang-Tsiatis (Jennison and Turnbull, 2000). 
For \code{Statistic="Wang-Tsiatis"}, the user has to choose a number in the (0, 0.5] interval for \code{Delta}.


\code{alphaSpend} is used for arbitrary cumulative type I error probability spending defined by the user.
Alternatively, the user can select among one of the four classical alpha spending shapes bellow:\cr
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

For \code{RR} the user must specify the target relative risks for calculation of statistical performance measures.
It can be a vector of positive numbers or a single number. 

For details on the calculation of signaling thresholds and alpha spending for Poisson data,
see the paper by Silva et al. (2018).


}
\value{
\item{cvs}{Signaling threshold in the scale of the \code{Statistic} selected for the user-specified alpha spending.}
\item{Performance}{A matrix with the following three performance measures for each target RR: statistical power, expected time to signal and expected sample size.}
}
\author{Ivair Ramos Silva, Martin Kulldorff.
}
\references{
Jennison C, Turnbull B. (2000). Group Sequential Methods with Applications to Clinical Trials, London: Chapman and Hall/CRC.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

Li L, Kulldorff M. (2010). A conditional maximized sequential probability ratio test for Pharmacovigilance. Statistics in Medicine, 29 (2), 284--295. 

O'Brien PC, Fleming TR. (1979). A multiple testing procedure for clinical trials. Biometrics. 35:549--556.

Pocock SJ. (1977). Group sequential methods in the design and analysis of clinical trials. Biometrika. 64:191--199.

Silva IR. (2018). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance With Poisson Data. Methodol Comput Appl Probab, 20(2), 739--750.

Silva IR, Kulldorff M. (2015). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva IR, Maro J, Kulldorff M. (2021). Exact sequential test for clinical trials and post-market drug and vaccine safety surveillance with Poisson and binary data. Statistics in Medicine, DOI: 10.1002/sim.9094.
}


\section{Acknowledgements}{
Development of the Performance.AlphaSpend.Poisson function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0,2.0 to 3.1).\cr
-	Federal University of Ouro Preto (UFOP), through contract under internal UFOP's resolution CEPE 4600 (v2.0 to 3.1).\cr
\cr
}

\keyword{Performance and threshold for given alpha spending with Poisson data}

\section{See also}{
\code{\link[Sequential]{Performance.Threshold.Poisson}}: for calculating performance and alpha spending for user-specified signaling threshold with Poisson data.\cr
\code{\link[Sequential]{CV.Poisson}}: for calculating Wald-type signaling thresholds for continuous sequential analysis with Poisson data based on the CMaxSPRT test statistic.\cr
\code{\link[Sequential]{Analyze.Poisson}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion for Poisson data based on the CMaxSPRT test statistic.
}


\examples{

#### Example
## Perfofmance and threshold for continuous sequential analysis
#  with a maximum sample size of 30 events for upper-tailed
#  testing, i.e. H0:RR<=1, with alpha spending of the
#  power-type and threshold delivered in the MaxSPRT
#  test statistic scale.  
#  The statistical performance is evaluated for three
#  target RR= 1.5:

# res<- Performance.AlphaSpend.Poisson(SampleSize=30, alpha=0.05,alphaSpend=1,
# RR=1.5,rho=0.5,gamma="n",Statistic="MaxSPRT")

}


