\name{Performance.AlphaSpend.CondPoisson}
\alias{Performance.AlphaSpend.CondPoisson}
\title{Calculates performance and signaling thresholds for user-defined alpha spending for sequential analysis with conditional Poisson data.}
\description{The function \code{Performance.AlphaSpend.CondPoisson} calculates power, expected time to signal and signaling threshold (critical values) associated to any user-specified alpha spending, flat or non-flat, for continuous or group sequential analysis with conditional Poisson data,
all for a pre-specified upper limit on the sample size.}
\usage{
Performance.AlphaSpend.CondPoisson(K,cc,alpha,AlphaSpend="n",
GroupSizes="n",rho=0.5,gamma="n",Tailed="upper",RR)
      }
\arguments{
\item{K}{The upper limit on the sample size (length of surveillance) expressed in terms of the number of events arriving in the surveillance period. There is no default value.}
\item{cc}{Number of events observed in the historical period. There is no default value.}
\item{alpha}{The overall significance level.}
\item{AlphaSpend}{A vector with the cummulative Type I error probabiliy to be spent up to each test. Alternatively, one can use an integer between 1 to 4. There is no default value. See Details.}
\item{GroupSizes}{Test-specific number of events between two looks at the data for group or continuous sequential analysis. There is no default value. See Details.}
\item{rho}{Positive number used for the power-type alpha spending function (\code{AlphaSpend=1}) only. The default value is "rho=0.5". See Details.}
\item{gamma}{Positive number used for the gamma-type alpha spending function (\code{AlphaSpend=4}) only. There is no default value. See Details.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
\item{RR}{Vector of relative risks for performance calculation. There is no default value.}
}
\details{
For continuous and group sequential analysis based on monitoring Poisson data conditioned on matched
historical Poisson data, the threshold impplied by user-specified alpha spending is calculated with
\code{Performance.AlphaSpend.CondPoisson}. The function delivers the threshold in two scales, the the Conditional Maximized Sequential
Probability Ratio Test statistic (CMaxSPRT) scale (Li and Kulldorff, 2010), and the surveillance versus historical person-time ratio (Silva et al., 2019a).
For the CMaxSPRT scale, the null hypothesis is rejected if the test statistic exceeds the critical value at some test. Regarding the person-time ratio scale,
using the notation by Silva et al. (2019a) and Silva et al. (2019b), let \eqn{V} denote the
total person-time from the historical data, where \code{cc} events were observed, and let \eqn{P_{k(i)}} denote the the cummulative
person-time from the surveillance data at the i-th test with a cummulative \eqn{k(i)} events. Suppose that for
a three-group sequential plan, with sample sizes of 20, 15, 25, the critival values were 0.1, 0.5, and 1. This way, H0 is rejected if: \eqn{P_20/V} <= 0.1 in the first test,
or \eqn{P_35/V} <= 0.5 in the second test, or \eqn{P_60/V} <= 1 in the third test. 

\code{AlphaSpend} is used for arbitrary cumulative type I error probability spending defined by the user.
Alternatively, the user can select among one of the four classical alpha spending shapes bellow:\cr
\eqn{F_{1}(t)=\alpha t^{\rho}}, where \eqn{\rho>0},\cr 
\eqn{F_{2}(t)=2-2\Phi(x_{\alpha}\sqrt{t^{-1}})}, where \eqn{x_{\alpha}=\Phi^{-1}(1-\alpha/2)},\cr 
\eqn{F_{3}(t)= \alpha \times log(1+[exp{1}-1]\times t) },\cr
\eqn{F_{4}(t)=\alpha[1-exp(-t\gamma)]/[1-exp(-\gamma)]} with \eqn{\gamma \in \Re},\cr
and \eqn{t} represents a fraction of the maximum length of surveillance. For more details on these alpha spending choices,
see the paper by Silva et al. (2019c), Section 2.6.


To select one of the four alpha spending types above, and using an integer \eqn{i} to indicate the type among
\eqn{i=} 1, 2, 3, and 4, for \eqn{F_{1}(t)}, \eqn{F_{2}(t)}, \eqn{F_{3}(t)} and \eqn{F_{4}(t)}, respectively,
one needs to set \code{AlphaSpend=i}. Specifically for \code{AlphaSpend=1}, it is necessary to choose a \code{rho} value,
or a \code{gamma} value if  \code{AlphaSpend=4} is used.


With \code{GroupSizes} the user informs the sample size of each test in the scale of the
number of events in the surveillance period. Therefore, only positive
numbers are accepted in \code{GroupSizes}. For irregular group sizes, a vector must be informed
with each test-specific number of events between two looks at the data, therefore the entries of
\code{GroupSizes} must sums up \code{K}. For regular group sizes, a single number can be informed
for the constant sample size of each test. For example, for continuous sequential analysis,
\code{GroupSizes=1}. For ten-group sequential analysis with \code{K=50}, \code{GroupSizes=5}. 

For \code{RR} the user must specify the target relative risks for calculation of statistical performance measures.
It can be a vector of positive numbers or a single number. 

For details on the calculation of signaling thresholds and alpha spending for Poisson data conditioned to historical data,
see the papers by Silva et al. (2019a) and Silva et al. (2019b), respectively.


}
\value{
\item{CV}{Signaling threshold in the scale of CMaxSPRT associated to the user-specified alpha spending.}
\item{Person-timeRatioH0}{Signaling threshold in the \eqn{P_k/V} scale.}
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

Silva IR, Li L, Kulldorff M. (2019a), Exact conditional maximized sequential probability ratio test adjusted for covariates. Sequential Analysis, 38(1), 115--133.

Silva IR., Lopes LM., Dias P., Yih WK. (2019b). Alpha Spending for Historical Versus Surveillance Poisson Data With CMaxSPRT. Statistics in Medicine, 38(12), 2126--2138.

Silva IR, Kulldorff M. (2015). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva IR, Maro J, Kulldorff M. (2019c). Exact Sequential Analysis Using R Sequential. Working Paper.
}


\section{Acknowledgements}{
Development of the Performance.AlphaSpend.CondPoisson function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0,2.0 to 3.1).\cr
-	Federal University of Ouro Preto (UFOP), through contract under internal UFOP's resolution CEPE 4600 (v2.0 to 3.1).\cr
\cr
}

\keyword{Performance and Threshold for given alpha spending with conditional Poisson data}

\section{See also}{
\code{\link[Sequential]{Performance.Threshold.CondPoisson}}: for calculating performance and alpha spending for user-specified signaling threshold with conditional Poisson data.\cr
\code{\link[Sequential]{CV.CondPoisson}}: for calculating Wald-type signaling thresholds for continuous sequential analysis with conditional Poisson data based on the CMaxSPRT test statistic.\cr
\code{\link[Sequential]{Analyze.CondPoisson}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion for condicional Poisson data based on the CMaxSPRT test statistic.
}


\examples{

#### Example 1
## Performance and threshold for three group CMaxSPRT sequential with
#  a maximum sample size of 30 events for upper-tailed
#  testing, i.e. H0:RR<=1, with regular groups of sizes 10
#  and alpha spending of 0.01, 0.02, and 0.05 for tests
#  1, 2 and 3, respectively. 
#  The statistical performance is evaluated for three
#  target RR= 1.2, 1.5, 2:

# res<- Performance.AlphaSpend.CondPoisson(K=30,cc=10,alpha=0.05,
# AlphaSpend=c(0.01,0.02,0.05),GroupSizes=c(10,10,10),RR=c(1.2,1.5,2))

#### Example 2
## Performance and threshold for three group CMaxSPRT with
#  a maximum sample size of 30 events for upper-tailed
#  testing, i.e. H0:RR<=1, with regular groups of sizes 10
#  and alpha spending of the power-type.
#  The statistical performance is evaluated for three
#  target RR= 1.2, 1.5, 2:

# res<- Performance.AlphaSpend.CondPoisson(K=30,cc=10,alpha=0.05,
# AlphaSpend=1,GroupSizes=c(10,10,10),rho=0.75,RR=c(1.2,1.5,2))

}


