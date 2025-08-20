\name{Performance.Threshold.CondPoisson}
\alias{Performance.Threshold.CondPoisson}
\title{Performance and alpha spending for user-defined signaling threshold for sequential analysis with conditional Poisson data.}
\description{The function \code{Performance.Threshold.CondPoisson} calculates the statistical power, expected time to signal, expected sample size and alpha spending associated to any user-specified signaling threshold, flat or non-flat, for continuous or group sequential analysis with conditional Poisson data,
all for a pre-specified upper limit on the sample size.}
\usage{
Performance.Threshold.CondPoisson(K,cc,CV.upper="n",
Person_timeRatioH0="n",GroupSizes="n",Tailed="upper",RR)
      }
\arguments{
\item{K}{The upper limit on the sample size (length of surveillance) expressed in terms of the number of events arriving in the surveillance period. There is no default value.}
\item{cc}{Number of events observed in the historical period. There is no default value.}
\item{CV.upper}{User-specified signaling threshold given in the scale of the CMaxSPRT test statistic. There is no default value.}
\item{Person_timeRatioH0}{Test-specific amount of information of each test given in terms of the ratio between the person-time in the surveillance period and the overall person-time of the historical period. See Details. There is no default value.}
\item{GroupSizes}{Test-specific number of events between two looks at the data for group or continuous sequential analysis. There is no default value. See Details.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
\item{RR}{Vector of relative risks for performance calculation. There is no default value.}
}
\details{
For continuous and group sequential analysis based on monitoring Poisson data conditioned on matched
historical Poisson data, the power, expected time to signal, expected sample size and alpha spending impplied by user-specified thresholds are calculated with
\code{Performance.Threshold.CondPoisson}. The user can select one between two scales to enter with the
threshold, the Conditional Maximized Sequential
Probability Ratio Test statistic (CMaxSPRT) scale (Li and Kulldorff, 2010), or the surveillance versus historical person-time ratio (Silva et al., 2019a).
For the CMaxSRT scale, the input is \code{CV.upper}. This can be entered as a vector for group sequential analysis. For example, for a three-group sequential test,
the i-th entry represents the critical value for the i-th test, with i=1, 2, 3. If a single number is informed in \code{CV.upper}, then
a flat critical value for all tests is used for both continuous or group sequential fashions. The number of tests is defined with the input \code{GroupSizes},
as shall be detailed here after the desciption of \code{Person_timeRatioH0}. 

An alternative way to inform the threshold is by using \code{Person_timeRatioH0}, which is in the scale of the ratio between the person-time
in the surveillance period and the
overall person-time of the historical period. Using the notation by Silva et al. (2019a) and Silva et al. (2019b), let \eqn{V} denote the
total person-time from the historical data, where \code{cc} events were observed, and let \eqn{P_{k(i)}} denote the the cummulative
person-time from the surveillance data at the i-th test with a cummulative \eqn{k(i)} events. With \code{Person_timeRatioH0},
the entries must have increasing numbers, from the first to the last. For example,
for a three-group sequential plan with sample sizes of 20, 15, 25, a hypothetical choice is
\code{Person_timeRatioH0=c(0.1, 0.5, 1)}. This way, H0 is rejected if: \eqn{P_20/V} <= 0.1 in the first test,
or \eqn{P_35/V} <= 0.5 in the second test, or \eqn{P_60/V} <= 1 in the third test. 

Note: only one of the inputs \code{CV.upper} or \code{Person_timeRatioH0} is to be used.
     
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
\item{Performance}{A matrix with the following three performance measures for each target RR: statistical power, expected time to signal and expected sample size.}
\item{AlphaSpend}{The alpha spending associated to the user-specified threshold.}
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

Silva IR, Maro J, Kulldorff M. (2019). Exact Sequential Analysis Using R Sequential. Working Paper.
}


\section{Acknowledgements}{
Development of the Performance.Threshold.CondPoisson function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0,2.0 to 3.1).\cr
-	Federal University of Ouro Preto (UFOP), through contract under internal UFOP's resolution CEPE 4600 (v2.0 to 3.1).\cr
\cr
}

\keyword{Performance and Alpha spending given threshold with conditional Poisson data}

\section{See also}{
\code{\link[Sequential]{Performance.AlphaSpend.CondPoisson}}: for calculating performance and signaling threshold for user-specified alpha spending with conditional Poisson data.\cr
\code{\link[Sequential]{CV.CondPoisson}}: for calculating Wald-type signaling thresholds for continuous sequential analysis with conditional Poisson data based on the CMaxSPRT test statistic.\cr
\code{\link[Sequential]{Analyze.CondPoisson}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion for condicional Poisson data based on the CMaxSPRT test statistic.
}


\examples{

#### Example 1
## Power, expected time to signal, expected sample size and
## alpha spending of three group CMaxSPRT sequential analysis with
#  a maximum sample size of 30 events for upper-tailed
#  testing, i.e. H0:RR<=1, with regular groups of sizes 10
#  and a flat threshold equal to 3.6. 
#  The statistical performance is evaluated for two
#  target RR= 1.5, 2:

# Performance.Threshold.CondPoisson(K=30,cc=10,CV.upper=3.6, Person_timeRatioH0="n",
# GroupSizes=10,RR=c(1.5,2))

#### Example 2
## Power, expected time to signal, expected sample size and
#  alpha spending of three group CMaxSPRT sequential analysis with
#  a maximum sample size of 30 events for upper-tailed
#  testing, i.e. H0:RR<=1, with regular groups of sizes 10
#  and thresholds given in the P_k/V scale:
#  "Person_timeRatioH0=c(0.1, 0.5, 1)". 
#  The statistical performance is evaluated for two
#  target RR= 1.5, 2:

# Performance.Threshold.CondPoisson(K=30,cc=10,CV.upper="n", Person_timeRatioH0=c(0.1, 0.5, 1),
# GroupSizes=10,RR=c(1.5,2))

}


