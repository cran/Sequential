\name{Performance.AlphaSpend.Binomial}
\alias{Performance.AlphaSpend.Binomial}
\title{Calculates performance and signaling threshold for user-defined alpha spending for sequential analysis with binomial data.}
\description{The function \code{Performance.AlphaSpend.Binomial} calculates power, expected time to signal, expected sample size and signaling threshold (critical values) associated to any user-specified alpha spending for continuous or group sequential analysis with binomial data.
The user can select the scale for the signaling threshold among MaxSPRT, Pocock, OBrien-Fleming, or Wang-Tsiatis test statistics, all for a pre-specified upper limit on the sample size.}
\usage{
Performance.AlphaSpend.Binomial(N,alpha,AlphaSpend=1,AlphaSpendLower="n",
AlphaSpendUpper="n",z="n",p="n",GroupSizes="n",
Tailed="upper",rho=0.5,gamma="n",RR=2,Statistic=c("MaxSPRT", "Pocock",
"OBrien-Fleming", "Wang-Tsiatis"),Delta="n")
      }
\arguments{
\item{N}{The upper limit on the sample size (length of surveillance) expressed in terms of the total number of events (cases plus controls). There is no default value.}
\item{alpha}{Overall significance level. Must be a number in the (0, 0.5] interval. There is no default value.}
\item{AlphaSpend}{A vector with the cummulative Type I error probabiliy to be spent up to each test when "Tailed=lower" or "Tailed=upper". Alternatively, one can use an integer between 1 to 4. Default is 1. See Details.}
\item{AlphaSpendLower}{A vector with the cummulative Type I error probabiliy regarding the lower signaling threshold to be spent up to each test when "Tailed=two". Alternatively, one can use an integer between 1 to 4. Default is "n", which means the specified "AlphaSpend" divided by 2.}
\item{AlphaSpendUpper}{A vector with the cummulative Type I error probabiliy regarding the upper signaling threshold to be spent up to each test when "Tailed=two". Alternatively, one can use an integer between 1 to 4. Default is "n", which means the specified "AlphaSpend" divided by 2.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case under the null hypothesis. There is no default value.}
\item{p}{The probability of having a case under the null hypothesis. If just a single number is given, then it will be used as a constant probability for all groups. Otherwise, the dimension of p must coincide with the dimension of GroupSizes. There is no default value.}
\item{GroupSizes}{Vector with the total number of events (cases+controls) between two looks at the data with regular and irregular group sizes. Important: Must sums up N. For continuos sequential analysis, specify GroupSizes=1. There is no default value.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
\item{rho}{Positive number used for the power-type alpha spending function (\code{AlphaSpend=1}) only. The default value is "rho=0.5". See Details.}
\item{gamma}{Positive number used for the gamma-type alpha spending function (\code{AlphaSpend=4}) only. There is no default value. See Details.}
\item{RR}{Vector of relative risks for performance calculation. There is no default value.}
\item{Statistic}{The test statistic scale used for "cvs.lower" and "cvs.upper", therefore, this input has no effect if "cases.lower" and "cases.upper" are used instead. There is no default.}
\item{Delta}{Parameter needed for calculation of Wang-Tsiatis test statistic if this is the option selected in "Statistic". Must be a number in the (0, 0.5] interval. There is no default value.}
}
\details{
For continuous and group sequential analysis with binomial data, signaling thresholds for
user-specified alpha spending are calculated with \code{Performance.AlphaSpend.Binomial}.

\code{N} must be a positive integer defining the maximum length of surveillance. To avoid very large computation times,
we suggest not using values greater than 1000.

\code{AlphaSpend} is used for arbitrary cumulative type I error probability spending defined by the user.
Alternatively, the user can select among one of the four classical alpha spending shapes bellow:\cr
\eqn{F_{1}(t)=\alpha t^{\rho}}, where \eqn{\rho>0},\cr 
\eqn{F_{2}(t)=2-2\Phi(x_{\alpha}\sqrt{t^{-1}})}, where \eqn{x_{\alpha}=\Phi^{-1}(1-\alpha/2)},\cr 
\eqn{F_{3}(t)= \alpha \times log(1+[exp{1}-1]\times t) },\cr
\eqn{F_{4}(t)=\alpha[1-exp(-t\gamma)]/[1-exp(-\gamma)]} with \eqn{\gamma \in \Re},\cr
and \eqn{t} represents a fraction of the maximum length of surveillance. For more details on these alpha spending choices,
see the paper by Silva et al. (2019), Section 2.6.


To select one of the four alpha spending types above, and using an integer \eqn{i} to indicate the type among
\eqn{i=} 1, 2, 3, and 4, for \eqn{F_{1}(t)}, \eqn{F_{2}(t)}, \eqn{F_{3}(t)} and \eqn{F_{4}(t)}, respectively,
one needs to set \code{AlphaSpend=i}. Specifically for \code{AlphaSpend=1}, it is necessary to choose a \code{rho} value,
or a \code{gamma} value if \code{AlphaSpend=4} is used.

\code{z} is a vector of positive numbers representing the matching ratios for each test (group). If a single number is given, then it will be used as a constant
matching ratio for all tests (groups). Otherwise, the dimension of \code{z} must coincide with the dimension of \code{GroupSizes}.
\code{z} represents the number of controls matched to each case. For example, if there are 3 controls matched to each case, \code{z}=3. 
In a self-control analysis, z is the ratio of the control interval to the risk interval.
For example, if the risk interval is 2 days long and the control interval is 7 days long, \code{z}=7/2.
In terms of \code{p}, the binomial probability under the null hypothesis, \code{p=1/(1+z)}, or equivalently, \code{z=1/p-1}.

Alternatively, instead of \code{z} the user can specify \code{p} directly.
Note that only one of these inputs, \code{z} or \code{p}, has to be specified, but if both are entered the code will only work if \code{z}
and \code{p} are such that \code{p=1/(1+z)}. Otherwise, an error message will appear to remind that such condition must be complied.
   
For details about the algorithm used to calculate the critical value, see the paper by Silva (2018).

With \code{GroupSizes} the user informs the sample size of each subsequent test. Therefore, only positive integers are accepted in \code{GroupSizes}. 

The input \code{Statistic} specifies the scale selected by the user to inform \code{cvs.lower} and \code{cvs.upper}among the classic methods:
MaxSPRT (Kulldorf et al., 2011), Pocock (Pocock, 1977), OBrien-Fleming (O'Brien and Fleming, 1979), or Wang-Tsiatis (Jennison and Turnbull, 2000). 
For \code{Statistic="Wang-Tsiatis"}, the user has to choose a number in the (0, 0.5] interval for \code{Delta}. 

For \code{RR} the user must specify the target relative risks for calculation of the statistical performance measures to be delivered in the output.
It can be a vector of positive number or a single number. 

}
\value{
\item{CV}{Signaling threshold in the scale of the selected \code{Statistic} for the user-defined alpha spending.}
\item{CV.cases}{Signaling threshold in the scale of the binomial data for the user-defined alpha spending.}
\item{ActualSpend}{The actual Type I error probability with the calculated threshold.}
\item{Performance}{A matrix with the following three performance measures for each target RR: statistical power, expected time to signal and expected sample size.}
}
\author{Ivair Ramos Silva, Martin Kulldorff.
}
\references{
Jennison C, Turnbull B. (2000). Group Sequential Methods with Applications to Clinical Trials, London: Chapman and Hall/CRC.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

O'Brien PC, Fleming TR. (1979). A multiple testing procedure for clinical trials. Biometrics. 35:549--556.

Pocock SJ. (1977).Group sequential methods in the design and analysis of clinical trials. Biometrika. 64:191--199.

Silva IR. (2018). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance with Binomial Data. Statistics in Medicine, 15;37(1), 107-118.

Silva IR, Kulldorff M. (2015). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva IR, Maro J, Kulldorff M. (2019). Exact Sequential Analysis Using R Sequential. Working Paper.
}


\section{Acknowledgements}{
Development of the Performance.AlphaSpend.Binomial function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0,2.0 to 3.1).\cr
-	Federal University of Ouro Preto (UFOP), through contract under internal UFOP's resolution CEPE 4600 (v2.0 to 3.1).\cr
\cr
}

\keyword{Signaling threshold given alpha spending with binomial data}

\section{See also}{
\code{\link[Sequential]{Performance.Threshold.Binomial}}: for calculating performance and alpha spending for user-specified signaling threshold with binomial data.\cr
\code{\link[Sequential]{CV.Binomial}}: for calculating Wald-type signaling thresholds for continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{Analyze.Binomial}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion with binomial data.
}


\examples{

## Performance and threshold for group binomial sequential analysis
#  with a maximum sample size of 50 events for upper-tailed testing
#  that is, H0:RR<=1, with irregular group sizes 15, 15, 10, and 10.
#  This is done for alpha spending of the power-type.
#  The statistical performance is evaluated for RR= 2:

# res<- Performance.AlphaSpend.Binomial(N=50,alpha=0.05,AlphaSpend=1,z=c(1,1.5,2,1.3),
# p="n",GroupSizes=c(15,15,10,10),Tailed="upper",rho=0.5, RR=2,Statistic="MaxSPRT")

}


