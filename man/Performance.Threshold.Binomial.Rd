\name{Performance.Threshold.Binomial}
\alias{Performance.Threshold.Binomial}
\title{Statistical Performance and Alpha Spending For User-defined Signaling Threshold With Binomial Data.}
\description{The function \code{Performance.Threshold.Binomial} calculates power, expected time to signal, expected sample size and alpha spending associated to any user-specified signaling threshold, flat or non-flat, for continuous or group sequential analysis with binomial data.
The user can select the scale for the signaling threshold among MaxSPRT, Pocock, OBrien-Fleming, or Wang-Tsiatis test statistics. Alternatively, the threshold can be informed also in the scale of the binomial data.}
\usage{
Performance.Threshold.Binomial(N,CV.lower="n",CV.upper="n",z="n",p="n",
GroupSizes="n",Tailed="upper",Statistic=c("MaxSPRT", "Pocock",
"OBrien-Fleming", "Wang-Tsiatis","Cases"),Delta="n",RR)
      }
\arguments{
\item{N}{The upper limit on the sample size (length of surveillance) expressed in terms of the total number of events (cases plus controls). There is no default value.}
\item{CV.lower}{Signaling threshold for evidence of "RR<1". It is given in the scale of the selected test statistic infomed in "Statistic". There is no default value.}
\item{CV.upper}{Signaling threshold for evidence of "RR>1". It is given in the scale of the selected test statistic infomed in "Statistic". There is no default value.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case under the null hypothesis. There is no default value.}
\item{p}{The probability of having a case under the null hypothesis. If just a single number is given, then it will be used as a constant probability for all groups. Otherwise, the dimension of p must coincide with the dimension of GroupSizes. There is no default value.}
\item{GroupSizes}{Vector with the total number of events (cases+controls) between two looks at the data with regular and irregular group sizes. Important: Must sums up N. For continuos sequential analysis, specify GroupSizes=1. There is no default value.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
\item{Statistic}{The test statistic scale used for "CV.lower" and "CV.upper". There is no default.}
\item{Delta}{Parameter needed when "Statistic=Wang-Tsiatis" is selected. Must be a number in the (0, 0.5] interval. There is no default value.}
\item{RR}{Vector of relative risks for performance calculation. There is no default value.}
}
\details{
For continuous and group sequential analysis with binomial data, alpha spending for
user-specified thresholds are calculated with \code{Performance.Threshold.Binomial}.

\code{N} must be a positive integer defining the maximum length of surveillance. To avoid very large computation times,
we suggest not using values greater than 1000.

For two-tailed testing (\code{Tailed="two"}), both lower and upper signaling thresholds must be informed through
\code{CV.lower} and \code{CV.upper}. If the user
desires a constant threshold (critical value) in the scale of a test statistic, then a single number can be informed.
For time-variable (non-constant) thresholds, the length of \code{CV.upper} and \code{CV.lower} must coincide with the length of \code{GroupSizes}.



\code{z} is a vector of positive numbers representing the matching ratios for each test (group). If a single number is given, then it will be used as a constant
matching ratio for all tests (groups). Otherwise, the dimension of \code{z} must coincide with the dimension of \code{GroupSizes}.
\code{z} represents the number of controls matched to each case. For example, if there are 3 controls matched to each case, \code{z}=3. 
In a self-control analysis, z is the ratio of the control interval to the risk interval.
For example, if the risk interval is 2 days long and the control interval is 7 days long, \code{z}=7/2.
In terms of \code{p}, the binomial probability under the null hypothesis, \code{p=1/(1+z)}, or equivalently, \code{z=1/p-1}.

Alternatively, instead of \code{z} the user can specify \code{p} directly.
Note that only one of these inputs, \code{z} or \code{p}, has to be specified, but if both are entered the code will only work if \code{z}
and \code{p} are such that \code{p=1/(1+z)}. Otherwise, an error message will appear to remind that such condition must be complied.
   

With \code{GroupSizes} the user informs the sample size of each subsequent test. Therefore, only positive integers are accepted in \code{GroupSizes}. 

The input \code{Statistic} specifies the scale selected by the user to inform \code{CV.lower} and \code{cvs.upper}among the classic methods:
MaxSPRT (Kulldorf et al., 2011), Pocock (Pocock, 1977), OBrien-Fleming (O'Brien and Fleming, 1979), or Wang-Tsiatis (Jennison and Turnbull, 2000). 
For \code{Statistic="Wang-Tsiatis"}, the user has to choose a number in the (0, 0.5] interval for \code{Delta}.

Important: for time-variable matching ratios (i.e. when \code{z} or \code{p} changes from a test to another), only the "Statistic=Cases" option works.
This is because the test statistic options are non-monotone with the number of cumulative cases under a variable \code{p} or \code{z} situation.

For \code{RR} the user must specify the target relative risks for calculation of the statistical performance measures to be delivered in the output.
It can be a vector of positive number or a single number. 

For details about the algorithm used to calculate the critical value, see the paper by Silva (2018).

}
\value{
\item{AlphaSpend}{The alpha spending associated to the user-specified threshold.}
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

Silva IR, Maro J, Kulldorff M. (2021), Exact sequential test for clinical trials and post-market drug and vaccine safety surveillance with Poisson and binary data. Statistics in Medicine, DOI: 10.1002/sim.9094.
}


\section{Acknowledgements}{
Development of the Performance.Threshold.Binomial function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0,2.0 to 3.1).\cr
-	Federal University of Ouro Preto (UFOP), through contract under internal UFOP's resolution CEPE 4600 (v2.0 to 3.1).\cr
\cr
}

\keyword{Alpha spending given threshold with binomial data}

\section{See also}{
\code{\link[Sequential]{Performance.AlphaSpend.Binomial}}: for calculating signaling threshold for user-specified alpha spending with binomial data.\cr
\code{\link[Sequential]{CV.Binomial}}: for calculating Wald-type signaling thresholds for continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{Analyze.Binomial}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion with binomial data.
}


\examples{

## Performance and Alpha spending of a four-group sequential
#  analysis with threshold informed in the scale of the
#  binomial data, i.e. Statistic="Cases". 
#  The analysis is for a maximum sample size of 50 events under 
#  upper-tailed testing, that is, H0:RR<=1, with irregular group
#  sizes of 12, 25, 35, and 45. 
#  The matching ratio also changes in time with z= 1, 1.5, 2, 1.3.
#  The statistical performance is evaluated for RR= 1.2, 1.5, 2:

# res<- Performance.Threshold.Binomial(N=50,CV.upper=c(12,25,35,45),
# z=c(1,1.5,2,1.3),GroupSizes=c(15,15,10,10),Tailed="upper",
# Statistic="Cases", RR=c(1.2,1.5,2))

}


