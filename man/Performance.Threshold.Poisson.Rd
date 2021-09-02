\name{Performance.Threshold.Poisson}
\alias{Performance.Threshold.Poisson}
\title{Performance and alpha spending for user-defined signaling threshold for sequential analysis with Poisson data.}
\description{The function \code{Performance.Threshold.Poisson} calculates power, expected time to signal and alpha spending associated to any user-specified signaling threshold, flat or non-flat, for continuous or group sequential analysis with Poisson data.
The user can select the scale for the signaling threshold among MaxSPRT, Pocock, OBrien-Fleming, or Wang-Tsiatis test statistics, all for a pre-specified upper limit on the sample size.}
\usage{
Performance.Threshold.Poisson(SampleSize,CV.lower="n",CV.upper="n",
CV.events.upper="n",M=1,D=0,GroupSizes="n",Tailed="upper",
Statistic=c("MaxSPRT", "Pocock", "OBrien-Fleming", "Wang-Tsiatis"),
Delta="n",RR)
      }
\arguments{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) expressed in terms of the expected number of events under the null hypothesis. There is no default value.}
\item{CV.lower}{Signaling threshold for evidence of "RR<1". It is given in the scale of the selected test statistic infomed in "Statistic". There is no default value.}
\item{CV.upper}{Signaling threshold for evidence of "RR>1". It is given in the scale of the selected test statistic infomed in "Statistic". There is no default value.}
\item{CV.events.upper}{Signaling threshold for evidence of "RR>1". It is given in the scale of the events. The is no default value unless when one uses "CV.upper" instead.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. It must be a positive integer.
A good rule of thumb is to set M=4 (Kulldorff and Silva, 2015). The default value is M=1, which means that even a single event
can reject the null hypothesis if it occurs sufficiently early.}
\item{D}{The expected number of events under the null hypothesis at the first look at the data.
This is used when there is an initial large chunk of data arriving, followed by continuous sequential analysis. The default value is D=0, which is
also the best choice. This means that there is no delay in the start of the sequential analyses. If D is very large, the maximum 
sample size will be set equal to D if a non-sequential analysis provides the desired power.}
\item{GroupSizes}{Vector with the test-specific expected number of events under the null hypothesis between two looks at the data. Important: Must sums up "SampleSize". There is no default value.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
\item{Statistic}{The test statistic scale used for "CV.lower" and "CV.upper". There is no default.}
\item{Delta}{Parameter needed for calculation of Wang-Tsiatis test statistic if this is the option selected in "Statistic". Must be a number in the (0, 0.5] interval. There is no default value.}
\item{RR}{Vector of relative risks for performance calculation. There is no default value.}
}
\details{
For continuous and group sequential analysis with Poisson data, the alpha spending for
user-specified thresholds are calculated with \code{Performance.Threshold.Poisson}. 


The inputs \code{CV.lower}, \code{CV.upper}, if two-tailed testing (code{Tailed="two"}),
both lower and upper signaling thresholds must be informed. If the user
desires a constant threshold (critical value), then a single number can be informed for each of these inputs.
For time-variable (non-constant) thresholds, the length of \code{CV.lower} and \code{CV.upper} must coincide with the length of \code{GroupSizes}.

   
For details about the algorithm used to calculate the critical value, see the paper by Silva (2018).

With \code{GroupSizes} the user informs the sample size of each subsequent test. Therefore, only positive numbers are accepted in \code{GroupSizes}. 

The input \code{Statistic} specifies the scale selected by the user to inform \code{CV.lower} and code{CV.upper}among the classic methods:
MaxSPRT (Kulldorf et al., 2011), Pocock (Pocock, 1977), OBrien-Fleming (O'Brien and Fleming, 1979), or Wang-Tsiatis (Jennison and Turnbull, 2000). 
For \code{Statistic="Wang-Tsiatis"}, the user has to choose a number in the (0, 0.5] interval for \code{Delta}. 

For \code{RR} the user must specify the target relative risks for calculation of the statistical performance measures to be delivered in the output.
It can be a vector of positive number or a single number. 

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

Pocock SJ. (1977). Group sequential methods in the design and analysis of clinical trials. Biometrika. 64:191--199.

Silva IR. (2018). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance With Poisson Data. Methodol Comput Appl Probab, 20(2), 739--750.

Silva IR, Kulldorff M. (2015). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva IR, Maro J, Kulldorff M. (2021). Exact sequential test for clinical trials and post-market drug and vaccine safety surveillance with Poisson and binary data. Statistics in Medicine, DOI: 10.1002/sim.9094.
}


\section{Acknowledgements}{
Development of the Performance.Threshold.Poisson function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0,2.0 to 3.1).\cr
-	Federal University of Ouro Preto (UFOP), through contract under internal UFOP's resolution CEPE 4600 (v2.0 to 3.1).\cr
\cr
}

\keyword{Performance and alpha spending given threshold with Poisson data}

\section{See also}{
\code{\link[Sequential]{Performance.AlphaSpend.Poisson}}: for calculating performance and signaling threshold for user-specified alpha spending with Poisson data.\cr
\code{\link[Sequential]{CV.Poisson}}: for calculating Wald-type signaling thresholds for continuous and group sequential analysis with Poisson data.\cr
\code{\link[Sequential]{Analyze.Poisson}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion for Poisson data.
}


\examples{

## Example 1
#  Performance and alpha spending of group Poisson sequential analysis
#  with a maximum sample size of 90 expected events for two-tailed
#  testing, i.e. H0:RR=1, with irregular group sizes and different
#  and lower and upper thresholds with irregular
#  sample sizes 25, 20, 20, and 25. 
#  The statistical performance is evaluated for four different
#  target RR= 1, 1.2, 2, 3:

# res<- Performance.Threshold.Poisson(SampleSize=90,CV.lower=c(2.5,2.6,2.7,2.8),
# CV.upper=c(3,3.1,3.2,3.3),GroupSizes=c(25,20,20,25),Tailed="two",
# Statistic="MaxSPRT",Delta="n",RR=c(1,1.2,2,3))

## Example 2
#  Suppose that the Analyze.Poisson function was used for an actual analysis.
#  For evaluating the cumulative power after a certain number of subsequent tests,
#  one can enter with the critical values delivered by Analyze.Poisson in the
#  Performance.Threshold.Poisson.
#  For example, suppose that the following thresholds in the scale of the events
#  were printed by Analyze.Poisson for the first three tests:
#     cv.events<- c(2,3,5)
#  which were obtained for the following specific sample sizes:
#     mus<- c(0.05,0.5,1.2)
#  Calculating the cumulative power, the expected time to signal, and 
#  the expected sample size for RR=2:

# res<-Performance.Threshold.Poisson(SampleSize=sum(mus),CV.events.upper=cv.events,
# GroupSizes=mus, Statistic="MaxSPRT",RR=2)


# This returns a power about 30%.

}


