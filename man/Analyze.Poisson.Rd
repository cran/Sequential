\name{Analyze.Poisson}
\alias{Analyze.Poisson}
\title{Function to conduct group sequential analyses for Poisson data without the need to know group sizes a priori.}
\description{The function \code{Analyze.Poisson} is used for either continuous or group sequential analysis, or for a  combination of the two. Unlike \code{CV.Poisson} and \code{CV.G.Poisson}, it is not necessary to pre-specify the group sizes before the sequential analysis starts. Moreover, under the null hypothesis, the expected number of events, mu0, can be different for different observations.  It is possible to use either a Wald type rejection boundary, which is flat with respect to the likelihood ratio, or a user defined alpha spending function. \code{Analyze.Poisson} is run at each look at the data. Before running it by the first time, it is necessary to run the \code{\link[Sequential]{AnalyzeSetUp.Poisson}} function.}
\usage{
Analyze.Poisson(name,test,mu0,events,AlphaSpend="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and it must be the same as the name given by the \code{AnalyzeSetup.Poisson} function. Should never be the same as another sequential analysis that is run simultaneously on the same computer.}
\item{test}{An integer indicating the number of hypothesis tests performed up to and including the current test. For example, if there were four prior looks at the data, and this is the fifth one, then "test=5". This number should be increased by one each time that the \code{Analyze.Poisson} function is run for a new group of data, when it is part of the same sequential analysis.  If not, there is an error message.}
\item{mu0}{The expected number of events under the null hypothesis. The parameter mu0 must be a positive number. There is no default value.}
\item{events}{The number of observed events.}
\item{AlphaSpend}{The alpha spending function is specified in the \code{AnalyzeSetUp.Poisson} function. At any look at the data, it is possible to over ride that pre-specified alpha spending plan by using the AlphaSpend parameter. AlphaSpend is a number representing the maximum amount of alpha (Type I error probabiliy) to be spent up to and including the current test. Because of the discrete nature of the Poisson distribution, the actual amount of alpha spent may be less than the maximum amount specified. It must be in the range (0,alpha]. The default value is no override, which means that, if AlphaSpend= "n", then the function will use the alpha spending plan specified in the \code{AnalyzeSetUp.Poisson} function.}
}


\details{
The function \code{Analyze.Poisson} performs continuous or group sequential analysis for Poisson data. It can also be used for mixed continuous-group sequential
analysis where some data arrives continuously while other data arrives in groups. Unlike \code{CV.Poisson} and \code{CV.G.Poisson}, there is (i) no need to pre-specify the
group sizes before the sequential analysis starts, (ii) a variety of alpha spending functions are available, and (iii)  it is possible to include an offset term where,
under the null hypothesis, different observations have different Poisson rates mu0.

In sequential analysis, data is formed by cumulative information, collected in separated chunks or groups, which are observed at different moments in time. \code{Analyze.Poisson}
is run each time a new group of data arrives at which time a new sequential test is conducted. When running \code{Analyze.Poisson}, only the data from the new group should be
included when calling the function. The prior data has been stored, and it will be automatically retrieved by \code{Analyze.Poisson}, with no need to reenter that data. Before
running \code{Analyze.Poisson}  for the first time, it is necessary to set up the sequential analysis using the \code{AnalyzeSetUp.Poisson} function, which is run once,
and just once, to define the sequential analysis parameters. For information about this, see the description of the \code{\link[Sequential]{AnalyzeSetUp.Poisson}} function. 

The function \code{Analyze.Poisson} calculates critical values to determine if the null hypothesis should be rejected or not at each analysis. Critical values are given in the
scale of the number of events. This is done for a pre-specified overall statistical significance level (alpha), and for an upper limit on the sample size (N).
The exact analytical solution is obtained through numerical calculations. Based on the data and the critical value, the function determines if the null hypothesis should
be rejected or not, and if subsequent tests should be conducted. After each test, the function also provides information about the amount of alpha that has been spent,
the cumulative number of events, and the estimated relative risk. 


Before running \code{Analyze.Poisson}, it is necessary to specify a planned default alpha spending function, which is done using the AlphaSpendType parameter in the
\code{AnalyzeSetUp.Poisson} function. The default alpha spending plan can be either, (i) the polynomial power-type alpha spending plan, which is parameterized with rho,
and the default is rho=0.5 as suggested by Silva (2018), or (ii) the alpha spending associated to the Wald-type rejection boundary, which is flat with respect to the likelihood ratio. See the \code{\link[Sequential]{AnalyzeSetUp.Poisson}}
for more details. 

In most cases, this pre-specified alpha spending function is used throughout the analysis, but if needed, it is possible to override it at any or each of the sequential tests.
This is done using the AlphaSpend parameter, which specifies the maximum amount of alpha to spend up to and including the current test. In this way, it is possible to use any
alpha spending function, and not only those available in \code{AnalyzeSetUp.Poisson}. It is also possible to use a flexible adaptive alpha spending plan that is not set in stone
before the sequential analysis starts. The only requirement is that for a particular test with a new group of data, AlphaSpend must be decided before knowing the number of events
in that group. To ensure a statistically valid sequential analysis, AlphaSpend can only depend on cumulative mu0 values at prior tests and of the mu0 value in the current test.
This is important.

The function \code{Analyze.Poisson} is meant to perform the Poisson sequential analysis with a certain level of autonomy. After running a test, the code offers a synthesis about
the general parameter settings, the main conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information from previous tests.
A table with the main analyses results is automatically printed in the R console. Each column of the table contains a historical characteristic,
including the information for the current test. Each line of the table corresponds to a specific test organized by calendar time. The table is titled with the title
input defined through the function \code{AnalyzeSetUp.Poisson}, and its columns are organized and labeled in the following way: "Test",
"mu0", "Events", "Cumulative mu0", "Cumulative Events", "RR", "LLR", "target", "actual", "CV", "Reject H0". Here follows a short description of each column:

-	"Test" shows the order of the analysis, i.e., the arrival order of each chunk of data. 

-     "mu0" is the expected number of events under the null hypothesis for the chunck of data to be analyzed at each test. 

-	"Events" presents the observed number of events from the Poisson counting entered at each test.

-     "Cumulative mu0" expected number of events under the null hypothesis up to the i-th test.

-     "Cumulative Events" observed number of events up to the i-th test.
 
-	"RR" is the estimated relative risk for test i. 

-	"LLR" is the observed log-likelihood ratio test statistic.

-	"target" is the target alpha spending for the i-th test.

-	"actual" is the actual alpha spent up to the i-th test. 

-	"CV" is the critical value in the scale of the number of events, showing how many events are needed to reject the null hypothesis at this test.

-	"Reject H0" is a logical variable that is "Yes" when the null hypothesis is rejected, and the label "No" when H0 is not to be rejected


Observe that, because the Poisson distribution is discrete, the target alpha spending will rarely be reached. The actual alpha spending is then shown to facilitate
a realistic interpretation of the results. 
   
The function \code{Analyze.Poisson} was designed to instruct the user with minimal information about bugs from the code, or about non-applicable parameter input usages. Some entries are not
applicable for the parameter inputs. For example, the input "mu0" must be a positive number, and then if the user sets "mu0= -1", the code will report an error with the message "the entry of "mu0" must be a number greater than zero". Thus, messages will appear when mistakes and inconsistencies are detected, and instructions about how to proceed to solve such problems will automatically appear. 
}

\value{
\item{result}{A table containing the main characteristics, conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information from previous tests.}
}
\author{ Ivair Ramos Silva,  Ned Lewis, Martin Kulldorff.
}


\references{
Fireman B, et al. (2013). Exact sequential analysis for Poisson data with time varying probabilities. Manuscript in preparation.

Jennison C, Turnbull B. (2000). Group Sequential Methods with Applications to Clinical Trials. London: Chapman and Hall/CRC.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. Sequential Analysis, 30, 58--78.

Kulldorff M, Silva IR. (2015). Continuous Post-market Sequential Safety Surveillance with Minimum Events to Signal. REVSTAT Statistical Journal, 15(3): 373--394.

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva IR. (2018). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance With Poisson Data. Methodol Comput Appl Probab, 20(2), 739-750.

}


\section{Acknowledgements}{
Development of the \code{Analyze.Poisson} function was funded by:
- National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999.
}
\keyword{Poisson sequential analysis}

\section{See also}{
\code{\link[Sequential]{AnalyzeSetUp.Poisson}}: for setting up sequential analysis with the \code{Analyze.Poisson} function, before the first look at the data.\cr 
\code{\link[Sequential]{Performance.G.Poisson}}: for calculating the statistical power, expected time to signal and expected sample size for group sequential analysis with Poisson data.\cr
\code{\link[Sequential]{SampleSize.Poisson}}: for calculating the needed sample size to achieve the desired statistical power for continuous sequential analysis with Poisson data.\cr 
\code{\link[Sequential]{CV.G.Poisson}}: for calculating critical values for group sequential analysis with Poisson data.
}


\examples{

### Example. Four chunks of data.

### Firstly, it is necessary to set up the input parameters.
##  Here we use the Wald type alpha spending.
##  Note: cut off the "#" symbol before running the two lines below, and,
##  very important, choose an actual "address" to save your set up information. 
#   AnalyzeSetUp.Poisson(name="VaccineA", SampleSize=100, alpha=0.05,
#   M=1,AlphaSpendType="power-type",rho=0.5,title="n",
#   address="C:/Users/Ivair/Documents")

### Now we can sequentially apply a test for each one of three chunks of data.
# -------------------------------------------------------------------------
  
## Test 1 - Situation where the expected number of events under H0 is equal to 2.
## This first test uses the default Wald type alpha spending (AlphaSpend="n").
## Note: cut off the "#" symbol before running the line below.
#  Analyze.Poisson(test=1,mu0=2,events=1,AlphaSpend="n")

## Test 2 - Situation where the expected number of events under H0 is equal to 0.8.
## Observe that here we use an arbitrary alpha spending of 0.02.
## Note: cut off the "#" symbol before running the line below.
#  Analyze.Poisson(test=2,mu0=0.8,events=2, AlphaSpend=0.02)

## Test 3 - Situation of elevated number of events, but now the
## arbitrary alpha spending is of 0.04.
## Note: cut off the "#" symbol before running the line below.
#  Analyze.Poisson(test=3,mu0=9,events=10, AlphaSpend=0.04)
 
}


