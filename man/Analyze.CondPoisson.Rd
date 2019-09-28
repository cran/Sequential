\name{Analyze.CondPoisson}
\alias{Analyze.CondPoisson}
\title{Function to conduct group sequential analyses for conditional Poisson data without the need to know group sizes a priori.}
\description{The function \code{Analyze.CondPoisson} is used for either continuous, group, or mixed continuous-group sequential analysis for Poisson data conditioned on observed historical data. Unlike \code{CV.CondPoisson}, it is not necessary that data arrives in a near-continuous fashion. It is possible to use either a Wald type rejection boundary, which is flat with respect to the likelihood ratio, or a user defined alpha spending function. \code{Analyze.CondPoisson} is run at each look at the data. Before running it by the first time, it is necessary to run the \code{\link[Sequential]{AnalyzeSetUp.CondPoisson}} function.}
\usage{
Analyze.CondPoisson(name,test,events,PersonTimeRatio,AlphaSpend="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and it must be the same as the name given by the \code{AnalyzeSetup.CondPoisson} function. Should never be the same as another sequential analysis that is run simultaneously on the same computer.}
\item{test}{An integer indicating the number of hypothesis tests performed up to and including the current test. For example, if there were four prior looks at the data, and this is the fifth one, then "test=5". This number should be increased by one each time that the \code{Analyze.CondPoisson} function is run for a new group of data, when it is part of the same sequential analysis.  If not, there is an error message.}
\item{events}{The number of observed events during the surveillance period.}
\item{PersonTimeRatio}{The observed ratio between the punctual, instead of cumulative from previous tests, person-time observed in the current test, by the total person-time observed in the historical period.}
\item{AlphaSpend}{The alpha spending function is specified in the \code{AnalyzeSetUp.CondPoisson} function. At any look at the data, it is possible to over ride that pre-specified alpha spending plan by using the AlphaSpend parameter. AlphaSpend is a number representing the maximum amount of alpha (Type I error probabiliy) to be spent up to and including the current test. Because of the discrete nature of the Poisson distribution, the actual amount of alpha spent may be less than the maximum amount specified. It must be in the range (0,alpha]. The default value is no override, which means that, if AlphaSpend= "n", then the function will use the alpha spending plan specified in the \code{AnalyzeSetUp.CondPoisson} function.}
}


\details{

The function \code{Analyze.CondPoisson} performs continuous or group sequential analysis for Poisson data conditioned on observed historical data, (Li and Kulldorff, 2010).
It can also be used for mixed continuous-group sequential analysis where some data arrives continuously while other data arrives in groups.
Unlike \code{CV.CondPoisson}, there is a variety of alpha spending functions are available.

In sequential analysis, data is formed by cumulative information, collected in separated chunks or groups, which are observed at different moments in time.
\code{Analyze.CondPoisson} is run each time a new group of data arrives at which time a new sequential test is conducted. When running \code{Analyze.CondPoisson},
only the data from the new group should be included when calling the function. The prior data has been stored, and it will be automatically retrieved by
\code{Analyze.CondPoisson}, with no need to reenter that data. Before
running \code{Analyze.CondPoisson} for the first time, it is necessary to set up the sequential analysis using the \code{AnalyzeSetUp.CondPoisson} function, which
is run once, and just once, to define the sequential analysis parameters. For information about this, see the description of the
\code{\link[Sequential]{AnalyzeSetUp.CondPoisson}} function. 

The function \code{Analyze.CondPoisson} calculates critical values to determine if the null hypothesis should be rejected or not at each analysis.
Critical values are given in the
scale of the number of events. This is done for a pre-specified overall statistical significance level (alpha), and for an upper limit on the sample size, which is
given by "T" or "K". Go to the documentation of \code{\link[Sequential]{AnalyzeSetUp.CondPoisson}} for more details about the choice between "T" or "K".

The exact analytical solution is obtained through numerical calculations. Based on the data and the critical value, the function determines if the null hypothesis should
be rejected or not, and if subsequent tests should be conducted. After each test, the function also provides information about the amount of alpha that has been spent,
the cumulative number of events, and the observed log-likelihood ratio statistic. 


Before running \code{Analyze.CondPoisson}, it is necessary to specify a planned default alpha spending function, which is done using the AlphaSpendType parameter in the
\code{AnalyzeSetUp.CondPoisson} function. The default alpha spending plan can be either, (i) the polynomial power-type alpha spending plan, which is parameterized with rho,
or (ii) the alpha spending associated to the Wald-type rejection boundary, which is flat with respect to the likelihood ratio. See the \code{\link[Sequential]{AnalyzeSetUp.CondPoisson}}
for more details. 

In most cases, this pre-specified alpha spending function is used throughout the analysis, but if needed, it is possible to override it at any or each of the sequential tests.
This is done using the AlphaSpend parameter, which specifies the maximum amount of alpha to spend up to and including the current test. In this way, it is possible to use any
alpha spending function, and not only those available in \code{AnalyzeSetUp.CondPoisson}. It is also possible to use a flexible adaptive alpha spending plan that is not set in stone
before the sequential analysis starts. The only requirement is that for a particular test with a new group of data, AlphaSpend must be decided before knowing the PersonTimeRatio
in that group. Hence, in order to ensure a statistically valid sequential analysis, AlphaSpend can only depend on the cumulative events.
This is important.

The function \code{Analyze.CondPoisson} is meant to perform the conditional Poisson sequential analysis with a certain level of autonomy. After running a test, the code offers a synthesis about
the general parameter settings, the main conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information from previous tests.
A table with the main analyses results is automatically printed in the R console. Each column of the table contains a historical characteristic,
including the information for the current test. Each line of the table corresponds to a specific test organized by calendar time. The table is titled with the title
input defined through the function \code{AnalyzeSetUp.CondPoisson}, and its columns are organized and labeled in the following way:
"Test","Person-timeR","events","Cumulative Person-timeR","Cumulative events","LLR", "target","actual","CV","Reject H0". Here follows a short description of each column:

-	"Test" shows the order of the analysis, i.e., the arrival order of each chunk of data. 

-	"Person-timeR" shows the observed ratio between the punctual person-time observed in the current test by the total person-time observed in the historical period.

-	"Events" presents the observed number of events from the Poisson counting entered at each test.

-	"Cumulative Person-timeR" shows the observed person-time ratio up to the current test.

-	"Cumulative events" presents the observed number of events from the Poisson counting up to the current test.

-	"LLR" is the observed log-likelihood ratio test statistic.

-	"target" is the target alpha spending for the i-th test.

-	"actual" is the actual alpha spent up to the i-th test. 

-	"CV" is the critical value in the scale of the log-likelihood ratio test statistic.

-	"Reject H0" is a logical variable that receives the label "Yes" when the null hypothesis is rejected, and the label "No" when H0 is not to be rejected


Observe that, depending on the choices of the input parameters M and alpha through the \code{\link[Sequential]{AnalyzeSetUp.CondPoisson}} function,
the actual alpha spending can differ from the target one. The actual alpha spending is then shown in order to favor a realistic interpretation of the results. 
   
The function \code{Analyze.CondPoisson} was designed to instruct the user with minimal information about bugs from the code, or about non-applicable input parameters
usage. Some entries are not applicable. For example, the input "Person-timeR" must be a positive number, and then if the user sets "Person-timeR= -1",
then the code will report an error with the message "the entry of "Person-timeR" must be a number greater than zero". Thus, messages will appear when mistakes
and inconsistencies are detected. Instructions about how to proceed to solve such problems will automatically appear too. 
}

\value{
\item{result}{A table containing the main characteristics, conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information from previous tests.}
}
\author{ Ivair Ramos Silva, Martin Kulldorff.
}


\references{
Fireman B, et al. (2013). Exact sequential analysis for Poisson data with time varying probabilities. Manuscript in preparation.

Jennison C, Turnbull B. (2000). Group Sequential Methods with Applications to Clinical Trials. London: Chapman and Hall/CRC.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. Sequential Analysis, 30, 58--78.

Kulldorff M, Silva IR. (2015). Continuous Post-market Sequential Safety Surveillance with Minimum Events to Signal. REVSTAT Statistical Journal, 15(3): 373--394.

Li L, Kulldorff M. (2010). A conditional maximized sequential probability ratio test for Pharmacovigilance. Statistics in Medicine, 29 (2), 284--295. 

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71(3), 851--858.

Silva IR, Li L, Kulldorff M. (2019). Exact conditional maximized sequential probability ratio test adjusted for covariates. Sequential Analysis, 38(1), 115--133.

Silva IR., Lopes LM., Dias P., Yih WK. (2019). Alpha Spending for Historical Versus Surveillance Poisson Data With CMaxSPRT. Statistics in Medicine, DOI: 10.1002/sim.8097, 1--13. 

}


\section{Acknowledgements}{
Development of the \code{Analyze.CondPoisson} function was funded by:
- National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999.
- Foundation for Research Support of Minas Gerais State (FAPEMIG), MG, Brazil, through the grant Demanda Universal.
}
\keyword{Poisson sequential analysis}

\section{See also}{
\code{\link[Sequential]{AnalyzeSetUp.CondPoisson}}: for setting up sequential analysis with the \code{Analyze.CondPoisson} function, before the first look at the data.\cr
\code{\link[Sequential]{AnalyzeSetUp.Poisson}}: for setting up sequential analysis with the \code{Analyze.Poisson} function, before the first look at the data.\cr 
\code{\link[Sequential]{SampleSize.Poisson}}: for calculating the needed sample size to achieve the desired statistical power for continuous sequential analysis with Poisson data.\cr
}


\examples{

### Example. Four chunks of data.

### Firstly, it is necessary to set up the input parameters.
##  Here we use the Wald type alpha spending.
##  Note: cut off the "#" symbol before running the two lines below, and,
##  Important: choose an actual "address" to save your set up information. 
#   AnalyzeSetUp.CondPoisson(name="TestA", SampleSizeType="Events", K=100,
#   cc=20,alpha=0.05, M=1,AlphaSpendType="power-type",rho=0.5,title="n",
#   address="C:/Users/Ivair/Documents")

### Now we apply a test for each one of three chunks of data.
# -------------------------------------------------------------------------
  
## Test 1 - Situation where the expected number of events is equal to 5.
## The observed hypothetical information is "PersonTimeRatio=0.5".
## Note: cut off the "#" symbol before running the line below.

#  Analyze.CondPoisson(name="TestA",test=1,events=5,PersonTimeRatio=0.5)


## Test 2 - Situation where the expected number of events is equal to 6.
## The observed hypothetical information is "PersonTimeRatio=0.3".

#Analyze.CondPoisson(name="TestA",test=2,events=6,PersonTimeRatio=0.3)


## Test 3 - Situation where the expected number of events is equal to 10.
## The observed hypothetical information is "PersonTimeRatio=0.1".

#Analyze.CondPoisson(name="TestA",test=3,events=10,PersonTimeRatio=0.1)
 
}


