\name{Analyze.wBinomial}
\alias{Analyze.wBinomial}
\title{Function for group sequential analyses of multiple weighted binomial endpoints, without the need to know group sizes a priori.}
\description{The function \code{Analyze.wBinomial} is used for either continuous or group sequential analysis, or for a  combination of the two. Unlike \code{CV.Binomial} and \code{CV.G.Binomial}, it is not necessary to pre-specify the group sizes before the sequential analysis starts. More important, this function is designed specifically when multiple outcomes with weights are analyzed. This is done using user defined alpha spending functions. For single binomial outcomes, please use \code{Analyze.Binomial}. \code{Analyze.wBinomial} is run at each look at the data. Before running it by the first time, it is necessary to run the \code{\link[Sequential]{AnalyzeSetUp.wBinomial}} function.}
\usage{
Analyze.wBinomial(name,test,z,w,ExposureA,ExposureB,AlphaSpend="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and it must be the same as the name given by the \code{AnalyzeSetup.wBinomial} function. Should never be the same as another sequential analysis that is run simultaneously on the same computer.}
\item{test}{An integer indicating the number of hypothesis tests performed up to and including the current test. For example, if there were four prior looks at the data, and this is the fifth one, then "test=5". This number should be increased by one each time that the \code{Analyze.wBinomial} function is run for a new group of data, when it is part of the same sequential analysis.  If not, there is an error message.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case. For example, if there are 3 controls matched to each case, "z=3". In a self-control analysis, z is the ratio of the length of the control interval to the length of the risk interval. For example, if the risk interval is 2 days long and the control interval is 7 days long, "z=7/2". In terms of p, the binomial probability under the null hypothesis, "p=1/(1+z)", or equivalently, "z=1/p-1". The parameter z must be a positive number. The default value is z=1 (p=0.5). If the ratio is the same for all observations, then z can be any positive number.}
\item{w}{A vector containing the weights associated to each outcome.}
\item{ExposureA}{A number or a vector of the same length as w containing the number of cases from an exposure A per outcome. The length of ExposureA equals to the number of outcomes.}
\item{ExposureB}{A number or a vector of the same length as w containing the number of cases from an exposure B per outcome. The length of ExposureB equals to the number of outcomes.}
\item{AlphaSpend}{The alpha spending function is specified in the \code{AnalyzeSetUp.wBinomial} function. At any look at the data, it is possible to over ride that pre-specified alpha spending plan by using the AlphaSpend parameter. AlphaSpend is a number representing the maximum amount of alpha (Type I error probabiliy) to be spent up to and including the current test. Because of the discrete nature of the binomial distribution, the actual amount of alpha spent may be less than the maximum amount specified. It must be in the range (0,alpha]. The default value is no override, which means that, if AlphaSpend= "n", then the function will use the alpha spending plan specified in the \code{AnalyzeSetUp.wBinomial} function.}
}


\details{
The function \code{Analyze.wBinomial} performs exact sequential testing for
multiple weighted binomial endpoints, that is, the analysis reflects the drugs combined benefit and safety profile.
It works with a variety of alpha spending functions for continuous, group or mixed group-continuous sequential analysis.
The binomial probabilities, given by \eqn{p=1/(1+z)}. 

The test statistic is based on the weighted sum of binomial endpoints introduced by Silva et al (2018). 


Unlike \code{CV.Binomial} and \code{CV.G.Binomial}, there is (i) no need to pre-specify the
group sizes before the sequential analysis starts, (ii) a variety of alpha spending functions are available,
and (iii) it is designed for multiple weighted binomial endpoints.

In sequential analysis, data is formed by cumulative information, collected in separated chunks or groups, which are observed at different moments in time. \code{Analyze.wBinomial}
is run each time a new group of data arrives at which time a new sequential test is conducted. When running \code{Analyze.wBinomial}, only the data from the new group should be
included when calling the function. The prior data has been stored, and it will be automatically retrieved by \code{Analyze.wBinomial}, with no need to reenter that data. Before
running \code{Analyze.wBinomial} for the first time, it is necessary to set up the sequential analysis using the \code{AnalyzeSetUp.wBionimial} function, which is run once,
and just once, to define the sequential analysis parameters. For information about this, see the description of the \code{\link[Sequential]{AnalyzeSetUp.wBinomial}} function. 

The function \code{Analyze.wBinomial} calculates critical values to determine if the null hypothesis should be rejected or not at each analysis.
The null hypothesis is that the relative risk of each outcome is equal to 1.

Critical values are given in the scale of the ratio \eqn{S_A/S_B}, where \eqn{S_A=w_1 ExposeA_1+w_2 ExposeA_2+...+w_k ExposeA_k},
and \eqn{S_B=w_1 ExposeB_1 +w_2 ExposeB_2+...+w_kExposeB_k}, and \eqn{k} is the length of w.  

Critical values for each test are elicited for a pre-specified overall statistical significance level (alpha), and for an upper limit on the sample size (N).
The exact analytical solution is obtained through numerical calculations. Based on the data and the critical value, the function determines if the null hypothesis
should be rejected or not, and if subsequent tests should be conducted. After each test, the function also provides information about the amount of alpha that
has been spent, the cumulative number of cases and controls, and the maximum likelihood estimate of the relative risk. 

For binomial and Bernoulli data, there are a number of 0/1 observations that can either be an ExposureA or an ExposureB. Under the null hypothesis,
the probability of being an ExposureA is p, and the probability of being an ExposureB is 1-p. If data comes from a self-control analysis,
the observation is an ExposureA if the event occurred in the risk interval, and it is an ExposureB
if the event occurred in the control interval. Under the null hypothesis, we then have that \eqn{p=1/(1+z)}, where z is the ratio of the length of
the control interval to the length of the risk interval. 

Before running \code{Analyze.wBinomial}, it is necessary to specify a planned default alpha spending function, which is done using the AlphaSpendType parameter in the
\code{AnalyzeSetUp.wBinomial} function. The default alpha spending is of the polynomial power-type
parameterized with rho, which, according to Silva (2018), 'rho=0.5' is indicated when expected time to signal is the design criterion,
hence the default in \code{AnalyzeSetUp.wBinomial}. See the \code{\link[Sequential]{AnalyzeSetUp.wBinomial}} for more details. 

In most cases, this pre-specified alpha spending function is used throughout the analysis, but if needed, it is possible to override it
at any or each of the sequential tests. This is done using the AlphaSpend parameter, which specifies the maximum amount of alpha to
spend up to and including the current test. In this way, it is possible to use any alpha spending function, and not only those available
in \code{AnalyzeSetUp.wBinomial}. It is also possible to use a flexible adaptive alpha spending plan that is not set in stone
before the sequential analysis starts. The only requirement is that for a particular test with a new group of data,
AlphaSpend must be decided before knowing the number of ExposureA and ExposureB in that group.
To ensure a statistically valid sequential analysis, AlphaSpend can only depend on the number of events (ExposureA + ExposureB)
at prior tests and the total number of events in the current test. This is important.

The function \code{Analyze.WBinomial} is meant to perform the binomial sequential analysis with a certain level of autonomy.
After running a test, the code offers a synthesis about the general parameter settings, the main conclusions concerning the
acceptance or rejection of the null hypothesis, and the historical information from previous tests.
A table with the main analyses results is automatically printed in the R console. Each column of the table contains a historical characteristic,
including the information for the current test. Each line of the table corresponds to a specific test organized by calendar time. The table is titled with the title
input defined through the function \code{AnalyzeSetUp.wBinomial}, and its columns are organized and labeled in the following way: "Test",
"#Events", "Relative Risk", "Test Statistic", "Critical Value", "Alpha", "Reject H0". Here follows a short description of each column:

-	"Test" shows the order of the analysis, i.e., the arrival order of each chunk of data. 

-	"#Events" present the total of events per outcome.

-	"Relative Risk" is the estimated relative risk per outcome. 

-	"Test Statistic" the ratio of weighted sum of binomial endpoints between ExposureA and ExposureB.

-	"Critical Value" is the signaling threshold for each test.

-	"Alpha" shows the target and actual alpha spending up to the i-th test. 

-	"Reject H0" is a logical variable that is "Yes" when the null hypothesis is rejected, and the label "No" when H0 is not to be rejected


Observe that, because the binomial distribution is discrete, the target alpha spending will rarely be reached. The actual alpha spending is then shown to facilitate
a realistic interpretation of the results. 
   
The function \code{Analyze.wBinomial} was designed to instruct the user with minimal information about bugs from the code,
or about non-applicable parameter input usages. Some entries are not applicable for the parameter inputs. For example,
the input "z" must be a positive number, and then if the user sets "z= -1", the code will report an error with the message "z must be a positive number".
Thus, messages will appear when mistakes and inconsistencies are detected, and instructions about how to proceed to solve such problems will automatically appear.
}

\value{
\item{result}{A table containing the main characteristics, conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information from previous tests.}
}
\author{ Ivair Ramos Silva, Martin Kulldorff.}


\references{
Silva IR, Gagne J, Najafzadeh M, Kulldorff M. (2020). Exact Sequential Analysis for Multiple Weighted Binomial Endpoints. Statistics in Medicine, 39(3), 340--351.

Silva IR. (2018). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance with Binomial Data. Statistics in Medicine, 15;37(1), 107--118.
}


\section{Acknowledgements}{
Development of the \code{Analyze.wBinomial} function was funded by:
- National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999: user-defined alpha spending functions, power-type alpha spending function,
increased computational speed, end of schedule analysis using left-over alpha, enhanced error handling and messages, improved documentation.
}
\keyword{Binomial sequential analysis}

\section{See also}{
\code{\link[Sequential]{AnalyzeSetUp.wBinomial}}: for setting up sequential analysis with the \code{Analyze.wBinomial} function, before the first look at the data.\cr 
}


\examples{

### Example. Four chunks of data.

### Firstly, it is necessary to set up the input parameters.
##  Note: cut off the "#" symbol before running the two lines below.

#   AnalyzeSetUp.wBinomial(name="Rofe_Naisds",N=1000,alpha=0.05,M=1,
#   rho=0.5,title="rofecoxib (Vioxx) vs. NSAID comparison",
#   address="C:/Users/Ivair/Documents",
#   Tailed=2)

### Now we apply sequential tests to each of two chunks of data.
# -------------------------------------------------------------------------

## This example is based on two outcomes, myocardinal
#  infarction (w1=2.2), and major bleeding (w2=0.04), obtained
#  from a study comparing risk of myocardial infarction and
#  gastrointestinal bleeding. See details in Silva et al (2019).
  
## Test 1
## Note: cut off the "#" symbol before running the line below.
#  Analyze.wBinomial(name="Rofe_Naisds",test=1,z=1,w=c(2.2,0.04),
#  ExposureA=c(11,12),ExposureB=c(13,10),AlphaSpend="n" )

## Test 2
## Note: cut off the "#" symbol before running the line below.
#  Analyze.wBinomial(name="Rofe_Naisds",test=2,z=c(1,1),
#  w=c(2.2,0.04),ExposureA=c(19,12),ExposureB=c(16,11),AlphaSpend="n")
}



