\name{Analyze.Binomial}
\alias{Analyze.Binomial}
\title{Function to Conduct Group Sequential Analyses for Binomial Data When the Goup Sizes are not Known a Priori.}
\description{The function \code{Analyze.Binomial} is used for either continuous or group sequential analysis, or for a  combination of the two. Unlike \code{CV.Binomial} and \code{CV.G.Binomial}, it is not necessary to pre-specify the group sizes before the sequential analysis starts. Moreover, under the null hypothesis, the binomial probability, p, can be different for different observations. In a matched case-control setting, this means that the matching ratios can be different for different matched sets.  It is possible to use either a Wald type rejection boundary, which is flat with respect to the likelihood ratio, or a user defined alpha spending function. \code{Analyze.Binomial} is run at each look at the data. Before running it by the first time, it is necessary to run the \code{\link[Sequential]{AnalyzeSetUp.Binomial}} function.}
\usage{
Analyze.Binomial(name,test,z="n",p="n",cases,controls,AlphaSpend="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and it must be the same as the name given by the \code{AnalyzeSetup.Binomial} function. Should never be the same as another sequential analysis that is run simultaneously on the same computer.}
\item{test}{An integer indicating the number of hypothesis tests performed up to and including the current test. For example, if there were four prior looks at the data, and this is the fifth one, then "test=5". This number should be increased by one each time that the \code{Analyze.Binomial} function is run for a new group of data, when it is part of the same sequential analysis.  If not, there is an error message.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case. For example, if there are 3 controls matched to each case, "z=3". In a self-control analysis, z is the ratio of the length of the control interval to the length of the risk interval. For example, if the risk interval is 2 days long and the control interval is 7 days long, "z=7/2". In terms of p, the binomial probability under the null hypothesis, "p=1/(1+z)", or equivalently, "z=1/p-1". The parameter z must be a positive number. The default value is z=1 (p=0.5). If the ratio is the same for all observations, then z can be any positive number. If the ratio is different for different observations, then z is a vector of positive numbers.}
\item{p}{The probability of having a case under the null hypothesis. There is no default value.}
\item{cases}{A number or a vector of the same length as z containing the number of cases.}
\item{controls}{A number or a vector of the same length as z containing the number of controls.}
\item{AlphaSpend}{The alpha spending function is specified in the \code{AnalyzeSetUp.Binomial} function. At any look at the data, it is possible to over ride that pre-specified alpha spending plan by using the AlphaSpend parameter. AlphaSpend is a number representing the maximum amount of alpha (Type I error probabiliy) to be spent up to and including the current test. Because of the discrete nature of the binomial distribution, the actual amount of alpha spent may be less than the maximum amount specified. It must be in the range (0,alpha]. The default value is no override, which means that, if AlphaSpend= "n", then the function will use the alpha spending plan specified in the \code{AnalyzeSetUp.Binomial} function.}
}


\details{
The function \code{Analyze.Binomial} performs continuous or group sequential analysis for Bernoulli or binomial data. It can also be used for mixed continuous-group sequential
analysis where some data arrives continuously while other data arrives in groups. Unlike \code{CV.Binomial} and \code{CV.G.Binomial}, there is (i) no need to pre-specify the
group sizes before the sequential analysis starts, (ii) a variety of alpha spending functions are available, and (iii)  it is possible to include an offset term where,
under the null hypothesis, different observations have different binomial probabilities p.

In sequential analysis, data is formed by cumulative information, collected in separated chunks or groups, which are observed at different moments in time. \code{Analyze.Binomial}
is run each time a new group of data arrives at which time a new sequential test is conducted. When running \code{Analyze.Binomial}, only the data from the new group should be
included when calling the function. The prior data has been stored, and it will be automatically retrieved by \code{Analyze.Binomial}, with no need to reenter that data. Before
running \code{Analyze.Binomial}  for the first time, it is necessary to set up the sequential analysis using the \code{AnalyzeSetUp.Bionimial} function, which is run once,
and just once, to define the sequential analysis parameters. For information about this, see the description of the \code{\link[Sequential]{AnalyzeSetUp.Binomial}} function. 

The function \code{Analyze.Binomial} calculates critical values to determine if the null hypothesis should be rejected or not at each analysis. Critical values are given in the
scale of the number of cases. This is done for a pre-specified overall statistical significance level (alpha), and for an upper limit on the sample size (N).
The exact analytical solution is obtained through numerical calculations. Based on the data and the critical value, the function determines if the null hypothesis should
be rejected or not, and if subsequent tests should be conducted. After each test, the function also provides information about the amount of alpha that has been spent,
the cumulative number of cases and controls, and the estimated relative risk. When the null hypothesis is rejected, or when the sequential analysis ends without rejecting the null,
then the function calculates a confidence interval for the relative risk. 

For binomial and Bernoulli data, there are a number of 0/1 observations that can either be a case or a control. Under the null hypothesis, the probability of being a case is p,
and the probability of being a control is 1-p. If data comes from a self-control analysis, the observation is a case if the event occurred in the risk interval, and it is a control
if the event occurred in the control interval. Under the null hypothesis, we then have that \eqn{p=1/(1+z)}, where z is the ratio of the length of the control interval to the length of the
risk interval. This ratio, and hence p, does not need to be the same for all observations. If data comes from a matched set of exposed and unexposed individuals, then the observation
is a case if the event occurred among one of the exposed, and it is a control if it occurred among one of the unexposed. Under the null hypothesis, \eqn{p=1/(1+z)}, where z is the number of
unexposed individuals divided by the number of exposed individuals in the matched set. Again, this ratio does not have to be the same for all matched sets.  The variable z can be any
positive number. 

If the ratio parameter z, and hence p, is the same for all observations in the same group of data, then z is just a positive number. On the other hand, if different observations in the
same group of data have different values for z, then z is a vector, representing multiple z values. For each value of z, it is necessary to specify the number of cases and the number
of controls. This means that for a group of data, the vector of zs has to be of the same length as the vector of cases and the vector of controls. The first entry of the vector
z is the matching ratio associated to the first entries of cases and of controls. The second entry of z is the matching ratio with respect to the second entries of cases and
of controls, and so on. For example, consider that each of five observations came from four different matching ratios. In this situation, the vectors cases, controls and z are
all of length four. For example, suppose "z=c(2,1,0.5,3)",  "cases=c(1,1,0,0)" and "controls=c(0,0,1,2)". The matching ratio for the first observation, which turned out as a case,
is equal to 2. For the second observation, also a case, the matching is equal to 1. With a matching ration of 0.5, the third observation turned out to be a control. The two last
observations both had a matching ratio of 3, and both of them were controls. If all observations in the same data group has the same ratio, the vectors are of size one, that is,
they are simple numbers. For example, if there were ten observations that all had a ratio of 2, with seven cases and three controls, we have "z=2", "cases=7", and "controls=3". 

Alternatively, instead of z the user can specify p directly.
Note that only one of these inputs, z or p, has to be specified, but if both are entered the code will only work if z and p are such that p=1/(1+z).
Otherwise, an error message will appear to remind that such condition must be complied.

Before running \code{Analyze.Binomial}, it is necessary to specify a planned default alpha spending function, which is done using the AlphaSpendType parameter in the
\code{AnalyzeSetUp.Binomial} function. The default alpha spending plan can be either, (i) the polynomial power-type alpha spending plan, which is parameterized with rho,
or (ii) the alpha spending associated to the Wald-type rejection boundary, which is flat with respect to the likelihood ratio. See the \code{\link[Sequential]{AnalyzeSetUp.Binomial}}
for more details. 

In most cases, this pre-specified alpha spending function is used throughout the analysis, but if needed, it is possible to override it at any or each of the sequential tests.
This is done using the AlphaSpend parameter, which specifies the maximum amount of alpha to spend up to and including the current test. In this way, it is possible to use any
alpha spending function, and not only those available in \code{AnalyzeSetUp.Binomial}. It is also possible to use a flexible adaptive alpha spending plan that is not set in stone
before the sequential analysis starts. The only requirement is that for a particular test with a new group of data, AlphaSpend must be decided before knowing the number of cases
and controls in that group. To ensure a statistically valid sequential analysis, AlphaSpend can only depend on the number of events (cases + controls) at prior tests and the total
number of events in the current test. This is important.

The function \code{Analyze.Binomial} is meant to perform the binomial sequential analysis with a certain level of autonomy. After running a test, the code offers a synthesis about
the general parameter settings, the main conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information from previous tests.
A table with the main analyses results is automatically printed in the R console. Each column of the table contains a historical characteristic,
including the information for the current test. Each line of the table corresponds to a specific test organized by calendar time. The table is titled with the title
input defined through the function \code{AnalyzeSetUp.Binomial}, and its columns are organized and labeled in the following way: "Test",
"Cases", "Controls", "Cumulative Cases", "Cumulative Controls", "Cumulative E[Cases]", "RR", "LLR", "target", "actual", "CV", "H0 rejected". Here follows a short description of each column:

-	"Test" shows the order of the analysis, i.e., the arrival order of each chunk of data. 

-	"Cases" and "Controls" present the total of cases and controls that entered at each test, respectively.

-	"Cumulative Cases" and "Cumulative Controls" in the i-th line have the cumulative counts of cases and controls up to the i-th test, respectively. 

-	"Cumulative E[Cases]" in line i is the expected cumulative number of cases for the i-th test under the null hypothesis.

-	"RR" is the estimated relative risk for test i. 

-	"LLR" is the observed log-likelihood ratio test statistic.

-	"target" is the target alpha spending for the i-th test.

-	"actual" is the actual alpha spent up to the i-th test. 

-	"CV" is the critical value in the scale of the number of cases, showing how many casesa re needed to reject the null hypothesis at this test.

-	"H0 rejected" is a logical variable that is "Yes" when the null hypothesis is rejected, and the label "No" when H0 is not to be rejected


Observe that, because the binomial distribution is discrete, the target alpha spending will rarely be reached. The actual alpha spending is then shown to facilitate
a realistic interpretation of the results. 
   
The function \code{Analyze.Binomial} was designed to instruct the user with minimal information about bugs from the code, or about non-applicable parameter input usages. Some entries are not
applicable for the parameter inputs. For example, the input "z" must be a positive number, and then if the user sets "z= -1", the code will report an error with the message "the entries of
the vector "z" must be positive numbers". Thus, messages will appear when mistakes and inconsistencies are detected, and instructions about how to proceed to solve such problems will automatically appear. 
}

\value{
\item{result}{A table containing the main characteristics, conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information from previous tests.}
}
\author{ Ivair Ramos Silva,  Ned Lewis, Martin Kulldorff.
}


\references{
Fireman B, et al. (2013). Exact sequential analysis for binomial data with time varying probabilities. Manuscript in preparation.

Jennison C, Turnbull B. (2000). Group Sequential Methods with Applications to Clinical Trials. London: Chapman and Hall/CRC.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. Sequential Analysis, 30, 58--78.

Kulldorff M, Silva IR. (2015). Continuous post-market sequential safety surveillance with minimum
events to signal. arxiv:1503.01978 [stat.ap].

Silva IR, Kulldorff M. (2015). Continuous versus Group Sequential Analysis for Vaccine and Drug
Safety Surveillance. Biometrics, in press. 
}


\section{Acknowledgements}{
Development of the \code{Analyze.Binomial} function was funded by:
- Food and Drug Administration, Center for Drug Evaluation and Research, through Mini-Sentinel Project: base version, documentation, unequal matching ratios;\cr
- National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999: user-defined alpha spending functions, power-type alpha spending function,
increased computational speed, confidence intervals for relative risks, end of schedule analysis using left-over alpha, enhanced error handling and messages, improved documentation.

We thank Claudia Coronel-Moreno for valuable editorial support, Bruce Fireman for general guidance, and Josh Gagne for important feedback on the unequal matching ratio feature. 
}
\keyword{Binomial sequential analysis}

\section{See also}{
\code{\link[Sequential]{AnalyzeSetUp.Binomial}}: for setting up sequential analysis with the \code{Analyze.Binomial} function, before the first look at the data.\cr 
\code{\link[Sequential]{Performance.G.Binomial}}: for calculating the statistical power, expected time to signal and expected sample size for group sequential analysis with binomial data.\cr
\code{\link[Sequential]{SampleSize.Binomial}}: for calculating the needed sample size to achieve the desired statistical power for continuous sequential analysis with binomial data.\cr 
\code{\link[Sequential]{CV.G.Binomial}}: for calculating critical values for group sequential analysis with binomial data.\cr
\code{\link[Sequential]{CV.G.Poisson}}: for calculating critical values for group sequential analysis with Poisson data.
}


\examples{

### Example. Four chunks of data.

### Firstly, it is necessary to set up the input parameters.
##  Here we use the Wald type alpha spending.
##  Note: cut off the "#" symbol before running the two lines below.
#     AnalyzeSetUp.Binomial(name="VaccineA",N=200,alpha=0.05,zp=1,M=3,AlphaSpendType="Wald",
#     title="Monitoring_vaccineA", address="C:/Users/Ivair/Documents")

### Now we can sequentially apply a test for each one of four chunks of data.
# -------------------------------------------------------------------------
  
## Test 1 - Situation where each individual event came from a different matching ratio.
## This first test uses the default Wald type alpha spending (AlphaSpend="n").
## Note: cut off the "#" symbol before running the line below.
#  Analyze.Binomial(name= "VaccineA",test=1,z=c(1.1,1.3,1.2,1),cases= c(1,0,0,0),
#  controls= c(0,1,1,1))

## Test 2 - Situation where some of the events came from the same matching ratio.
## Observe that here we use an arbitrary alpha spending of 0.02.
## Note: cut off the "#" symbol before running the line below.
#  Analyze.Binomial(name= "VaccineA",test=2,z=c(1,1.5),cases= c(12,1),controls= c(0,10),
#  AlphaSpend=0.02)

## Test 3 - Situation of elevated number of events, but now the
## arbitrary alpha spending is of 0.04, and p is entered instead of z.
## Note: cut off the "#" symbol before running the line below.
#  Analyze.Binomial(name= "VaccineA",test=3,p=c(0.4,0.5),cases= c(12,10),controls= c(10,14),
#  AlphaSpend=0.04)
 
## Test 4 - Situation where all the events came from the same matching ratio.
## Here the original target alpha spending is used.
## Note: cut off the "#" symbol before running the line below.
#  Analyze.Binomial(name= "VaccineA",test=4,z=2,cases= 20,controls= 10)
}


