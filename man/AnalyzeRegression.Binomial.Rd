\name{AnalyzeRegression.Binomial}
\alias{AnalyzeRegression.Binomial}
\title{Function for MaxSPRT regression analyses with binary/binomial data, without the need to know group sizes a priori.}
\description{The function \code{AnalyzeRegression.Binomial} is used for either continuous or group sequential analysis, or for a  combination of the two. Unlike \code{CV.Binomial}, it is not necessary to pre-specify the group sizes before the sequential analysis starts. Moreover, under the null hypothesis, the binomial probability, p, can be different for different observations. In a matched case-control setting, this means that the matching ratios can be different for different matched sets. \code{AnalyzeRegression.Binomial} is run at each look at the data. Before running it by the first time, it is necessary to run the \code{\link[Sequential]{AnalyzeSetUpRegression.Binomial}} function.}
\usage{
AnalyzeRegression.Binomial(name,test,z="n",p="n",cases,controls,covariates,AlphaSpend="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and it must be the same as the name given by the \code{AnalyzeSetupRegression.Binomial} function. Should never be the same as another sequential analysis that is run simultaneously on the same computer.}
\item{test}{An integer indicating the number of hypothesis tests performed up to and including the current test. For example, if there were four prior looks at the data, and this is the fifth one, then "test=5". This number should be increased by one each time that the \code{AnalyzeRegression.Binomial} function is run for a new group of data, when it is part of the same sequential analysis.  If not, there is an error message.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case. For example, if there are 3 controls matched to each case, "z=3". In a self-control analysis, z is the ratio of the length of the control interval to the length of the risk interval. For example, if the risk interval is 2 days long and the control interval is 7 days long, "z=7/2". In terms of p, the binomial probability under the null hypothesis, "p=1/(1+z)", or equivalently, "z=1/p-1". The parameter z must be a positive number. The default value is z=1 (p=0.5). If the ratio is the same for all observations, then z can be any positive number. If the ratio is different for different observations, then z is a vector of positive numbers.}
\item{p}{The probability of having a case under the null hypothesis. There is no default value.}
\item{cases}{A number or a vector of the same length as z containing the number of cases per stratum of individuals defined by the matrix covariates.}
\item{controls}{A number or a vector of the same length as z containing the number of controls per stratum of individuals defined by the matrix covariates.}
\item{covariates}{Matrix with the covariates for the regression model. The i-th line of covariates has the information related to the i-th entry of cases and controls. Each column of covariates corresponds to a different explanatory covariate.}
\item{AlphaSpend}{The alpha spending function is specified in the \code{AnalyzeSetUpRegression.Binomial} function. At any look at the data, it is possible to over ride that pre-specified alpha spending plan by using the AlphaSpend parameter. AlphaSpend is a number representing the maximum amount of alpha (Type I error probabiliy) to be spent up to and including the current test. Because of the discrete nature of the binomial distribution, the actual amount of alpha spent may be less than the maximum amount specified. It must be in the range (0,alpha]. The default value is no override, which means that, if AlphaSpend= "n", then the function will use the alpha spending plan specified in the \code{AnalyzeSetUpRegression.Binomial} function.}
}


\details{
The function \code{AnalyzeRegression.Binomial} performs continuous or group MaxSPRT regression analysis for Bernoulli or binomial data based on the method
proposed by Silva et al.(2025).

It can also be used for mixed continuous-group sequential analysis where some data arrives continuously while other data arrives in groups. 

Unlike \code{CV.Binomial}, there is (i) no need to pre-specify the
group sizes before the sequential analysis starts, (ii) a variety of alpha spending functions are available, and (iii)  it is possible to include an offset term z where,
under the null hypothesis, different observations have different binomial probabilities p.

In sequential analysis, data is formed by cumulative information, collected in separated chunks or groups, which are observed at different moments in time. \code{AnalyzeRegression.Binomial}
is run each time a new group of data arrives at which time a new sequential test is conducted. When running \code{AnalyzeRegression.Binomial}, only the data from the new group should be
included when calling the function. The prior data has been stored, and it will be automatically retrieved by \code{AnalyzeRegression.Binomial}, with no need to reenter that data. Before
running \code{AnalyzeRegression.Binomial} for the first time, it is necessary to set up the sequential analysis using the \code{AnalyzeSetUpRegression.Binomial} function, which is run once,
and just once, to define the sequential analysis parameters. For information about this, see the description of the \code{\link[Sequential]{AnalyzeSetUpRegression.Binomial}} function. 

The function \code{AnalyzeRegression.Binomial} calculates critical values to determine if the null hypothesis should be rejected or not at each analysis.
Critical values are given by the value of the AlphaSpending function, and the test statistic is the Monte Carlo p-value calculated as proposed by Silva et al.(2025).
This way, the null hypothesis H0:RR<= R0, where R0 is the testing margin given by the analyst in the \code{AnalyzeSetUpRegression.Binomial},
and RR is the true unknow relative risk to test. 


After each test, the function also provides information about the amount of alpha that has been spent,
the cumulative number of cases and controls, regression coefficient estimates, confidence intervals,
and the maximum likelihood estimate of the relative risk per stratum of individuals observed during the
sequential analysis.



For binomial and Bernoulli data, there are a number of 0/1 observations that can either be a case or a control. Under the null hypothesis,
the probability of being a case is p, and the probability of being a control is 1-p. If data comes from a self-control analysis,
the observation is a case if the event occurred in the risk interval, and it is a control
if the event occurred in the control interval. Under the null hypothesis, we then have that \eqn{p=1/(1+z)}, where z is the ratio of the
length of the control interval to the length of the
risk interval. This ratio, and hence p, does not need to be the same for all observations.

If data comes from a matched set of exposed and unexposed individuals, then the observation
is a case if the event occurred among one of the exposed, and it is a control if it occurred
among one of the unexposed. Under the null hypothesis, \eqn{p=1/(1+z)}, where z is the number of
unexposed individuals divided by the number of exposed individuals in the matched set. Again,
this ratio does not have to be the same for all matched sets.  The variable z can be any
positive number. 

The ratio parameter z is a vector, representing multiple z values. For each value of z, it is necessary to specify the number of cases and the number
of controls. This means that for a chunk of data, the vector of zs has to be of the same length as the vector of cases and the vector of controls.
Therefore, the first entry of the vector z is the matching ratio associated to the first entries of cases and of controls. The second entry of z is the matching ratio with respect to the second entries of
cases and of controls, and so on. For example, consider that each of five observations came from four different matching ratios.
In this situation, the vectors cases, controls and z are all of length four. For example, suppose "z=c(2,1,0.5,3)", "cases=c(1,1,0,0)" and "controls=c(0,0,1,2)".
The matching ratio for the first observation, which turned out as a case,
is equal to 2. For the second observation, also a case, the matching is equal to 1. With a matching ration of 0.5, the third observation turned out to be a control. The two last
observations both had a matching ratio of 3, and both of them were controls. To complete this example, the "covariates" input has four lines (one line per "cases" entry). For this example, suppose that two explanatory continuous variables
are available, covariates=matrix(c(0.1,1.01, 0.2,0, 0.3,1, 0.15,0.9),4,2).

If all observations in the same data group has the same ratio, the vectors are of size one, that is,
they are simple numbers. For example, if there were ten observations that all had a ratio of 2, with seven cases and three controls, we have "z=2",
"cases=7", and "controls=3". In this case, the the "covariates" matrix has only one line. For example, covariates=matrix(c(0.1,1.01),1,2), which represents
the situation with two explanatory variables (covariates has two columns).

Alternatively, instead of z the user can specify p directly.
Note that only one of these inputs, z or p, has to be specified, but if both are entered the code will only work if z and p are such that p=1/(1+z).
Otherwise, an error message will appear to remind that such condition must be complied.




Before running \code{AnalyzeRegression.Binomial}, it is necessary to specify a planned default alpha spending function,
which is done using the AlphaSpendType parameter in the
\code{AnalyzeSetUpRegression.Binomial} function. The default alpha spending plan is the polynomial power-type alpha spending plan
parameterized with "rho=1". Different alpha spending plans can be obtained by selecting different values for rho (Silva, 2018).

In most cases, this pre-specified alpha spending function is used throughout the analysis, but if needed,
it is possible to override it at any or each of the sequential tests.
This is done using the AlphaSpend parameter in \code{AnalyzeRegression.Binomial}, which specifies the
maximum amount of alpha to spend up to and including the current test.
In this way, it is possible to use any
alpha spending function, and not only the power-type available in \code{AnalyzeSetUpRegression.Binomial}.
It means that the AlphaSpend parameter can be used to promote a flexible adaptive alpha spending plan that is not set in stone
before the sequential analysis starts. The only requirement is that for a particular test with a new group of data,
AlphaSpend must be decided before knowing the number of cases
and controls in that group. To ensure a statistically valid sequential analysis, AlphaSpend can only depend on the number of events (cases + controls)
at prior tests and the total number of events in the current test. This is important.

The function \code{AnalyzeRegression.Binomial} is meant to perform the binomial sequential analysis with a certain level of autonomy.
After running a test, the code offers a synthesis about
the general parameter settings, the main conclusions concerning the acceptance or rejection of the null hypothesis,
and the historical information from previous tests.


Observe that, because the binomial distribution is discrete, the target alpha spending will rarely be reached. The actual alpha spending is then shown to facilitate
a realistic interpretation of the results. 
   
The function \code{AnalyzeRegression.Binomial} was designed to instruct the user with minimal information about bugs from the code,
or about non-applicable parameter input usages. Some entries are not
applicable for the parameter inputs. For example,
the input "z" must be a positive number, and then if the user sets "z= c(-1,0,2)",
the code will report an error with the message "the entries of
the vector "z" must be positive numbers". Thus, messages will appear when mistakes and inconsistencies are detected,
and instructions about how to proceed to solve such problems will automatically appear. 

}

\value{
\item{result}{Four data.frames (Decision_table, Relative_risk_estimates, Coefficients, Confidence_Intervals) with the main information concerning the tuning parameterization for the planned surveillance and the historical information about the performed tests.}
}
\author{ Ivair Ramos Silva.
}


\references{
Fireman B, et al. (2013). Exact sequential analysis for binomial data with time varying probabilities. Manuscript in preparation.

Jennison C, Turnbull B. (2000). Group Sequential Methods with Applications to Clinical Trials. London: Chapman and Hall/CRC.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. Sequential Analysis, 30, 58--78.

Kulldorff M, Silva IR. (2015). Continuous post-market sequential safety surveillance with minimum
events to signal. arxiv:1503.01978 [stat.ap].

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71(3), 851--858. 

Silva IR, Kulldorff M, Yih W. Katherine. (2020), Optimal alpha spending for sequential analysis with binomial data. Journal of the Royal Statistical Society Series B, 82(4) p. 1141--1164.

Silva IR. (2018). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance with Binomial Data. Statistics in Medicine, 15;37(1), 107-118.

Silva IR, Zhuang, Y. (2022), Bounded-width confidence interval following optimal sequential analysis of adverse events with binary data, Statistical Methods in Medical Research, 31(12), 2323--2337.

Silva IR, Montalban, J., Oliveira, F. (2025), Maximized Sequential Probability Ratio Test Regression. Working paper - Sentinel (TIDE) project, Department of Population Medicine, Harvard Medical School and Harvard Pilgrim Health Care Institute.
}


\section{Acknowledgements}{
Development of the \code{AnalyzeRegression.Binomial} function was funded by:
-	Food and Drug Administration, Center for Drug Evaluation and Research, through the Mini-Sentinel Project (base version, documentation);\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (user defined alpha spending functions, improved documentation);\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil, process number 302882/2022-7. (v4.3.1).\cr
-	Support Foundation to Minas Gerais State Research-Fapemig, Brazil, grant numbers PQ-00787-21 and RED-00133-21. (v3.1 to v4.3). \cr 
}
\keyword{Binomial sequential analysis}

\section{See also}{
\code{\link[Sequential]{AnalyzeSetUpRegression.Binomial}}: for setting up sequential analysis with the \code{AnalyzeRegression.Binomial} function, before the first look at the data.\cr  
}


\examples{

### Example. Three chunks of data.

### Firstly, it is necessary to set up the input parameters.
##  Note: cut off the "#" symbol before running the lines below.
#     AnalyzeSetUpRegression.Binomial(name="VaccineA",N=200,alpha=0.05,
#     R0=1, rho=1,mref=999,title="Monitoring_vaccineA",
#     address="C:/Users/Ivair/Documents")

### Now we apply sequential tests to each of three chunks of data.
# -------------------------------------------------------------------------
  
## Test 1 - Situation where each stratum came from a different
## matching ratio.
## This first test uses the default power-type (rho=1) alpha spending (AlphaSpend="n").
## Note: cut off the "#" symbol before running the lines below.
#  AnalyzeRegression.Binomial(name= "VaccineA",test=1,z=c(1.1,1.3,1.2,1),
#  cases= c(1,0,0,0), controls= c(0,1,1,1),  
#  covariates=matrix(c(0.1,1.01, 0.2,0, 0.3,1, 0.15,0.9),4,2)  )

## Test 2 - Situation where some of the strata came from the same matching
## ratio.
## Observe that here we use an arbitrary alpha spending of 0.02.
## Note: cut off the "#" symbol before running the line below.
#  AnalyzeRegression.Binomial(name= "VaccineA",test=2,z=c(1,1,1.5),cases= c(12,1,4),
#  controls= c(0,10,5), covariates=matrix(c(0.3,0.9, 0.5,0.4, 0.55,1.1),3,2),
#  AlphaSpend=0.02)

## Test 3 - Situation of elevated number of events, but now the
## arbitrary alpha spending is of 0.04, and p is entered instead of z.
## Note: cut off the "#" symbol before running the line below.
#  AnalyzeRegression.Binomial(name= "VaccineA",test=3,p=c(0.6,0.4),cases= c(15,12),
#  controls= c(13,16), covariates=matrix(c(0.1,1.1, 0.6,0.35),2,2), 
#  AlphaSpend=0.04)
 

}


