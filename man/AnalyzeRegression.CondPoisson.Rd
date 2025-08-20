\name{AnalyzeRegression.CondPoisson}
\alias{AnalyzeRegression.CondPoisson}
\title{Function for MaxSPRT regression analyses with conditional Poisson data, without the need to know group sizes a priori.}
\description{The function \code{AnalyzeRegression.CondPoisson} is used for either continuous or group sequential analysis, or for a  combination of the two. Unlike \code{CV.CondPoisson}, it is not necessary to pre-specify the group sizes before the sequential analysis starts. \code{AnalyzeRegression.CondPoisson} is run at each look at the data. Before running it by the first time, it is necessary to run the \code{\link[Sequential]{AnalyzeSetUpRegression.CondPoisson}} function.}
\usage{
AnalyzeRegression.CondPoisson(name,test,events,PersonTimeRatio="n",
covariates,AlphaSpend="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and it must be the same as the name given by the \code{AnalyzeSetupRegression.CondPoisson} function. Should never be the same as another sequential analysis that is run simultaneously on the same computer.}
\item{test}{An integer indicating the number of hypothesis tests performed up to and including the current test. For example, if there were four prior looks at the data, and this is the fifth one, then "test=5". This number should be increased by one each time that the \code{AnalyzeRegression.CondPoisson} function is run for a new group of data, when it is part of the same sequential analysis.  If not, there is an error message.}
\item{events}{The test specific vector with the number of observed events per covariates combination.}
\item{PersonTimeRatio}{Vector with the ratio between the punctual, instead of cumulative from previous
tests, person-time observed in the current test, by the total person-time observed in the historical period, per covariates combination.}
\item{covariates}{Matrix with the covariates for the regression model. The i-th line of covariates has the information related to the i-th entry of the vector events. Each column of covariates corresponds to a different explanatory covariate.}
\item{AlphaSpend}{The alpha spending function is specified in the \code{AnalyzeSetUpRegression.CondPoisson} function. At any look at the data, it is possible to over ride that pre-specified alpha spending plan by using the AlphaSpend parameter. AlphaSpend is a number representing the maximum amount of alpha (Type I error probabiliy) to be spent up to and including the current test. Because of the discrete nature of the Poisson distribution, the actual amount of alpha spent may be less than the maximum amount specified. It must be in the range (0,alpha]. The default value is no override, which means that, with the default AlphaSpend= "n", the function will use the alpha spending plan specified in the \code{AnalyzeSetUpRegression.CondPoisson} function.}
}



\details{
The function \code{AnalyzeRegression.CondPoisson} performs continuous or group MaxSPRT regression analysis for conditional Poisson data based on the method
proposed by Silva et al.(2025).

It can also be used for mixed continuous-group sequential analysis where some data arrives continuously while other data arrives in groups. 

Unlike \code{CV.CondPoisson}, there is (i) no need to pre-specify the
group sizes before the sequential analysis starts, (ii) a variety of alpha spending functions are available, and (iii) it is possible to include an offset term where,
under the null hypothesis, different observations have different Poisson person-time ratios between the surveillance and the historical information.

In sequential analysis, data is formed by cumulative information, collected in separated chunks or groups, which are observed at different moments in time. \code{AnalyzeRegression.CondPoisson}
is run each time a new group of data arrives at which time a new sequential test is conducted. When running \code{AnalyzeRegression.CondPoisson}, only the data from the new group should be
included when calling the function. The prior data has been stored, and it will be automatically retrieved by \code{AnalyzeRegression.CondPoisson}, with no need to reenter that data. Before
running \code{AnalyzeRegression.CondPoisson} for the first time, it is necessary to set up the sequential analysis using the \code{AnalyzeSetUpRegression.CondPoisson} function, which is run once,
and just once, to define the sequential analysis parameters. For information about this, see the description of the \code{\link[Sequential]{AnalyzeSetUpRegression.CondPoisson}} function. 



The function \code{AnalyzeRegression.CondPoisson} calculates critical values to determine if the null hypothesis should be rejected or not at each analysis.
Critical values are given by the value of the AlphaSpending function, and the test statistic is the Monte Carlo p-value calculated as proposed by Silva et al.(2025).
This way, the null hypothesis H0:RR<= R0, where R0 is the testing margin given by the analyst in the \code{AnalyzeSetUpRegression.CondPoisson},
and RR is the true unknow relative risk to test. 


After each test, the function also provides information about the amount of alpha that has been spent,
the cumulative number of cases and controls, regression coefficient estimates, confidence intervals,
and the maximum likelihood estimate of the relative risk per stratum of individuals observed during the
sequential analysis.



Before running \code{AnalyzeRegression.CondPoisson}, it is necessary to specify a planned default alpha spending function,
which is done using the AlphaSpendType parameter in the
\code{AnalyzeSetUpRegression.CondPoisson} function. The default alpha spending plan is the polynomial power-type alpha spending plan
parameterized with "rho=1". Different alpha spending plans can be obtained by selecting different values for rho (Silva et al, 2019).

In most cases, this pre-specified alpha spending function is used throughout the analysis, but if needed,
it is possible to override it at any or each of the sequential tests.
This is done using the AlphaSpend parameter in \code{AnalyzeRegression.CondPoisson}, which specifies the
maximum amount of alpha to spend up to and including the current test.
In this way, it is possible to use any
alpha spending function, and not only the power-type available in \code{AnalyzeSetUpRegression.CondPoisson}.
It means that the AlphaSpend parameter can be used to promote a flexible adaptive alpha spending plan that is not set in stone
before the sequential analysis starts. The only requirement is that for a particular test with a new group of data,
AlphaSpend must be decided before knowing the number of events in that group.
To ensure a statistically valid sequential analysis, AlphaSpend can only depend on the values of PersonTimeRatio. This is important.

The function \code{AnalyzeRegression.CondPoisson} is meant to perform the Poisson sequential analysis with a certain level of autonomy.
After running a test, the code offers a synthesis about
the general parameter settings, the main conclusions concerning the acceptance or rejection of the null hypothesis,
and the historical information from previous tests.


Observe that, because the Poisson distribution is discrete, the target alpha spending will rarely be reached. The actual alpha spending is then shown to facilitate
a realistic interpretation of the results. 
   
The function \code{AnalyzeRegression.CondPoisson} was designed to instruct the user with minimal information about bugs from the code,
or about non-applicable parameter input usages. Some entries are not
applicable for the parameter inputs.  

}

\value{
\item{result}{Four data.frames (Decision_table, Relative_risk_estimates, Coefficients, Confidence_Intervals) with the main information concerning the tuning parameterization for the planned surveillance and the historical information about the performed tests.}
}
\author{ Ivair Ramos Silva.
}


\references{
Jennison C, Turnbull B. (2000). Group Sequential Methods with Applications to Clinical Trials. London: Chapman and Hall/CRC.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. Sequential Analysis, 30, 58--78.

Kulldorff M, Silva IR. (2015). Continuous post-market sequential safety surveillance with minimum
events to signal. arxiv:1503.01978 [stat.ap].

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71(3), 851--858. 

Silva IR, Kulldorff M, Yih W. Katherine. (2020), Optimal alpha spending for sequential analysis with binomial data. Journal of the Royal Statistical Society Series B, 82(4) p. 1141--1164.

Silva IR, Lopes LM, Dias P, Yih WK. (2019), Alpha Spending for Historical Versus Surveillance Poisson Data With CMaxSPRT. Statistics in Medicine, 38(12), 2126--2138.

Silva IR, Montalban, J., Oliveira, F. (2025), Maximized Sequential Probability Ratio Test Regression. Working paper - Sentinel (TIDE) project, Department of Population Medicine, Harvard Medical School and Harvard Pilgrim Health Care Institute.
}


\section{Acknowledgements}{
Development of the \code{AnalyzeRegression.CondPoisson} function was funded by:
-	Food and Drug Administration, Center for Drug Evaluation and Research, through the Mini-Sentinel Project (base version, documentation);\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (user defined alpha spending functions, improved documentation);\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil, process number 302882/2022-7. (v4.3.1).\cr
-	Support Foundation to Minas Gerais State Research-Fapemig, Brazil, grant numbers PQ-00787-21 and RED-00133-21. (v3.1 to v4.3). \cr 
}
\keyword{Poisson sequential analysis}

\section{See also}{
\code{\link[Sequential]{AnalyzeSetUpRegression.CondPoisson}}: for setting up sequential analysis with the \code{AnalyzeRegression.CondPoisson} function, before the first look at the data.\cr  
}


\examples{

### Example. Three chunks of data.

### Firstly, it is necessary to set up the input parameters.
##  Note: cut off the "#" symbol before running the lines below.
#     AnalyzeSetUpRegression.CondPoisson(name="VaccineA",N=100,cc=20,alpha=0.05,
#     R0=1, rho=1,mref=999,title="Monitoring_vaccineA",
#     address="C:/Users/Ivair/Documents")


### Now we apply sequential tests to each of three chunks of data.
# -------------------------------------------------------------------------
  
## Test 1 - Situation of three covariates combinations, i.e.,
## three lines in the covariates matrix, and two covariates.
## Regarding the PersonTimeRatio parameter, in this first test
## we illustrate the situation with the same information 0.05 per
## covariates combination. For example, if the total historical person-time
## is V= 5000, and if the test-specific person-time (per covariates combination)
## is P=250, then PersonTimeRatio= c(250,250,250)/5000,
## which gives PersonTimeRatio=c(0.05,0.05,0.05).  
## Note: cut off the "#" symbol before running the lines below.
#  AnalyzeRegression.CondPoisson(name= "VaccineA",test=1,events= c(3,2,4),
#  PersonTimeRatio= c(0.05,0.05,0.05),
#  covariates=matrix(c(0.1,1,0.2,0,0.3,1),3,2)  )

## Test 2 - Situation where there is only one event per covariates combination:
## Note: cut off the "#" symbol before running the line below.
#  AnalyzeRegression.CondPoisson(name= "VaccineA",test=2,events= c(1,1,1),
#  PersonTimeRatio= c(0.05,0.05,0.05),
#  covariates=matrix(c(0.3,0.9,0.25,0.1,0.35,15),3,2) )

## Test 3 - Situation of elevated number of events and only two groups:
## Note: cut off the "#" symbol before running the line below.
#  AnalyzeRegression.CondPoisson(name= "VaccineA",test=3,events= c(10,13),
#  PersonTimeRatio= c(0.25,0.3),
#  covariates=matrix(c(0.5,1,0.35,0.2),2,2) )
 

}


