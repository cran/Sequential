\name{Analyze.Multinomial}
\alias{Analyze.Multinomial}
\title{Function for group sequential analyses for multinomial data, without the need to know group sizes a priori.}
\description{The function \code{Analyze.Multinomial} is used for either continuous or group sequential analysis, or for a  combination of the two. \code{Analyze.Multinomial} is run at each look at the data. Before running it by the first time, it is necessary to run the \code{\link[Sequential]{AnalyzeSetUp.Multinomial}} function.}
\usage{
Analyze.Multinomial(name,test,cases,controls,N_exposures,N_controls,AlphaSpend="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and it must be the same as the name given by the \code{AnalyzeSetup.Multinomial} function. Should never be the same as another sequential analysis that is run simultaneously on the same computer.}
\item{test}{An integer indicating the number of hypothesis tests performed up to and including the current test. For example, if there were four prior looks at the data, and this is the fifth one, then "test=5". This number should be increased by one each time that the \code{Analyze.Multinomial} function is run for a new group of data, when it is part of the same sequential analysis. If not, there is an error message.}
\item{cases}{A vector with the number of events per multinomial entry from the exposed individuals.}
\item{controls}{A vector with the number of events per multinomial entry from the unexposed individuals.}
\item{N_exposures}{A vector with the number of individuals in the same risk window per multinomial entry excluding the entry for the control group. Must have the same dimension as cases.}
\item{N_controls}{The number of individuals in the control (unexposed) group.}
\item{AlphaSpend}{The alpha spending function is specified in the \code{AnalyzeSetUp.Multinomial} function. At any look at the data, it is possible to over ride that pre-specified alpha spending plan by using the AlphaSpend parameter. AlphaSpend is a number representing the maximum amount of alpha (Type I error probabiliy) to be spent up to and including the current test. Because of the discrete nature of the multinomial distribution, the actual amount of alpha spent may be less than the maximum amount specified. It must be in the range (0, alpha]. The default value is no override, which means that, if AlphaSpend= "n", then the function will use the alpha spending plan specified in the \code{AnalyzeSetUp.Multinomial} function.}
}


\details{
The function \code{Analyze.Multinomial} performs continuous or group sequential analysis for multinomial data.
It can also be used for mixed continuous-group sequential
analysis where some data arrives continuously while other data arrives in groups. There is (i) no need to pre-specify the
group sizes before the sequential analysis starts, and (ii) a variety of alpha spending functions are available.

This function runs the adaptive sequential multiple hypotheses testing for seasonal concomitant vaccines safety surveillance introduced by Silva and Maro(2025),
which was developed for rapidly detecting increased risks of adverse events from one or more combinations of simultaneous vaccine exposures.
As explained by Silva and Maro(2025), multiple seasonal vaccines may be recommended for some strata of the population defined by covariates, such as gender
and age, for example. Such vaccines can be administered simultaneously or at least within the same vaccination season. Thus, monitoring adverse events origination from
multiple vaccines may be more informative and statistically powerful than monitoring each vaccine separately.   

The method is an extension of the binomial MaxSPRT method to a multinomial exposure model. With an exact analytical alpha spending approach,
the computationally feasible limit of multiple exposures is likely limited to no more than two vaccines.
For more complex situations with several vaccines, which can multiple types of endpoints (i.e. designated adverse events),
a valid Monte Carlo decision rule is constructed. 

The alpha spending plan is designed in a way to ensure that a target statistical power, gamma, is endured for detecting a target increased risk, r,
for each multinomial entry. This is possible by means of a new concept introduced by Silva and Maro(2025), the robust alpha spending plan.
For more details about the method and the construction of the critical values, see Silva and Maro(2025). 

Before running \code{Analyze.Multinomial}, it is necessary to specify a planned default alpha spending function,
which is done using the AlphaSpendType parameter in the
\code{AnalyzeSetUp.Multinomial} function. As detailed by Silva and Maro(2025), there are two robust alpha spending plan choices, AlphaSpendType=1 or AlphaSpendType=2.
See \code{\link[Sequential]{AnalyzeSetUp.Multinomial}} for more details. 

The function \code{Analyze.Multinomial} is meant to perform the multinomial sequential analysis with a certain level
of autonomy. After running a test, the code offers
a synthesis about the general parameter settings, the main conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information
from previous tests. Tables with the main analyses results is automatically printed in the R console. Each column of the table contains a historical characteristic,
including the information for the current test. Each line of the table corresponds to a specific test organized by calendar time. 
Observe that, because the multinomial distribution is discrete, the target alpha spending will rarely be reached. The actual alpha spending is then shown to facilitate
a realistic interpretation of the results. 
   
}

\value{
\item{result}{Table containing the main characteristics, conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information from previous tests.}
}
\author{ Ivair Ramos Silva, Judith Maro.
}


\references{
Silva IR, Maro J. (2025), Adaptive Sequential Multiple Hypotheses Testing for Seasonal Concomitant Vaccines Safety Surveillance. Working paper, Department of Population Medicine, Harvard Medical School and Harvard Pilgrim Health Care Institute. 
}


\section{Acknowledgements}{
Development of the \code{Analyze.Multinomial} function was funded by:
- Department of Population Medicine, Harvard Medical School and Harvard Pilgrim Health Care Institute.
- Conselho Nacional de Desenvolvimento Cientifico e Tecnologico (CNPq), Brazil, and Fundacao de Amparo a Pesquisa do Estado de Minas Gerais (FAPEMIG), Brazil. 
}
\keyword{Multinomial sequential analysis}

\section{See also}{
\code{\link[Sequential]{AnalyzeSetUp.Binomial}}: for setting up sequential analysis with the \code{Analyze.Binomial} function, before the first look at the data.\cr 
\code{\link[Sequential]{SampleSize.Binomial}}: for calculating the needed sample size to achieve the desired statistical power for continuous sequential analysis with binomial data.\cr 
}


\examples{

### Example. Six chunks of data.

### Firstly, it is necessary to set up the input parameters.
##  Here we use the Wald type alpha spending.
##  Note: cut off the "#" symbol before running the two lines below.

#   AnalyzeSetUp.Multinomial(name="testee",N=200,alpha=0.05,k=7,R0=1,rho=1,
#   m=10000, #title="tituto da tabela",
#   ExposuresNames=c("A","B","C","AB","AC","AD","ABC"),
#   address="C:/Users")


### Now we apply sequential tests to each of six chunks of data.
# -------------------------------------------------------------------------
  
# Test 1
#res1<- Analyze.Multinomial(name="testee",test=1,cases=c(8,6,1,4,1,0,0),controls= 3,
#N_exposures=c(1000,1500,300,300,320,280,240),N_controls=500,AlphaSpend="n")
#res1

# Test 2
#res2<- Analyze.Multinomial(name="testee",test=2,cases=c(13,1,4,6,4,2,1),controls= 5,
#N_exposures=c(900,1200,800,280,250,250,275),N_controls=1500,AlphaSpend="n")
#res2

# Test 3
#res3<- Analyze.Multinomial(name="testee",test=3,cases=c(8,5,3,3,3,1,7),controls= 3,
#N_exposures=c(800,1000,1000,230,280,310,300),N_controls=2000,AlphaSpend="n")
#res3

# Test 4
#res4<- Analyze.Multinomial(name="testee",test=4,cases=c(6,6,5,3,6,3,8),controls= 8,
# N_exposures=c(700,900,1000,200,300,350,310),N_controls=3000,AlphaSpend="n")
#res4

# Test 5
#res5<- Analyze.Multinomial(name="testee",test=5,cases=c(4,2,3,4,6,2,5),controls= 7,
#N_exposures=c(800, 850, 900, 200, 250, 300, 200),N_controls=1000,AlphaSpend="n")
#res5

# Test 6
#res6<- Analyze.Multinomial(name="testee",test=6,cases=c(6,3,3,6,7,3,6),controls= 8,
#N_exposures=c(815, 790, 880, 210, 240, 310, 208),N_controls=950,AlphaSpend="n")
#res6
}


