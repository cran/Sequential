\name{Analyze.Multinomial}
\alias{Analyze.Multinomial}
\title{Function for group sequential analyses for multinomial data, without the need to know group sizes a priori.}
\description{The function \code{Analyze.Multinomial} is used for either continuous or group sequential analysis, or for a  combination of the two. \code{Analyze.Multinomial} is run at each look at the data. Before running it by the first time, it is necessary to run the \code{\link[Sequential]{AnalyzeSetUp.Multinomial}} function.}
\usage{
Analyze.Multinomial(name,test,cases,controls,N_exposures,N_controls,
exposure_group,strata_group_cases="n",strata_group_controls="n",
AlphaSpend="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and it must be the same as the name given by the \code{AnalyzeSetup.Multinomial} function. Should never be the same as another sequential analysis that is run simultaneously on the same computer.}
\item{test}{An integer indicating the number of hypothesis tests performed up to and including the current test. For example, if there were four prior looks at the data, and this is the fifth one, then "test=5". This number should be increased by one each time that the \code{Analyze.Multinomial} function is run for a new group of data, when it is part of the same sequential analysis. If not, there is an error message.}
\item{cases}{A vector with the number of events per combination of exposure_group and strata_group_cases. There is no default.}
\item{controls}{A vector with the number of events per per strata_group_controls. There is no default.}
\item{N_exposures}{A vector with the number of individuals in the same risk window per combination of exposure_group and strata_group_cases excluding the entry for the control group. Must have the same dimension as cases.}
\item{N_controls}{The number of individuals in the control (unexposed) group per strata_group_controls. There is no default.}
\item{exposure_group}{Labels to indicate the exposure group of each entry in "cases". Each entry in "exposure_group" must belong to the set of labels in "ExposuresNames", and "ExposuresNames" is one of the inputs of the AnalyzeSet.Multinomial function. There is no default.}
\item{strata_group_cases}{Labels to indicate the strata group, e.g., by combining the categorical covariates for age and gender, of each entry in "cases". The default is strata_group_cases="n" to indicate that there are no covariates to adjust for.}
\item{strata_group_controls}{Labels to indicate the strata group, e.g., by combining the categorical covariates for age and gender, of each entry in "controls". The default is strata_group_controls="n" to indicate that there are no covariates to adjust for.}
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

The alpha spending plan is designed in a way to ensure that a target statistical power, gamma, is ensured for detecting a target increased risk, r1,
for each multinomial entry. This is possible by means of a new concept introduced by Silva and Maro(2025), the robust alpha spending plan.
For more details about the method and the construction of the critical values, see Silva and Maro(2025). 

Before running \code{Analyze.Multinomial}, it is necessary to specify a planned default alpha spending function,
which is done using the AlphaSpendType parameter in the
\code{AnalyzeSetUp.Multinomial} function. As detailed by Silva and Maro(2025), there are two robust alpha spending plan choices, AlphaSpendType=1 or AlphaSpendType=2.
See \code{\link[Sequential]{AnalyzeSetUp.Multinomial}} for more details. 

The function \code{Analyze.Multinomial} is meant to perform the multinomial sequential analysis with a certain level
of autonomy. After running a test, the code offers
a synthesis about the general parameter settings, the main conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information
from previous tests. Tables with the main analyses results are automatically printed in the R console. The columns of the tables show historical characteristics,
including the information for the current test. Each line of the tables corresponds to a specific test organized by calendar time. 
Observe that, because the multinomial distribution is discrete, the target alpha spending will rarely be reached. The actual alpha spending is then shown to facilitate
a realistic interpretation of the results.

To adjust for covariates confounders, \code{Analyze.Multinomial} utilizes the maxLik package (Henningsen et al., 2011) to compute the MLE of the baseline covariates
effects as derived by Silva and Maro(2025). 
   
}

\value{
\item{result}{Table containing the main characteristics, conclusions concerning the acceptance or rejection of the null hypothesis, and the historical information from previous tests.}
}
\author{Ivair Ramos Silva, Judith Maro.
}


\references{
Silva IR, Maro J. (2025), Adaptive Sequential Multiple Hypotheses Testing for Seasonal Concomitant Vaccines Safety Surveillance. Working paper, Department of Population Medicine, Harvard Medical School and Harvard Pilgrim Health Care Institute. 

Henningsen A, Toomet O (2011), maxLik: A package for maximum likelihood estimation in R. Computational Statistics, 26(3), 443--458.
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

### Example. Two chunks of data.

### Firstly, it is necessary to set up the input parameters.
##  Note: cut off the "#" symbol before running the four lines below.

#   AnalyzeSetUp.Multinomial(name="VaccinesABC",N=200,alpha=0.05,R0=1,rho=1,
#   m=10000, title="Title of the output table",
#   ExposuresNames=c("A","B","C","AB","AC","BC","ABC"),
#   address="C:/Users")


### Running sequential tests for two chunks of data.
# -------------------------------------------------------------------------
  
## Test 1: in this first test, events were observed for each of
#          the exposures "A","B","C","AB","AC","BC","ABC", and
#          two binary covariates are present to adjust for gender
#          and age with categories "F-Y", "F-O", "M-Y", "M-O".             
# Note: cut off the "#" symbol before running the lines below.


#res1<- Analyze.Multinomial(name="VaccinesABC", test=1,
#cases= c(2,1,0,1,0,1,0,
#         1,2,0,1,0,1,1,
#         3,3,1,1,0,0,1,
#         2,0,1,1,1,0,0),

#controls= c(3,3,3,3),

#N_exposures= c(1000,1500,300,300,320,280,240,
#               1000,1500,300,300,320,280,240,
#               1000,1500,300,300,320,280,240,
#               1000,1500,300,300,320,280,240),

#N_controls= c(500,500,500,500),

#exposure_group= c("A","B","C","AB","AC","BC","ABC",
#                  "A","B","C","AB","AC","BC","ABC",
#                  "A","B","C","AB","AC","BC","ABC", 
#                  "A","B","C","AB","AC","BC","ABC"),

#strata_group_cases= c("F-Y","F-Y","F-Y","F-Y","F-Y","F-Y","F-Y",
#                      "F-O","F-O","F-O","F-O","F-O","F-O","F-O",
#                      "M-Y","M-Y","M-Y","M-Y","M-Y","M-Y","M-Y",
#                      "M-O","M-O","M-O","M-O","M-O","M-O","M-O"),

#strata_group_controls= c("F-Y","F-O","M-Y","M-O"),

#AlphaSpend= "n"            )



## Test 2: in this second test, events were observed for each of
#          the exposures "A","B","C","AB","AC","BC","ABC", but
#          there are no covariates to adjust for.             
# Note: cut off the "#" symbol before running the lines below.

#res2<- Analyze.Multinomial(name="VaccinesABC", test=2,
#cases= c(7,7,3,5,2,3,0),
#controls= 2,
#N_exposures= c(1200,1600,400,250,300,300,280),
#N_controls= 600,
#exposure_group= c("A","B","C","AB","AC","BC","ABC"),
#strata_group_cases= "n",
#strata_group_controls= "n",
#AlphaSpend= "n"
#                          )


}


