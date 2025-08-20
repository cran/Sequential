\name{AnalyzeSetUp.Multinomial}
\alias{AnalyzeSetUp.Multinomial}
\title{Function to set up input parameters before using the \code{Analyze.Multinomial} function for the first time.}
\description{The function \code{AnalyzeSetUp.Multinomial} must be run ahead of \code{Analyze.Multinomial} in order to set up the sequential analysis before the first group of data is analyzed. The function obtains the main parameter settings and performs basic calculations that are necessary for the subsequent sequential analysis.}
\usage{
AnalyzeSetUp.Multinomial(name,N=200,alpha=0.05,AlphaSpendType=1,k,
R0=1,R1=2,rho=1,pmin=0.05, pmax=0.95,target_power=0.8,Rmin=1,Rmax=10,
gamma=0.9,m=100000,title="n",ExposuresNames="n",address="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and the same as the name given in the subsequent calls to the \code{Analyze.Multinomial} function. It cannot be the same as for another sequential analysis that is run simultaneously on the same computer. There is no default.}
\item{N}{The maximum sample size, at which the sequential analysis stops without rejecting the null hypothesis.}
\item{alpha}{The overall significance level. Must be in the range (0,0.5]. The default is "alpha=0.05".}
\item{AlphaSpendType}{Type of robust alpha spending. The possible values are 1 and 2 according to Silva and Maro(2025a). The default is equal to 1.}
\item{k}{Length of the multinomial vector.}
\item{R0}{Test margin under the null hypothesis. Must be a positive number. Default is R0=1.}
\item{R1}{Relative risk under the alternative hypothesis given the target_power.}
\item{rho}{The parameter rho is used to build the target alpha spending function according to a power-type function. The variable rho must be a positive number. The default value is "rho=1".}
\item{pmin}{Minimum value of the multinomial probability to activate its portion of alpha spending at each test. See Silva and Maro(2025) for more details. The default is pmin=0.05}
\item{pmax}{Maximum value of the multinomial probability to activate its portion of alpha spending at each test. See Silva and Maro(2025) for more details. The default is pmin=0.05}
\item{target_power}{The target power for detecting an elevated risk (>R0) for each entry of the multinomial vector.}
\item{Rmin}{Minimum value for the relative risks of adjacent exposures when constructing the confidence interval for RR of a given exposure.}
\item{Rmax}{Maximum value for the relative risks of adjacent exposures when constructing the confidence interval for RR of a given exposure.}
\item{gamma}{Confidence coefficient for the interval estimators of the relative risk of each multinomial entry. Default is 0.9.}
\item{m}{Monte Carlo replications of the multinomial for critical values calculations.}
\item{title}{Title for the results shown in the output tables and the illustrative graphics. It can be any text string. The default is that there is no title.}
\item{ExposuresNames}{Optional. This is to inform the name of the exposures related to each entry of the multinomial vector. For example, it can be c("A","B","AB") indicating that the counting entries in the "cases" vector are related to populations exposed to vaccines A, B, and AB, respectively.}
\item{address}{The address of the directory where the settings information of this sequential analysis is saved.}
}
\details{
The function \code{AnalyzeSetUp.Multinomial} has to be executed once, but just once, to set up the general statistical characteristics of the intended
sequential analysis, which is performed using the companion \code{Analyze.Multinomial} function. 

Sequential analysis methods are devoted to analyze data sets that accrue cumulatively over time, by conducting multiple statistical tests sequentially
as more data accrues. In such a setting, it is important to carefully plan the sequential analysis before the first data arrives. For example,
it is important to maintain certain analysis parameter values over time to avoid counting the same data twice, and to make sure that there
are no changes in the past data that has already been included in a prior test. To avoid these kinds of problems, the \code{AnalyzeSetUp.Multinomial}
function is used to set the analysis parameters a priori and to create a place to save the data as it accumulates over time. At the time of each sequential test,
this information is then automatically imported by the \code{Analyze.Multinomial} function, to ensure the correct concatenation of old and new information.

At each test, the function \code{Analyze.Multinomial} makes this concatenation automatically, but it will only work if the function \code{AnalyzeSetUp.Multinomial}
is executed before performing the very first test. 

When running \code{AnalyzeSetUp.Multinomial}, the user has the opportunity to choose the directory where the file with the general setup information and
the historical data are to be saved. Important: The location of this parameter and data file is saved in the temporary directory, so that directory cannot be cleaned until
the sequential analysis has been completed. Each sequential analysis needs a different identifier, which is set using the "name" parameter. Once a name is chosen,
it has to be written exactly the same way when running the function \code{Analyze.Multinomial}. 

The alpha spending plan is designed in a way to ensure that a target statistical power, gamma, is endured
for detecting a target increased risk, R1, for each multinomial entry. This is possible by means of a new concept introduced by
Silva and Maro(2025), the robust alpha spending plan. For more details about the method and the construction of the critical values,
see Silva and Maro(2025).

In addition to selecting the alpha spending plan, it is necessary to specify the overall alpha, or maximum Type I error probability, for the sequential analysis as a whole.
It is also necessary to specify the maximum length of the sequential analysis, N, so that the sequential analysis stops without rejecting the null hypothesis when
a total of N observations are obtained.  
  
}

\value{
\item{inputSetUp}{The \code{AnalyzeSetUp.Multinomial} function creates a data.frame with the main information concerning the tuning parameterization for the planned surveillance and the historical information about the performed tests.}
}

\author{ Ivair Ramos Silva, Judith Maro.
}

\section{Acknowledgements}{
Development of the AnalyzeSetUp.Multinomial function was funded by:\cr
- Department of Population Medicine, Harvard Medical School and Harvard Pilgrim Health Care Institute.\cr
- Conselho Nacional de Desenvolvimento Cientifico e Tecnologico (CNPq), Brazil, and Fundacao de Amparo a Pesquisa do Estado de Minas Gerais (FAPEMIG), Brazil.\cr
}

\references{
Silva IR, Maro J. (2025), Adaptive Sequential Multiple Hypotheses Testing for Seasonal Concomitant Vaccines Safety Surveillance. Working paper, Department of Population Medicine, Harvard Medical School and Harvard Pilgrim Health Care Institute.
}

\keyword{Multinomial sequential analysis}

\section{See also}{
\code{\link[Sequential]{Analyze.Multinomial}}: for running the sequential analysis that was set up using the \code{AnalyzeSetUp.Multinomial} function.
}


\examples{

# See example in the description of the Analyze.Multinomial function. 

}


