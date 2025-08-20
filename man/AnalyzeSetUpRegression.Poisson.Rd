\name{AnalyzeSetUpRegression.Poisson}
\alias{AnalyzeSetUpRegression.Poisson}
\title{Function to set up input parameters before using the \code{AnalyzeRegression.Poisson} function for the first time.}
\description{The function \code{AnalyzeSetUpRegression.Poisson} must be run ahead of \code{AnalyzeRegression.Poisson} in order to set up the sequential analysis before the first group of data is analyzed. The function obtains the main parameter settings and performs basic calculations that are necessary for the subsequent sequential analysis.}
\usage{
AnalyzeSetUpRegression.Poisson(name,N="n",alpha=0.05,
R0=1,rho=1,mref=999,title="n",address="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and the same as the name given in the subsequent calls to the \code{AnalyzeRegression.Poisson} function. It cannot be the same as for another sequential analysis that is run simultaneously on the same computer. There is no default.}
\item{N}{The maximum sample size, at which the sequential analysis stops without rejecting the null hypothesis.}
\item{alpha}{The overall significance level. Must be in the range (0,0.5]. The default is "alpha=0.05".}
\item{R0}{A positive real-valued number for the relative risk under R<=R0. Default is R0=1.}
\item{rho}{The parameter rho is used to build the target alpha spending function according to a power-type function. See below for details. The variable rho must be a positive number. The default value is "rho=1".}
\item{mref}{The parameter mref is the minimum number of Monte Carlo replications to calculate the p-value per test following the method by Silva at al(2025).}
\item{title}{Title for the results shown in the output tables and the illustrative graphics. It can be any text string. The default is that there is no title.}
\item{address}{The address of the directory where the settings information of this sequential analysis is saved.}
}




\details{
The function \code{AnalyzeSetUpRegression.Poisson} has to be executed once, but just once, to set up the general statistical characteristics of the intended
sequential analysis, which is performed using the companion \code{AnalyzeRegression.Poisson} function. 

Sequential analysis methods are devoted to analyze data sets that accrue cumulatively over time, by conducting multiple statistical tests sequentially
as more data accrues. In such a setting, it is important to carefully plan the sequential analysis before the first data arrives. For example,
it is important to maintain certain analysis parameter values over time to avoid counting the same data twice, and to make sure that there
are no changes in the past data that has already been included in a prior test. To avoid these kinds of problems, the \code{AnalyzeSetUpRegression.Poisson}
function is used to set the analysis parameters a priori and to create a place to save the data as it accumulates over time. At the time of each sequential test,
this information is then automatically imported by the \code{AnalyzeRegression.Poisson} function, to ensure the correct concatenation of old and new information.

At each test, the function \code{AnalyzeRegression.Poisson} makes this concatenation automatically, but it will only work if the function \code{AnalyzeSetUpRegression.Poisson}
is executed before performing the very first test. 

When running \code{AnalyzeSetUpRegression.Poisson}, the user has the opportunity to choose the directory where the file with the general setup information and
the historical data are to be saved. Important: The location of this parameter and data file is saved in the temporary directory, so that directory cannot be cleaned until
the sequential analysis has been completed. Each sequential analysis needs a different identifier, which is set using the "name" parameter. Once a name is chosen,
it has to be written exactly the same way when running the function \code{AnalyzeRegression.Poisson}. 

\code{AnalyzeSetUpRegression.Poisson} and \code{AnalyzeRegression.Poisson} work for the power-type alpha spending plan (Kim and DeMetz 1987, p150; Jennison and Turnbull 2000, p148), with parameter rho: \eqn{F(t)= alpha*t^{rho}}, 
where \eqn{alpha} is the overall significance level and \eqn{t} is a fraction of N, the maximum length of sequential analysis.
According to Silva (2018), 'rho=0.5' is indicated when expected time to signal is the design criterion. The choice 'rho=1' is a proper option that provides
a balance between time to signal and maximum lenght and surveillance, hence this is the default in \code{AnalyzeSetUpRegression.Poisson}. 

In addition to selecting the alpha spending plan, it is necessary to specify the overall alpha, or maximum Type I error probability, for the sequential analysis as a whole.
It is also necessary to specify the maximum length of the sequential analysis, N, so that the sequential analysis stops without rejecting the null hypothesis when
a total of N observations are obtained.  
  
}

\value{
\item{inputSetUp}{The \code{AnalyzeSetUpRegression.Poisson} function creates four data.frames (Decision_table, Relative_risk_estimates, Coefficients, Confidence_Intervals) with the main information concerning the tuning parameterization for the planned surveillance and the historical information about the performed tests.
The 'inputSetUp' list is used by \code{AnalyzeRegression.Poisson}, then it must be available when running \code{AnalyzeRegression.Poisson}, but there is no need to manually look at it.}
}

\author{ Ivair Ramos Silva.
}

\section{Acknowledgements}{
Development of the AnalyzeSetUpRegression.Poisson function was funded by:\cr
-	Food and Drug Administration, Center for Drug Evaluation and Research, through the Mini-Sentinel Project (base version, documentation);\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (user defined alpha spending functions, improved documentation);\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil, process number 302882/2022-7. (v4.3.1).\cr
-	Support Foundation to Minas Gerais State Research-Fapemig, Brazil, grant numbers PQ-00787-21 and RED-00133-21. (v3.1 to v4.3). \cr
}

\references{
Jennison C, Turnbull B. (2000), Group Sequential Methods with Applications to Clinical
Trials, \emph{no. ISBN 0-8493-0316-8, London: Chapman and Hall/CRC}.

Kim K, DeMets DL. (1987), Design and Analysis of Group Sequential Tests Based on the Type I Error Spending Rate Function. Biometrika, \bold{74}, n.1: 149--154.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva IR. (2018b), Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance With Poisson Data. Methodol Comput Appl Probab, 20(2), 739--750.

Silva IR, Montalban, J., Oliveira, F. (2025), Maximized Sequential Probability Ratio Test Regression. Working paper - Sentinel (TIDE) project, Department of Population Medicine, Harvard Medical School and Harvard Pilgrim Health Care Institute.
}

\keyword{Poisson sequential analysis}

\section{See also}{
\code{\link[Sequential]{AnalyzeRegression.Poisson}}: for running the sequential analysis that was set up using the \code{AnalyzeSetUpRegression.Poisson} function.
}


\examples{

# See example in the description of the AnalyzeRegression.Poisson function. 

}


