\name{AnalyzeSetUp.Poisson}
\alias{AnalyzeSetUp.Poisson}
\title{Function to set up input parameters before using the \code{Analyze.Poisson} function for the first time.}
\description{The function \code{AnalyzeSetUp.Poisson} must be run ahead of \code{Analyze.Poisson} in order to set up the sequential analysis before the first group of data is analyzed. The function obtains the main parameter settings and performs basic calculations that are necessary for the subsequent sequential analysis.}
\usage{
AnalyzeSetUp.Poisson(name,SampleSize,alpha=0.05,M=1,
AlphaSpendType="Wald",rho="n",title="n",address="n",Tailed="upper")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and the same as the name given in the subsequent calls to the \code{Analyze.Poisson} function. It cannot be the same as for another sequential analysis that is run simultaneously on the same computer. There is no default.}
\item{SampleSize}{The maximum length of surveillance at which the sequential analysis stops without rejecting the null hypothesis. It is defined in terms of the expected sample size under the null hypothesis. There is no default value.}
\item{alpha}{The overall significance level. Must be in the range (0,0.5]. The default is "alpha=0.05".}
\item{M}{The minimum number of events required before the null hypothesis can be rejected. It must be a positive integer. The default value is "M=1".}
\item{AlphaSpendType}{The type of alpha spending function to be used. The options are AlphaSpendType="Wald" and AlphaSpendType="power-type". With the 'Wald' option, the Wald type upper rejection boundary is used, which is flat with respect to the likelihood ratio. With the power-type option, the alpha spending uses a power function with parameter rho, with rho defined by the user. This alpha spending setting is automatically used when the \code{Analyze.Poisson} function is run, but, during the sequential analysis, and before each test, the user can always specify an arbitrary amount of alpha spending to be used up until and including that test. See below for details.}
\item{rho}{The parameter rho is used to build the target alpha spending function according to a power-type function. See below for details. It is not used for other alpha spending options. The variable rho must be a positive number. The default value is "rho=1".}
\item{title}{Title for the results shown in the output tables and the illustrative graphics. It can be any text string. The default is that there is no title.}
\item{address}{The address of the directory where the settings information of this sequential analysis is saved.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
}


\details{
The function \code{AnalyzeSetUp.Poisson} has to be executed once, but just once, to set up the general statistical characteristics of the intended
sequential analysis, which is performed using the companion \code{Analyze.Poisson} function. 

Sequential analysis methods are devoted to analyze data sets that accrue cumulatively over time, by conducting multiple statistical tests sequentially
as more data accrues. In such a setting, it is important to carefully plan the sequential analysis before the first data arrives. For example,
it is important to maintain certain analysis parameter values over time to avoid counting the same data twice, and to make sure that there
are no changes in the past data that has already been included in a prior test. To avoid these kinds of problems, the \code{AnalyzeSetUp.Poisson}
function is used to set the analysis parameters a priori and to create a place to save the data as it accumulates over time. At the time of each sequential test,
this information is then automatically imported by the \code{Analyze.Poisson} function, to ensure the correct concatenation of old and new information.

At each test, the function \code{Analyze.Poisson} makes this concatenation automatically, but it will only work if the function \code{AnalyzeSetUp.Poisson}
is executed before performing the very first test. 

When running \code{AnalyzeSetUp.Poisson}, the user has to choose the directory where the file with the general setup information and
the historical data are to be saved. This step is mandatory and error messages are reported if non-valid address is informed. Important: The location of this parameter and data file is saved in the temporary directory, so that directory cannot be cleaned until
the sequential analysis has been completed. Each sequential analysis needs a different identifier, which is set using the "name" parameter. Once a name is chosen,
it has to be written exactly the same way when running the function \code{Analyze.Poisson}. 

\code{AnalyzeSetUp.Poisson} and \code{Analyze.Poisson} work for different types of alpha spending plans (\eqn{F(t)}). One option is to use the classical
Wald type upper rejection boundary, which is flat with respect to the likelihood function. This is the same boundary used by the \code{CV.Poisson} and
\code{CV.G.Poisson} functions.

Another alpha spending option is the power-type alpha spending plan (Kim and DeMetz 1987, p150; Jennison and Turnbull 2000, p148), with parameter rho: \eqn{F(t)= alpha*t^{rho}}, 
where \eqn{alpha} is the overall significance level and \eqn{t} is a fraction of SampleSize, the maximum length of sequential analysis.
According to Silva (2018), the choice 'rho=1' is indicated when minimization of expected time to signal is a design criterion,
which is then the default in \code{AnalyzeSetUp.Poisson.}. 

}

\value{
\item{inputSetUp}{The \code{AnalyzeSetUp.Poisson} function creates a data.frame with the main information concerning the tuning parameterization for the planned surveillance and the historical information about the performed tests. The 'inputSetUp' data.frame is used by \code{Analyze.Poisson}, then it must be available when running \code{Analyze.Poisson}, but there is no need to manually look at it.}
}

\author{ Ivair Ramos Silva, Martin Kulldorff.
}

\section{Acknowledgements}{
Development of the AnalyzeSetUp.Poisson function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999.
}

\references{
Jennison C, Turnbull B. (2000), Group Sequential Methods with Applications to Clinical
Trials, \emph{no. ISBN 0-8493-0316-8, London: Chapman and Hall/CRC}.

Kim K, DeMets DL. (1987), Design and Analysis of Group Sequential Tests Based on the Type I Error Spending Rate Function. Biometrika, \bold{74}, n.1: 149--154.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

Kulldorff M, Silva IR. (2015). Continuous Post-market Sequential Safety Surveillance with Minimum Events to Signal. REVSTAT Statistical Journal, 15(3): 373--394.

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva IR. (2018). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance With Poisson Data. Methodol Comput Appl Probab, 20(2), 739-750.
}

\keyword{Poisson sequential analysis}

\section{See also}{
\code{\link[Sequential]{Analyze.Poisson}}: for running the sequential analysis that was set up using the \code{AnalyzeSetUp.Poisson} function.
}


\examples{

# See example in the description of the Analyze.Poisson function. 

}


