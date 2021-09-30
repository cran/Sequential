\name{AnalyzeSetUp.wBinomial}
\alias{AnalyzeSetUp.wBinomial}
\title{Function to set up input parameters before using the \code{Analyze.wBinomial} function for the first time.}
\description{The function \code{AnalyzeSetUp.wBinomial} must be run ahead of \code{Analyze.wBinomial} in order to set up the sequential analysis before the first group of data is analyzed. The function obtains the main parameter settings and performs basic calculations that are necessary for the subsequent sequential analysis.}
\usage{
AnalyzeSetUp.wBinomial(name,N,alpha=0.05,M=1,rho=0.5,
title="n",address="n",Tailed="upper")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and the same as the name given in the subsequent calls to the \code{Analyze.wBinomial} function. It cannot be the same as for another sequential analysis that is run simultaneously on the same computer. There is no default.}
\item{N}{The maximum sample size, at which the sequential analysis stops without rejecting the null hypothesis.}
\item{alpha}{The overall significance level. Must be in the range (0,0.5]. The default is "alpha=0.05".}
\item{M}{The minimum number of events required before the null hypothesis can be rejected. It must be a positive integer. The default value is "M=1".}
\item{rho}{The parameter rho is used to build the target alpha spending function according to a power-type function. See below for details. It is not used for other alpha spending options. The variable rho must be a positive number. The default value is "rho=0.5".}
\item{title}{Title for the results shown in the output tables and the illustrative graphics. It can be any text string. The default is that there is no title.}
\item{address}{The address of the directory where the settings information of this sequential analysis is saved.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="two" for H0:RR=1.}
}


\details{
The function \code{AnalyzeSetUp.wBinomial} has to be executed once, but just once, to set up the general statistical characteristics of the intended
sequential analysis, which is performed using the companion \code{Analyze.wBinomial} function. 

Sequential analysis methods are devoted to analyze data sets that accrue cumulatively over time, by conducting multiple statistical tests sequentially
as more data accrues. In such a setting, it is important to carefully plan the sequential analysis before the first data arrives. For example,
it is important to maintain certain analysis parameter values over time to avoid counting the same data twice, and to make sure that there
are no changes in the past data that has already been included in a prior test. To avoid these kinds of problems, the \code{AnalyzeSetUp.wBinomial}
function is used to set the analysis parameters a priori and to create a place to save the data as it accumulates over time. At the time of each sequential test,
this information is then automatically imported by the \code{Analyze.wBinomial} function, to ensure the correct concatenation of old and new information.

At each test, the function \code{Analyze.wBinomial} makes this concatenation automatically, but it will only work if the function \code{AnalyzeSetUp.wBinomial}
is executed before performing the very first test. 

When running \code{AnalyzeSetUp.wBinomial}, the user has the opportunity to choose the directory where the file with the general setup information and
the historical data are to be saved. Important: The location of this parameter and data file is saved in the temporary directory, so that directory cannot be cleaned until
the sequential analysis has been completed. Each sequential analysis needs a different identifier, which is set using the "name" parameter. Once a name is chosen,
it has to be written exactly the same way when running the function \code{Analyze.wBinomial}. 

\code{AnalyzeSetUp.wBinomial} and \code{Analyze.wBinomial} works for different types of alpha spending plans (\eqn{F(t)}).
The alpha spending option is the power-type alpha spending plan (Kim and DeMetz 1987, p150; Jennison and Turnbull 2000, p148), with parameter rho: \eqn{F(t)= alpha*t^{rho}}, 
where \eqn{alpha} is the overall significance level and \eqn{t} is a fraction of N, the maximum length of sequential analysis.
According to Silva (2018), 'rho=0.5' is indicated when expected time to signal is the design criterion, hence this is the default in \code{AnalyzeSetUp.wBinomial}. 

In addition to selecting the alpha spending plan, it is necessary to specify the overall alpha, or maximum Type I error probability, for the sequential analysis as a whole.
It is also necessary to specify the maximum length of the sequential analysis, N, so that the sequential analysis stops without rejecting the null hypothesis when
a total of N observations are obtained.  
}

\value{
\item{inputSetUp}{The \code{AnalyzeSetUp.wBinomial} function creates a data.frame with the main information concerning the tuning parameterization for the planned surveillance and the historical information about the performed tests. The 'inputSetUp' data.frame is used by \code{Analyze.wBinomial}, then it must be available when running \code{Analyze.wBinomial}, but there is no need to manually look at it.}
}

\author{ Ivair Ramos Silva, Martin Kulldorff.
}

\section{Acknowledgements}{
Development of the AnalyzeSetUp.Binomial function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (user defined alpha spending functions, improved documentation);\cr
}

\references{
Jennison C, Turnbull B. (2000), Group Sequential Methods with Applications to Clinical
Trials, \emph{no. ISBN 0-8493-0316-8, London: Chapman and Hall/CRC}.

Kim K, DeMets DL. (1987), Design and Analysis of Group Sequential Tests Based on the Type I Error Spending Rate Function. Biometrika, \bold{74}, n.1: 149--154.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

Kulldorff M, Silva IR. (2015). Continuous post-market sequential safety surveillance with minimum events to signal. arxiv:1503.01978 [stat.ap].

Silva IR, Gagne J, Najafzadeh M, Kulldorff M. (2020). Exact Sequential Analysis for Multiple Weighted Binomial Endpoints. Statistics in Medicine, 39(3), 340--351.

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva, IR and Kulldorff, M. and Yih, W. Katherine (2020). Optimal alpha spending for sequential analysis with binomial data. Journal of the Royal Statistical Society Series B, 82(4) p. 1141--1164.

Silva IR. (2018). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance with Binomial Data. Statistics in Medicine, 15;37(1), 107-118.

}

\keyword{Binomial sequential analysis}

\section{See also}{
\code{\link[Sequential]{Analyze.wBinomial}}: for running the sequential analysis that was set up using the \code{AnalyzeSetUp.wBinomial} function.
}


\examples{

# See example in the description of the Analyze.wBinomial function. 

}


