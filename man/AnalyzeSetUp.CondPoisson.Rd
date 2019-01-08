\name{AnalyzeSetUp.CondPoisson}
\alias{AnalyzeSetUp.CondPoisson}
\title{Function to set up input parameters before using the \code{Analyze.CondPoisson} function for the first time.}
\description{The function \code{AnalyzeSetUp.CondPoisson} must be run ahead of \code{Analyze.CondPoisson} in order to set up the sequential analysis before the first group of data is analyzed. The function obtains the main parameter settings and performs basic calculations that are necessary for the subsequent sequential analysis.}
\usage{
AnalyzeSetUp.CondPoisson(name,SampleSizeType="Events",T="n",K="n",cc,
alpha=0.05,M=1, AlphaSpendType="Wald",rho="n",title="n",address="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and the same as the name given in the subsequent calls to the \code{Analyze.CondPoisson} function. It cannot be the same as for another sequential analysis that is run simultaneously on the same computer. There is no default.}
\item{SampleSizeType}{It is a string specifying the scale of the maximum sample size at which the sequential analysis stops without rejecting the null hypothesis. The only two possibilities are "SampleSizeType=PersonTimeRatio" or "SampleSizeType=Events". The default is "SampleSizeType=Events". See details.}
\item{T}{Maximum sample size defined in the scale of the ratio between surveillance and historical person-time. This only produces effects when "SampleSizeType=PersonTimeRatio".}
\item{K}{Maximum sample size defined in the scale of the number of events observed in the surveillance period. This only produces effects when "SampleSizeType=Events".}
\item{cc}{Number of events observed in the historical period.}
\item{alpha}{The overall significance level. Must be in the range (0,0.5]. The default is "alpha=0.05".}
\item{M}{The minimum number of events required before the null hypothesis can be rejected. It must be a positive integer. The default value is "M=1".}
\item{AlphaSpendType}{The type of alpha spending function to be used. The options are AlphaSpendType="Wald" and AlphaSpendType="power-type". With the 'Wald' option, the Wald type upper rejection boundary is used, which is flat with respect to the likelihood ratio. With the power-type option, the alpha spending uses a power function with parameter rho, with rho defined by the user. This alpha spending setting is automatically used when the \code{Analyze.CondPoisson} function is run, but, during the sequential analysis, and before each test, the user can always specify an arbitrary amount of alpha spending to be used up until and including that test. See below for details.}
\item{rho}{The parameter rho is used to build the target alpha spending function according to a power-type function. See below for details. It is not used for other alpha spending options. The variable rho must be a positive number. The default value is "rho=1".}
\item{title}{Title for the results shown in the output tables and the illustrative graphics. It can be any text string. The default is that there is no title.}
\item{address}{The address of the directory where the settings information of this sequential analysis is saved.}
}


\details{
The function \code{AnalyzeSetUp.CondPoisson} has to be executed once, but just once, to set up the general statistical characteristics of the intended
sequential analysis, which is performed using the companion \code{Analyze.CondPoisson} function. 

Sequential analysis methods are devoted to analyze data sets that accrue cumulatively over time, by conducting multiple statistical tests sequentially
as more data accrues. In such a setting, it is important to carefully plan the sequential analysis before the first data arrives. For example,
it is important to maintain certain analysis parameter values over time to avoid counting the same data twice, and to make sure that there
are no changes in the past data that has already been included in a prior test. To avoid these kinds of problems, the \code{AnalyzeSetUp.CondPoisson}
function is used to set the analysis parameters a priori and to create a place to save the data as it accumulates over time. At the time of each sequential test,
this information is then automatically imported by the \code{Analyze.CondPoisson} function, to ensure the correct concatenation of old and new information.

At each test, the function \code{Analyze.CondPoisson} makes this concatenation automatically, but it will only work if the function \code{AnalyzeSetUp.CondPoisson}
is executed before performing the very first test. 

When running \code{AnalyzeSetUp.CondPoisson}, the user has to choose the directory where the file with the general setup information and
the historical data are to be saved. This step is mandatory and error messages are reported if a non-valid address is informed. Important:
The location of this parameter and data file is saved in the temporary directory, so that directory cannot be cleaned until
the sequential analysis has been completed. Each sequential analysis needs a different identifier, which is set using the "name" parameter. Once a name is chosen,
it has to be written exactly in the same way when running the function \code{Analyze.CondPoisson}. 

\code{AnalyzeSetUp.CondPoisson} and \code{Analyze.CondPoisson} work for different types of alpha spending plans (\eqn{F(t)}). One option is to use the classical
Wald type upper rejection boundary, which is flat with respect to the likelihood function. This is the same boundary used by the \code{CV.CondPoisson} function.

Another alpha spending option is the power-type alpha spending plan (Kim and DeMetz 1987, p150; Jennison and Turnbull 2000, p148), with parameter rho:
\eqn{F(t)= alpha*t^{rho}}, 
where \eqn{alpha} is the overall significance level and \eqn{t} is a fraction of SampleSize, the maximum length of sequential analysis.

Attention is required for the input parameter "SampleSizeType". With this parameter, the user can choose the scale of the maximum sample size
at which the surveillance is stoped without rejecting the null hypothesis. The idea of having two options for defining the scale of the maximum sample size,
"SampleSizeType=PersonTimeRatio" or "SampleSizeType=Events",
was introduced by Silva et al. (2016). With SampleSizeType="PersonTimeRatio", the upper limit on the time of surveillance is given in the scale of the ratio
between the cumulative person-time from the surveillance data up to the kth event, \eqn{P_{k}} in the notation of Silva et al. (2016), by the person-time from
the historical data, which was denoted with \eqn{V} by Silva et al. (2016).

If SampleSizeType="PersonTimeRatio", then the user has to inform a positive value for the input parameter "T". Usually, choices between 2 and 5 are adequate.
With SampleSizeType="Events", the upper limit is given in the scale of the number of events observed during the surveillance, and hence the user must
specify a positive integer for the input parameter "K", such as e.g. "K=50" or "K=150". For more details on the exact calculations of critical values
implemented in this package, and of all the other statistical performance measures also available, see the work of Silva et al. (2016). 

}

\value{
\item{inputSetUp}{The \code{AnalyzeSetUp.CondPoisson} function creates a data.frame with the main information concerning the tuning parameterization for the planned surveillance and the historical information about the performed tests. The 'inputSetUp' data.frame is used by \code{Analyze.CondPoisson}, then it must be available when running \code{Analyze.CondPoisson}, but there is no need to manually look at it.}
}

\author{ Ivair Ramos Silva, Martin Kulldorff.
}

\section{Acknowledgements}{
Development of the AnalyzeSetUp.CondPoisson function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999.
- Foundation for Research Support of Minas Gerais State (FAPEMIG), MG, Brazil, through the grant Demanda Universal.
}

\references{
Jennison C, Turnbull B. (2000), Group Sequential Methods with Applications to Clinical
Trials, \emph{no. ISBN 0-8493-0316-8, London: Chapman and Hall/CRC}.

Kim K, DeMets DL. (1987), Design and Analysis of Group Sequential Tests Based on the Type I Error Spending Rate Function. Biometrika, \bold{74}, n.1: 149--154.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

KKulldorff M, Silva IR. (2015). Continuous Post-market Sequential Safety Surveillance with Minimum Events to Signal. REVSTAT Statistical Journal, 15(3): 373--394.

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva IR, Li L, Kulldorff M. (2016). Exact Conditional Sequential Testing for Poisson Data. Working paper.
}

\keyword{Poisson sequential analysis}

\section{See also}{
\code{\link[Sequential]{Analyze.CondPoisson}}: for running the sequential analysis that was set up using the \code{AnalyzeSetUp.CondPoisson} function.
}


\examples{

# See example in the description of the Analyze.CondPoisson function. 

}


