\name{AnalyzeSetUp.Binomial}
\alias{AnalyzeSetUp.Binomial}
\title{Function to Set Up the Input Parameters Before Using the Analyze.Binomial Function for the First Time.}
\description{The function \code{AnalyzeSetUp.Binomial} must be run ahead of \code{Analyze.Binomial} in order to set up the sequential analysis before the first group of data is analyzed. The function obtains the main parameter settings and performs basic calculations that are necessary for the subsequent sequential analysis.}
\usage{
AnalyzeSetUp.Binomial(name,N,alpha=0.05,zp="n",pp="n",M=1,AlphaSpendType="Wald",rho="n",
title="n",address="n")
      }
\arguments{
\item{name}{The name of the sequential analysis. Must be identical for all looks at the data, and the same as the name given in the subsequent calls to the \code{Analyze.Binomial} function. It cannot be the same as for another sequential analysis that is run simultaneously on the same computer. There is no default.}
\item{N}{The maximum sample size, at which the sequential analysis stops without rejecting the null hypothesis. There is no default value.}
\item{alpha}{The overall significance level. Must be in the range (0,0.5]. The default is "alpha=0.05".}
\item{zp}{The prediction for z, the expected ratio between cases and controls under the null hypothesis that will be specified in the \code{Analyze.Binomial} function. This variable is only needed when AlphaSpendType= "Wald", and it is used to calculate the appropriate rejection boundary. If the z used in \code{Analyze.Binomial} during the actual sequential analysis is different from zp, that is okay, and the sequential analysis will still maintain the correct alpha level. The default value is "z=1".}
\item{pp}{The prediction for p, the expected probability under the null hypothesis that will be specified in the \code{Analyze.Binomial} function. This variable is only needed when AlphaSpendType= "Wald", and it is used to calculate the appropriate rejection boundary. If the p used in \code{Analyze.Binomial} during the actual sequential analysis is different from pp, that is okay, and the sequential analysis will still maintain the correct alpha level. There is no default value.}
\item{M}{The minimum number of events required before the null hypothesis can be rejected. It must be a positive integer. The default value is "M=1".}
\item{AlphaSpendType}{The type of alpha spending function to be used. The options are AlphaSpendType="Wald" and AlphaSpendType="power-type". With the 'Wald' option, the Wald type upper rejection boundary is used, which is flat with respect to the likelihood ratio. With the power-type option, the alpha spending uses a power function with parameter rho, with rho defined by the user. This alpha spending setting is automatically used when the \code{Analyze.Binomial} function is run, but, during the sequential analysis, and before each test, the user can always specify an arbitrary amount of alpha spending to be used up until and including that test. See below for details.}
\item{rho}{The parameter rho is used to build the target alpha spending function according to a power-type function. See below for details. It is not used for other alpha spending options. The variable rho must be a positive number. The default value is "rho=1".}
\item{title}{Title for the results shown in the output tables and the illustrative graphics. It can be any text string. The default is that there is no title.}
\item{address}{The address of the directory where the settings information of this sequential analysis is saved.}
}


\details{
The function \code{AnalyzeSetUp.Binomial} has to be executed once, but just once, to set up the general statistical characteristics of the intended
sequential analysis, which is performed using the companion \code{Analyze.Binomial} function. 

Sequential analysis methods are devoted to analyze data sets that accrue cumulatively over time, by conducting multiple statistical tests sequentially
as more data accrues. In such a setting, it is important to carefully plan the sequential analysis before the first data arrives. For example,
it is important to maintain certain analysis parameter values over time to avoid counting the same data twice, and to make sure that there
are no changes in the past data that has already been included in a prior test. To avoid these kinds of problems, the \code{AnalyzeSetUp.Binomial}
function is used to set the analysis parameters a priori and to create a place to save the data as it accumulates over time. At the time of each sequential test,
this information is then automatically imported by the \code{Analyze.Binomial} function, to ensure the correct concatenation of old and new information.

At each test, the function \code{Analyze.Binomial} makes this concatenation automatically, but it will only work if the function \code{AnalyzeSetUp.Binomial}
is executed before performing the very first test. 

When running \code{AnalyzeSetUp.Binomial}, the user has the opportunity to choose the directory where the file with the general setup information and
the historical data are to be saved. Alternatively, the user can click in the "Cancel" button, in which case this file will be saved in the temporary
directory. Important: The location of this parameter and data file is saved in the temporary directory, so that directory cannot be cleaned until
the sequential analysis has been completed. Each sequential analysis needs a different identifier, which is set using the "name" parameter. Once a name is chosen,
it has to be written exactly the same way when running the function \code{Analyze.Binomial}. 

\code{AnalyzeSetUp.Binomial} and \code{Analyze.Binomial} works for different types of alpha spending plans (\eqn{F(t)}). One option is to use the classical
Wald type upper rejection boundary, which is flat with respect to the likelihood function. This is the same boundary used by the \code{CV.Binomial} and
\code{CV.G.Binomial} functions. In order to use this boundary, one should pre-specify the binomial probability p under the null hypothesis, or,
equivalently, the ratio \eqn{z=1/p-1}, which is the number of controls matched to each case in a matched analysis. For example, if the probability of having
a case (instead of a control) is \eqn{p = 1=(1 + z) = 0.5}, then we have "z=1" (1:1 matching ratio), and, if p = 0.25, we have "z=3" (1:3 matching ratio).

In \code{AnalyzeSetUp.Binomial}, the predicted z is specified (the input zp), but if it turns out that the actual z is different, that is okay,
since the actual z that is specified in \code{Analyze.Binomial} does not have to be the same as the predicted zp that is specified
in \code{AnalyzeSetUp.Binomial}. The latter is only used to set the alpha spending plan. The former, the actual z, is used to calculate the
likelihood function which in turn determines whether the null hypothesis should be rejected or not. If the actual z is variable, so that
it is different for different observations, we recommend setting the predicted z to be our best guess about the average of the actual zs.
Alternatively, instead of zp the user can specify pp, the best guess about the average of the actual ps.
Note that only one of these parameters has to be specified, but if both are entered the code will only work if zp and pp are such that pp=1/(1+zp).
Otherwise, an error message will appear to remind that such condition must be complied. 

Another alpha spending option is the power-type alpha spending plan (Kim and DeMetz 1987, p150; Jennison and Turnbull 2000, p148), with parameter rho: \eqn{F(t)= alpha*t^{rho}}, 
where \eqn{alpha} is the overall significance level and \eqn{t} is a fraction of N, the maximum length of sequential analysis.

In addition to selecting the alpha spending plan, it is necessary to specify the overall alpha, or maximum Type I error probability, for the sequential analysis as a whole.
It is also necessary to specify the maximum length of the sequential analysis, N, so that the sequential analysis stops without rejecting the null hypothesis when
N total observations have been obtained.  
}

\value{
\item{inputSetUp}{The \code{AnalyzeSetUp.Binomial} function creates a data.frame with the main information concerning the tuning parameterization for the planned surveillance and the historical information about the performed tests. The 'inputSetUp' data.frame is used by \code{Analyze.Binomial}, then it must be available when running \code{Analyze.Binomial}, but there is no need to manually look at it.}
}

\author{ Ivair Ramos Silva, Martin Kulldorff.
}

\section{Acknowledgements}{
Development of the AnalyzeSetUp.Binomial function was funded by:\cr
-	Food and Drug Administration, Center for Drug Evaluation and Research, through the Mini-Sentinel Project (base version, documentation);\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (user defined alpha spending functions, improved documentation);\cr
We thank Claudia Coronel-Moreno for valuable editorial support.
}

\references{
Jennison C, Turnbull B. (2000), Group Sequential Methods with Applications to Clinical
Trials, \emph{no. ISBN 0-8493-0316-8, London: Chapman and Hall/CRC}.

Kim K, DeMets DL. (1987), Design and Analysis of Group Sequential Tests Based on the Type I Error Spending Rate Function. Biometrika, \bold{74}, n.1: 149--154.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

Kulldorff M, Silva IR. (2015). Continuous post-market sequential safety surveillance with minimum events to signal. arxiv:1503.01978 [stat.ap].

Silva IR, Kulldorff M. (2015). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, in press.
}

\keyword{Binomial sequential analysis}

\section{See also}{
\code{\link[Sequential]{Analyze.Binomial}}: for running the sequential analysis that was set up using the \code{AnalyzeSetUp.Binomial} function.
}


\examples{

# See example in the description of the Analyze.Binomial function. 

}


