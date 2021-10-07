\name{ConfidenceInterval.Binomial}
\alias{ConfidenceInterval.Binomial}
\title{Confidence Interval for the Relative Risk Following a Sequential Test with Binary Data.}
\description{The function \code{ConfidenceInterval.Binomial} is used for constructing confidence interval for the relative risk in either continuous or group sequential analysis, or for a  combination of the two, on termination of the sequential surveillance. This function is useful, for example, in combination with the \code{Analyze.Binomial} function.}
\usage{
ConfidenceInterval.Binomial(Gamma=0.9,CV.lower="n",CV.upper="n",
GroupSizes,z="n",p="n",Cum.cases,Tailed="upper")
      }
\arguments{
\item{Gamma}{Confidence coefficient for the interval estimator of the relative risk. The default is Gamma=0.9.}
\item{CV.lower}{Signaling threshold for the evidence of "RR<R0", where R0 is a user-defined positive value when constructing this lower signaling threshold. It is given in the scale of the binomial cumulative data. Put NA for initial looks at the data when tests were not applicable. The default CV.lower="n" means that the sequential test was not designed to detect RR<R0 by the time of the confidence interval construction.}
\item{CV.upper}{Signaling threshold for evidence of "RR>R0", where R0 is a user-defined positive value when constructing this upper signaling threshold. It is given in the scale of the the binomial cumulative data. Put NA for initial looks at the data when tests were not applicable. The default CV.upper="n" means that the sequential test was not designed to detect RR>R0 by the time of the confidence interval construction.}
\item{GroupSizes}{Vector with the total number of events (cases+controls) between two looks at the data with regular and irregular group sizes. There is no default value.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case under the null hypothesis. If just a single number is given, then it will be used as a constant matching ratio for all groups. Otherwise, the dimension of z must coincide with the dimension of GroupSizes. The default z="n" means that the input p will be used instead.}
\item{p}{The probability of having a case under the null hypothesis H0:RR<=1. If just a single number is given, then it will be used as a constant probability for all groups. Otherwise, the dimension of p must coincide with the dimension of GroupSizes. The default p="n" means that the input z will be used instead.}
\item{Cum.cases}{Total number of cumulative cases on termination of the analysis. There is no default.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=R0, Tailed="lower" for H0:RR>=R0 or Tailed="two" for H0:R01<=RR<=R02. Important: R0, R01, and R02 are not parameters of this function. It is supposed that they were used when the user somehow constructed CV.lower and CV.upper.}
}
\details{
For continuous and group sequential analysis with binomial data, the confidence interval for the relative risk, RR, is constructed with \code{ConfidenceInterval.Binomial}.

For two-tailed testing (code{Tailed="two"}), both lower and upper signaling thresholds must be informed through
\code{CV.lower} and \code{CV.upper}. See details in Silva et al (2021).

\code{z} is a vector of positive numbers representing the matching ratios for each test (group). If a single number is given, then it will be used as a constant
matching ratio for all tests (groups). Otherwise, the dimension of \code{z} must coincide with the dimension of \code{GroupSizes}.
\code{z} represents the number of controls matched to each case. For example, if there are 3 controls matched to each case, \code{z}=3. 
In a self-control analysis, z is the ratio of the control interval to the risk interval.
For example, if the risk interval is 2 days long and the control interval is 7 days long, \code{z}=7/2.
In terms of \code{p}, the binomial probability under the null hypothesis, \code{p=1/(1+z)}, or equivalently, \code{z=1/p-1}.

Alternatively, instead of \code{z} the user can specify \code{p} directly.
Note that only one of these inputs, \code{z} or \code{p}, has to be specified, but if both are entered the code will only work if \code{z}
and \code{p} are such that \code{p=1/(1+z)}. Otherwise, an error message will appear to remind that such condition must be complied.
   

With \code{GroupSizes} the user informs the sample size of each subsequent test. Therefore, only positive integers are accepted in \code{GroupSizes}. 

The confidence interval is calculated by pivoting the probaility of rejecting the null hypothesis with one of the thresholds CV.lower or CV.upper. This is equivalent to inverting the two-tailed testing as described by Jennison and Turnbull (2000), page 181 of Section 8.5. See also Silva and Zhuang (2021) for more details on the calculations.

}
\value{
\item{RRl}{The lower limit of the confidence interval for RR.}
\item{RRu}{The upper limit of the confidence interval for RR.}
}
\author{Ivair Ramos Silva and Martin Kulldorf.
}
\references{

Jennison C, Turnbull B. (2000), Group Sequential Methods with Applications to Clinical Trials, London: Chapman and Hall/CRC.

Silva IR, Zhuang Y. (2021), Bounded-width Confidence Interval Following Optimal Sequential Test with Binary Data. Working paper.

Silva IR, Maro J, Kulldorff M. (2021), Exact sequential test for clinical trials and post-market drug and vaccine safety surveillance with Poisson and binary data. Statistics in Medicine, DOI: 10.1002/sim.9094.

}


\section{Acknowledgements}{
Development of the ConfidenceInterval.Binomial function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0,2.0 to 3.1).\cr
-     National Council of Scientific and Technological Development (CNPq), Brazil, process number 301391/2019-0. (v1.0).\cr
-     Research Support Foundation of the State of Minas Gerais (FAPEMIG), Brazil, grant number PQ-00787-21.\cr  
-	Federal University of Ouro Preto (UFOP), through contract under internal UFOP's resolution CEPE 4600 (v2.0 to 3.1).\cr
}

\keyword{Confidence Interval}

\section{See also}{
\code{\link[Sequential]{Performance.AlphaSpend.Binomial}}: for calculating signaling threshold for user-specified alpha spending with binomial data.\cr
\code{\link[Sequential]{CV.Binomial}}: for calculating Wald-type signaling thresholds for continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{Analyze.Binomial}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion with binomial data.
\code{\link[Sequential]{Analyze.wBinomial}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion with multiple weighted binary endpoints.\cr
}


\examples{


# ConfidenceInterval.Binomial(Gamma=0.9,CV.lower=c(NA,1,1,2),
# CV.upper=c(8,13,15,18), GroupSizes=c(8,5,5,6),z=1,
# Cum.cases=18,Tailed="two")

}


