\name{Optimal.Binomial}
\alias{Optimal.Binomial}
\title{Optimal alpha spending for minimizing expected time to signal for continuous and group sequential analysis with binomial data.}
\description{The function \code{Optimal.Binomial} obtains the optimal alpha spending function that minimizes the expected time to signal for target statistical power and fixed relative risk, when doing continuous or group sequential analysis for binomial data.}
\usage{
Optimal.Binomial(Objective="ETimeToSignal",N="n",
z="n",p="n",alpha,power,RR,GroupSizes="n",Tailed= "upper")
      }
\arguments{
\item{Objective}{Statistical performance measure to minimize. Options are "ETimeToSignal", for minimizing expected time to signal, and "ESampleSize", for minimizing expected sample size. Default is "ETimeToSignal".}
\item{N}{The upper limit on the sample size (length of surveillance) expressed in terms of the total number of events (cases plus controls). "N" must be a positive integer. To avoid very large computation times, we restrict the usage of values greater than 240. Typically, this is not a major restriction. For example, for "RR=2.5", "alpha=0.05" and "z=1", the statistical power is approximately 1 for "N>=150". Default is 'n', which means that the optimal 'N' will also be delivered by this function.}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case under the null hypothesis. There is no default value.}
\item{p}{The probability of having a case under the null hypothesis. There is no default value.}
\item{alpha}{The significance level. The default value is "alpha=0.05". Must be in the range (0, 0.5].}
\item{power}{The target statistical power to detect an increased risk of the relative risk (RR). There is no default value.}
\item{RR}{A target relative risk to be detected with the requested statistical power.}
\item{GroupSizes}{Vector with the number of events (exposed+unexposed) between two looks at the data, i.e, irregular group sizes. Important: Must sums up N. The default 'n' means continuous sequential testing.}
\item{Tailed}{Defines between one-tailed and two-tailed testing. Possible entries are "lower", "upper", and "two".  Default is "upper".}
}
\details{
The function \code{Optimal.Binomial} elicits the optimal alpha spending for continuous and group binomial sequential testing
in order to provide the desired statistical power for a user-specified relative risk RR. The alpha spending provided minimizes the expected time to signal,
for Objective="ETimeToSignal", which is expected number of events when the null hypothesis is rejected, or expected length of surveillance if Objective="ESampleSize". 

The sample size, N, is given in the scale of the total number of events (cases+controls). The input z represents the number of controls matched to each case.
For example, if there are 3 controls matched to each case, "z=3". 

In a self-control analysis, z is the ratio of the control interval to the risk interval. For example, if the risk interval is 2 days long and the control
interval is 7 days long, z=7/2. In terms of p, the binomial probability under the null hypothesis, p=1/(1+z), or equivalently, z=1/p-1.
The parameter z must be a positive number.

Alternatively, instead of z the user can specify p directly.
Note that only one of these inputs, z or p, has to be specified, but if both are entered the code will only work if z and p are such that p=1/(1+z).
Otherwise, an error message will appear to remind that such condition must be complied.

The optimal alpha spending solution is obtained by means of linear programing, which is possible following the exact derivations introduced by Silva and Kulldorff (2018).
For the linear programing part, the code uses the function \code{simplex}.

\code{Optimal.Binomial} works for large sample sizes such as 300 in the continuous sequential fashion, but it can take very long time to run in such cases.
Thus, for continuous sequential analysis, the usage is restricted for N values of at most 240. The computation time for N=240 under continuous fashion can
take one day or more. But, for smaller values, like e.g. N=150, the execution time is around 2 hours.
For N=120 this time reduces to something around 10 minuts, and, for N<100, \code{Optimal.Binomial} takes only a few seconds to run.
But, processing time is much smaller for group sequential analysis. For example, take N=240. In this case, the optimum solution for two-stage (G=2) group
sequential analysis takes 16 minutes to run. The cases of N values of 200 or less, for G values of 10 or less, it will take just a few seconds to run.
These execution times were estimated using a regular PC(Windows 7, Intel(R) Core(TM) i7-2675QM CPU, 2.20GHz). 
 
}


\value{
\item{optimal_alpha_spending}{The optimal cumulative alpha spending. In case of Tailed="two", there are 'optimal_alpha_spending_lower' and 'optimal_alpha_spending_upper', for the lower and upper signaling threshold, respectively.}
\item{minTimeToSignal}{The minimum expected time to signal under RR. This is provided only if Objective="ETimeToSignal".}
\item{minESampleSize}{The minimum expected time to signal under RR. This is provided only if Objective="ESampleSize".}
\item{ETimeToSignal}{The expected time to signal associated to the alpha spending solution irrespectively to the content of the input 'Objective' under RR.}
\item{EsampleSize}{The expected sample size associated to the alpha spending solution irrespectively to the content of the input 'Objective' under RR.}
\item{Power}{Statistical power, obtained by usage of the optimal alpha spending with relative risk equal to RRtrue.}
\item{solved}{Logical variable. It is equal to 1 if the linear programing procedure reached the optimal solution, 0 for inconclusive solution, and -1 if the method fails to find the solution.}
\item{Optimal_N}{The maximum length of surveillance at which the optimal alpha spending reaches 'alpha' by a precision of 10^-6. This is returned only for N="n".}
}

\author{Ivair Ramos Silva, Martin Kulldorff.
}
\references{
Silva, I.R. and Kulldorff, M. (2018). Optimal Alpha Spending for Sequential Analysis With Binary Data. Working paper.
}
\section{Acknowledgements}{
Development of the Optimal.Binomial function was funded by:\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999.
}

\keyword{Optimal binomial alpha spending}

\section{See also}{
\code{\link[Sequential]{CV.Binomial}}: for calculating critical values in continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{Performance.Binomial}}: for calculating the statistical power, expected time to signal and expected time of analysis for continuous sequential analysis with binomial data.\cr
\code{\link[Sequential]{SampleSize.Poisson}}: sample size calculation for continuous sequential analysis with Poisson data.  
}


\examples{

#system.time(resESS<- Optimal.Binomial(Objective="ESampleSize",
# N=120,z=1,p="n",alpha=0.05,power=0.9,RR=2,GroupSizes="n",Tailed= "upper"))

}


