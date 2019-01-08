\name{CV.G.Binomial}
\alias{CV.G.Binomial}
\title{Calculates exact critical values for group sequential analysis with binomial data.}
\description{The function \code{CV.G.Binomial} calculates the critical values for the group sequential MaxSPRT test with Binomial data, using a Wald-type
upper boundary, which is flat with respect to the likelihood ratio function, and with a pre-specified upper limit on the sample size (N).}
\usage{
CV.G.Binomial(N,alpha=0.05,M=1,z="n",p="n",GroupSizes)
      }
\arguments{
\item{N}{The upper limit on the sample size (length of surveillance) expressed in terms of the total number of events (cases plus controls). 
"N" must be a positive integer. To avoid very large computation times, we suggest not using values greater than 1000. Typically,
this is not a major restriction. For example, for "RR=1.1" and "alpha=0.01" and "z=1", the statistical power is approximately 1 for "N>500". There is no default value.}
\item{alpha}{The significance level. The "alpha" level must be in the range (0,0.5]. The default value is "alpha=0.05".}
\item{M}{The minimum number of cases needed before the null hypothesis can be rejected. It must be a positive integer, and the default value is 'M=1'}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case under the null hypothesis. There is no default value.}
\item{p}{The probability of having a case under the null hypothesis. There is no default value.}
\item{GroupSizes}{Vector with the number of events (cases+controls) between two consecutive looks (tests) at the data, i.e, the group sizes. The
length of this vector is equal to the maximum number of tests. The entries do not have to be the same, but they must sum up "N".
If the group sizes is an integer instead of a vector, then that integer is the group size for all looks at the data, and the number of looks is "N/GroupSizes".}
}
\details{
For group sequential analysis with binomial data, the function \code{CV.G.Binomial} calculates critical value for the decision of accepting/rejecting
the null hypothesis. This is done for a pre-defined significance level (alpha) and an upper limit on the sample size (N), as well as other parameter settings.
When the length of 'GroupSizes' is equal to one, N must be a multiple of GroupSizes.

The input z represents the number of controls matched to each case. For example, if there are 3 controls matched to each case, "z=3". 
In a self-control analysis, z is the ratio of the control interval to the risk interval.
For example, if the risk interval is 2 days long and the control interval is 7 days long, z=7/2.
In terms of p, the binomial probability under the null hypothesis, p=1/(1+z), or equivalently, z=1/p-1. The parameter z must be a positive number.

Alternatively, instead of z the user can specify p directly.
Note that only one of these inputs, z or p, has to be specified, but if both are entered the code will only work if z and p are such that p=1/(1+z).
Otherwise, an error message will appear to remind that such condition must be complied.

For several configurations of N, GroupSizes and M, there are no critical values for a true Type I error probability strictly equal to alpha. The is
because of the discrete nature of binomial data. For such cases, the function \code{CV.G.Binomial} returns the greatest critical value that will guarantee a
type I error probability that is smaller than 'alpha', so that the sequential analysis is conservative.  
}
\value{
\item{cv}{The critical value for a significance level equal to alpha. The largest conservative value is provided when
it is not possible to have an Type I error exactly equal to alpha.}
\item{Type_I_Error}{The exact Type I error probability given cv, and it is always less than or equal to alpha.}
}
\author{ Ivair Ramos Silva, Ned Lewis, Martin Kulldorff.
}
\references{
Fireman B, et al. (2013) Exact sequential analysis for binomial data with timevarying probabilities. Manuscript in Preparation.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

Kulldorff M, Silva IR. (2015). Continuous Post-market Sequential Safety Surveillance with Minimum Events to Signal. REVSTAT Statistical Journal, 15(3): 373--394.
}
\section{Acknowledgements}{
Development of the CV.G.Binomial function was funded by:\cr
-	Food and Drug Administration, Center for Drug Evaluation and Research, through the Mini-Sentinel Project (base version, documentation);\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (code revisions, increased computational speed,  improved documentation).
}

\keyword{Group MaxSPRT analysis}

\section{See also}{
\code{\link[Sequential]{SampleSize.Binomial}}: for calculating the minimum sample size given a target power in continuous sequential analysis with binomial data.\cr 
\code{\link[Sequential]{Analyze.Binomial}}: for performing sequential analysis with a user defined alpha spending function.\cr  
\code{\link[Sequential]{CV.Poisson}}: for calculating critical values for continuous sequential analysis with Poisson data. 
}

\examples{

# Example 1:
result<- CV.G.Binomial(N=25,alpha=0.05,M=1,z=7/2,GroupSizes=5)
# if you type:
result
# then you will get the following output:
# [[1]]
# [1] 1.9852

# [[2]]
# [1] 0.04775995

# Example 2:

result<- CV.G.Binomial(N=50,alpha=0.05,M=1,z=7/2,GroupSizes=c(10,10,15,15))
cv<- as.numeric(result[1])
# if you type:
cv
# then you will get the following output:
# [1] 1.99202
}


