\name{CV.Binomial}
\alias{CV.Binomial}
\title{Calculates exact critical values for group and continuous sequential analysis with binomial data.}
\description{The function \code{CV.Binomial} obtains critical values for the group continuous sequential MaxSPRT test with binomial data, using a Wald-type upper boundary,
which is flat with respect to the likelihood ratio function, and an pre-specified upper limit on the sample size.}
\usage{
CV.Binomial(N,alpha=0.05,M=1,z="n",p="n",GroupSizes=1,Tailed="upper")
      }
\arguments{
\item{N}{The upper limit on the sample size (length of surveillance) expressed in terms of the total number of events (cases plus controls). 
"N" must be a positive integer. To avoid very large computation times, we suggest not using values greater than 1000. Typically,
this is not a major restriction. For example, for "RR=1.1", "alpha=0.01" and "z=1", the statistical power is approximately 1 for "N>500". There is no default value.}
\item{alpha}{The significance level. The "alpha" level must be in the range (0,0.5]. The default value is "alpha=0.05".}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. "M" must be a positive integer, and the default value is "M=1".}
\item{z}{For a matched case-control analysis, z is the number of controls matched to each case under the null hypothesis. There is no default value.}
\item{p}{The probability of having a case under the null hypothesis. There is no default value.}
\item{GroupSizes}{Vector with the number of events (cases+controls) between two consecutive looks (tests) at the data, i.e, the group sizes. The
length of this vector is equal to the maximum number of tests. The entries do not have to be the same, but they must sum up "N".
If the group sizes is an integer instead of a vector, then that integer is the group size for all looks at the data, and the number of looks is "N/GroupSizes". The default is GroupSizes=1 for continuous sequential analysis.}
\item{Tailed}{Tailed="upper" (default) for H0:RR<=1, and Tailed="lower" for H0:RR>=1 or Tailed="two" for H0:RR=1.}
}
\details{
For the continuous and group binomial MaxSPRT, \code{CV.Binomial} calculates the upper boundary used to 
determine if the null hypothesis is to be rejected at each analysis. This is done for pre-specified values of the statistical significance level (alpha) and 
an upper limit on the sample size equal to N. 

The input z represents the number of controls matched to each case. For example, if there are 3 controls matched to each case, "z=3". 
In a self-control analysis, z is the ratio of the control interval to the risk interval.
For example, if the risk interval is 2 days long and the control interval is 7 days long, z=7/2.
In terms of p, the binomial probability under the null hypothesis, p=1/(1+z), or equivalently, z=1/p-1. The parameter z must be a positive number.

Alternatively, instead of z the user can specify p directly.
Note that only one of these inputs, z or p, has to be specified, but if both are entered the code will only work if z and p are such that p=1/(1+z).
Otherwise, an error message will appear to remind that such condition must be complied.

For details about the algorithm used to calculate the critical value, see the paper by Kulldorff et al. (2011).

For some configurations of N and alpha and GroupSizes, there is no critical value that gives a Type I error probability that is exactly equal to the
requested "alpha". This is because of the discrete nature of binomial data. In such situations, \code{CV.Binomial} returns the greatest critical value that guarantees a Type I error probability smaller than "alpha". Thus 
the critical value for the binomial sequential analysis is conservative in this sense. 
}
\value{
\item{cv}{The critical value for a significance level equal to alpha. The largest conservative value is provided when
it is not possible to have an Type I error exactly equal to alpha.}
\item{Type_I_Error}{The exact Type I error probability given cv. Always less than or equal to alpha.}
}
\author{ Ivair Ramos Silva, Ned Lewis, Ron Berman, Martin Kulldorff.
}
\references{
Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.
}


\section{Acknowledgements}{
Development of the CV.Binomial function was funded by:\cr
-	Food and Drug Administration, Center for Drug Evaluation and Research, through the Mini-Sentinel Project; base version, documentation;\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999; code revisions, increased computational speed, improved documentation.\cr
\cr
We thank Ron Berman, University of California, Berkeley, for a key suggestion to speed up the calculations, and Bruce Fireman for helpful discussions.
}

\keyword{Continuous and group MaxSPRT analysis}

\section{See also}{
\code{\link[Sequential]{Analyze.Binomial}}: for performing sequential analysis with group, continuous or unpredictable sequential fashion.
}


\examples{

# Example 1:
## Critical value for continuous binomial sequential analysis with
#  a maximum sample size of 20 events, requiring at 
#  least 3 events to reject the null, and with a significance level of 0.05:

CV.Binomial(N=20,alpha=0.05,M=3,z=1.1)

# Example 2:
## Critical value for five-group sequential analysis with
#  a maximum sample size of 25 events, requiring at 
#  least 1 event to reject the null, and with a significance level of 0.05:
result<- CV.Binomial(N=25,alpha=0.05,M=1,z=7/2,GroupSizes=5)
# if you type:
result
# then you will get the following output:
# [[1]]
# [1] 1.9852

# [[2]]
# [1] 0.04775995

# Example 3:
## Critical value for four-group sequential analysis with
#  a maximum sample size of 50 events, requiring at 
#  least 1 event to reject the null, and with a significance level of 0.05:
result<- CV.Binomial(N=50,alpha=0.05,M=1,z=7/2,GroupSizes=c(10,10,15,15))
cv<- as.numeric(result[1])
# if you type:
cv
# then you will get the following output:
# [1] 1.99202

}


