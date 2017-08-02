\name{CV.Poisson}
\alias{CV.Poisson}
\title{Critical Values for Continuous Sequential Analysis with Poisson Data.}
\description{The function \code{CV.Poisson} obtains critical values for the continuous sequential MaxSPRT test with Poisson data, using a Wald 
type upper boundary, which is flat with respect to the likelihood ratio function, and with a pre-specified upper limit on the sample size.}
\usage{
CV.Poisson(SampleSize,D=0,M=1,alpha=0.05)
      }
\arguments{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) expressed in terms of the expected number of events under the null hypothesis.
The SampleSize must be greater than 0. To avoid very large computation times, we suggest not using values greater than 1000. Typically,
this is not a major restriction. For example, for RR=1.1 and alpha=0.01, the statistical power is approximately 1 for a maximum 
sample size greater than 500. There is no default value.}
\item{D}{The expected number of events under the null hypothesis before the first look at the data. The default is D=0, which is
also the best choice. This means that there is no delay in the start of the sequential analyses. It is required that D<=SampleSize. }
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. The default value is M=1, 
which means that even a single event can reject the null hypothesis if it occurs sufficiently early. A good rule of thumb is to set M=4
(Kulldorff and Silva, 2015).}
\item{alpha}{The significance level, or the type 1 error probability, which is the probability of rejecting the null hypothesis when it is true. The alpha level  must be in the range (0,0.5]. The default value is alpha=0.05.}
}
\details{
For the continuous sequential analysis with Poisson data, using the maximized sequential probability ratio test (MaxSPRT), \code{CV.Poisson} calculates the upper boundary used to 
determine if the null hypothesis should be rejected. This is done for pre-specified values on the statistical significance level (alpha) and
the upper limit on the sample size, determining the maximum length of surveillance. The
algorithm used to calculate the critical value is described by Kulldorff et al. (2011). 

For some configurations of SampleSize, D and alpha, there is no critical value that gives a significance level that is exactly equal to the
requested alpha. In such situations, \code{CV.Poisson} returns the greatest critical value that will guarantee an alpha level less than the alpha specified, 
so that sequential analysis is conservative.

For large values of SampleSize, such as 200 or more, 
the computational requirements can be high. To speed things up, the function will sometimes use one of two look-up tables that contain
pre-calculated critical values for a pre-selected set of parameter values (\code{TableCV.PoissonD} and \code{TableCV.PoissonM}). 
}
\value{
\item{cv}{The critical value for a significance level equal to alpha. The largest conservative value is provided when
it is not possible to have a Type I error exactly equal to alpha.}
}
\author{ Ivair Ramos Silva, Martin Kulldorff.
}
\references{
Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

Kulldorff M, Silva IR. (2015). Continuous Post-market Sequential Safety Surveillance with Minimum Events to Signal. REVSTAT Statistical Journal, 15(3): 373--394.

}

\section{Acknowledgements}{
Development of the CV.Poisson function was funded by:\cr
-	Food and Drug Administration, Center for Biologics Evaluation and Research, through the Mini-Sentinel Post-Rapid Immunization Safety Monitoring (PRISM) program (v1.0).\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil (v1.0).\cr
-	Bank for Development of the Minas Gerais State (BDMG), Brazil (v1.0).\cr
-     National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0.1, 2.0.2).
}

\keyword{Continuous MaxSPRT analysis}

\section{See also}{
\code{\link[Sequential]{SampleSize.Poisson}}: for calculating the sample size needed for Continuous Sequential Analysis with Poisson Data.\cr
\code{\link[Sequential]{Performance.Poisson}}: for calculating the statistical power, expected time to signal and expected time of analysis for continuous sequential analysis with Poisson data.\cr
\code{\link[Sequential]{CV.G.Poisson}}: for calculating critical values for group sequential analysis with Poisson data.\cr
\code{\link[Sequential]{CV.Binomial}}: for calculating critical values in continuous sequential analysis with binomial data.
}


\examples{

## Calculates the critical value for continuous sequential analysis with
## a maximum sample size of ten expected cases under the null hypothesis,
## requiring at least 3 events to reject the null, and with a significance
## level of 0.05:

CV.Poisson(SampleSize=10,D=0,M=3,alpha=0.05)

## In this example, no critical value exist that will give the desired 0.05
## alpha level exactly. Instead, the function produces the critical value
## that makes the alpha as large as possible without exceeding 0.05.

CV.Poisson(SampleSize=3,D=1,M=1,alpha=0.05)

}


