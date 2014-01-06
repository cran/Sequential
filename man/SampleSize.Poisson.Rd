\name{SampleSize.Poisson}
\alias{SampleSize.Poisson}
\title{Sample Size Calculation for Continuous Sequential Analysis with Poisson Data.}
\description{The function \code{SampleSize.poisson} obtains the maximum sample size (length of surveillance) needed to guarantee a desired statistical power
when applying the MaxSPRT test with a Wald type upper boundary, which is flat with respect to the log-likelihood ratio.}
\usage{
SampleSize.Poisson(alpha=0.05,power=0.9,M=1,D=0,RR=2,precision=0.000001)
      }
\arguments{
\item{alpha}{The significance level. The default value is 'alpha=0.05'. To avoid numerical complications on the calculation 
of the critical value for values of alpha close to zero or one, 'alpha' must be in the range (0,0.5].}
\item{power}{The target overall power to detect an increased risk of the relative risk (RR). The default value is 'power=0.90'.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. The default value is 'M=1', 
which means that even a single event can reject the null hypothesis if it occurs sufficiently early. A good rule of thumb is to set 'M=4'
(Kulldorff and Silva, 2012).}
\item{D}{The expected number of events under the null hypothesis before the first look at the data. The default value is "D=0", which is
also the best choice. This means that there is no delay in the start of the sequential analyses. If D is very large, the maximum 
sample size will be set equal to D if a non-sequential analysis provides the desired power.}
\item{RR}{A target relative risk to be detected with the requested statistical power. The default value is 'RR=2'}
\item{precision}{The tolerance for how much larger is the allowed power with respect to the target one. The default value is 'precision=0.000001'.}
}
\details{
The function \code{SampleSize.Poisson} calculates the upper limit on the sample size (length of surveillance) that needs to be set for the
continuous Poisson based MaxSPRT in order to provide the user specified statistical power for the user specified relative risk (RR). 
The solution is exact using iterative numerical calculations. For details, see Kulldorff et al. (2011).}
\value{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) needed by the continuous Poisson based MaxSPRT 
in order to achive ethe desired statistical "power" for the user specified relative risk "RR".}
}
\author{ Ivair Ramos Silva, Martin Kulldorff.

Maintainer: Ivair Ramos Silva <jamesivair@yahoo.com.br>

}
\references{

Ivair R. Silva, Martin Kulldorff. (2012). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Manuscript, to be submitted.

Martin Kulldorff, Ivair R. Silva. (2012). Continuous Sequential Analysis with a Delayed Start. Munuscript, under review.

Martin Kulldorff, Robert L. Davis, Margarette Kolczak, Edwin Lewis, Tracy Lieu, and Richard Platt. (2011). A Maximized Sequential Probability Ratio Test 
for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

}

\keyword{Continuous MaxSPRT analysis}

\section{See also}{
\code{\link[Sequential]{Performance.Poisson}}, \code{\link[Sequential]{Performance.G.Poisson}}, 

\code{\link[Sequential]{CV.G.Poisson}}, \code{\link[Sequential]{CV.Poisson}},  \code{\link[stats]{poisson.test}}.
}


\examples{

#  Upper limit on the sample size required to obtain a power of 90%, for a relati-
#  ve risk of 3, no delay for starting the surveillance (D=0), and the null hypothe-
#  sis can be rejected with only one event (M=1).

# SampleSize.Poisson(alpha=0.05,power=0.8,M=1,D=0,RR=3)

}


