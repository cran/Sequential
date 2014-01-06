\name{CV.G.Poisson}
\alias{CV.G.Poisson}
\title{Critical Values for Group Sequential Analysis with Poisson Data.}
\description{The function \code{CV.G.Poisson} calculates the critical values for the group sequential MaxSPRT test with Poisson data, using a Wald 
type upper boundary, which is flat with respect to the likelihood ratio function, and with a pre-specified upper limit on the sample size.}
\usage{
CV.G.Poisson(SampleSize,alpha=0.05,Looks=1,M=1)
      }
\arguments{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) expressed in terms of the expected number of events under the null hypothesis. 
"SampleSize" must be greater than 0. To avoid very large computation times, we suggest not using values greater than 1000. Typically,
this is not a major restriction. For example, for "RR=1.1" and "alpha=0.01", the statistical power is approximately 1 for a maximum 
sample size greater than 500. There is no default value.}
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. The default value is 'M=1'. 
If there are frequent looks at the data, so that the group sizes are extremely small, a value of 'M=1' means that even a single event can reject the null hypothesis if it occurs sufficiently early.}
\item{alpha}{The significance level, or the type 1 error probability, which is the probability of rejecting the null hypothesis when it is true. The alpha level  must be in the range (0,0.5]. The default value is 'alpha=0.05'.}
\item{Looks}{The maximum number of tests to be performed, which is, the number of equally spaced looks at the data if the null hypothesis is never rejected. The value must be an integer greater or equal to 1. The default value is 'Looks=1'.}
}
\details{
For group sequential analysis with Poisson data, \code{CV.G.Poisson} calculates the critical value that constitutes the upper boundary used to 
determine if the null hypothesis should be rejected. This is done for pre-specified values of the statistical significance level (alpha) and 
an upper limit on the sample size, determining the maximum length of surveillance, as well as other parameter settings. Details about the
algorithm used to calculate the critical value is provided by Silva and Kulldorff (2012). 

For several configurations of 'SampleSize', 'Looks' and 'M' there is no critical value that gives a probability of Type I error that is exactly equal to "alpha".
In such cases, the function \code{CV.G.Poisson} returns the largest critical value that will guarantee a
type I error probability that is smaller than 'alpha', so that the sequential analysis is conservative.      

For large values of the maximum 'SampleSize', such as 200 or more, 
the computational requirements can be high. 
}
\value{
\item{cv}{The critical value for a significance level equal to alpha. The largest conservative value is provided when
it is not possible to have a type I error that is exactly equal to 'alpha'.}
}
\author{ Ivair Ramos Silva, Martin Kulldorff.

Maintainer: Ivair Ramos Silva <jamesivair@yahoo.com.br>

}
\references{

Ivair R. Silva, Martin Kulldorff. (2012). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Manuscript, to be submitted.

Jennison, V. and Turnbull, B. (1999), Group Sequential Methods with Applications to Clinical
Trials, \emph{no. ISBN 0-8493-0316-8, London: Chapman and Hall/CRC}.

Martin Kulldorff, Ivair R. Silva. (2012). Continuous Sequential Analysis with a Delayed Start. Munuscript, under review.

Martin Kulldorff, Robert L. Davis, Margarette Kolczak, Edwin Lewis, Tracy Lieu, and Richard Platt. (2011). A Maximized Sequential Probability Ratio Test 
for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.


}

\keyword{Group MaxSPRT analysis}

\section{See also}{
\code{\link[Sequential]{Performance.Poisson}}, \code{\link[Sequential]{Performance.G.Poisson}}, \code{\link[Sequential]{CV.Poisson}},

\code{\link[Sequential]{SampleSize.Poisson}},

\code{\link[stats]{poisson.test}}.
}


\examples{

#  Calculates the critical value for five equally spaced group sequential looks, at 6, 12,
#  18, 24 and 30 expected events under the null hypothesis, and for a statistical signifi-
#  cance level of 0.05. 

CV.G.Poisson(SampleSize=30,alpha=0.05,Looks=5,M=1)

}


