\name{CV.Poisson}
\alias{CV.Poisson}
\title{Critical Values for Continuous Sequential Analysis with Poisson Data.}
\description{The function \code{CV.Poisson} obtains critical values for the continuous sequential MaxSPRT test with Poisson data, using a Wald 
type upper boundary, which is flat with respect to the likelihood ratio function, and an pre-specified upper limit on the sample size.}
\usage{
CV.Poisson(SampleSize,D=0,M=1,alpha=0.05)
      }
\arguments{
\item{SampleSize}{The upper limit on the sample size (length of surveillance) expressed in terms of the epected number of events under the null hypothesis. 
"SampleSize" must be greater than 0. To avoid very large computation times, we suggest not using values greater than 1000. Typically,
this is not a major restriction. For example, for "RR=1.1" and "alpha=0.01", the statistical power is approximately 1 for a maximum 
sample size greater than 500. There is no default value.}
\item{D}{The expected number of events under the null hypothesis before the first look at the data. The default is "D=0", which is
also the best choice. This means that there is no delay in the start of the sequential analyses. It is required that 'D<=SampleSize'. }
\item{M}{The minimum number of events needed before the null hypothesis can be rejected. The default value is 'M=1', 
which means that even a single event can reject the null hypothesis if it occurs sufficiently early. A good rule of thumb is to set 'M=4'
(Kulldorff and Silva, 2012).}
\item{alpha}{The significance level, or the type 1 error probability, which is the probability of rejecting the null hypothesis when it is true. The alpha level  must be in the range (0,0.5]. The default value is 'alpha=0.05'.}
}
\details{
For the continuous Poisson based MaxSPRT, \code{CV.Poisson} calculates the upper boundary used to 
determine if the null hypothesis should be rejected. This is done for pre-specified values of the statistical significance level (alpha) and 
an upper limit on the sample size, determining the maximum length of surveillance. For details about the
algorithm used to calculate the critical value see Kulldorff et al. (2011). 

For some configurations of 'SampleSize', 'D' and 'alpha', there is no critical value that gives a significance level that is exactly equal to the
requested alpha. In such situations, \code{CV.Poisson} returns the greatest critical value that will guarantee an alpha level less than the alpha specified, 
so that sequential analysis is conservative.

For large values of 'SampleSize', such as 200 or more, 
the computational requirements can be high. For selected values, the function will use one of two look-up tables, \code{TableCV.PoissonD} 
and \code{TableCV.PoissonM}, to speed things up. 
}
\value{
\item{cv}{The critical value for a significance level equal to alpha. The largest conservative value is provided when
it is not possible to have an Type I error exactly equal to 'alpha'.}
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

\keyword{Continuous MaxSPRT analysis}

\section{See also}{
\code{\link[Sequential]{Performance.Poisson}}, \code{\link[Sequential]{Performance.G.Poisson}}, \code{\link[Sequential]{CV.G.Poisson}},

\code{\link[Sequential]{SampleSize.Poisson}}, \code{\link[stats]{poisson.test}}.
}


\examples{

## Critical value for continuous sequential analysis with
#  a maximum sample size of ten expected cases under the null hypothesis, requiring at 
#  least 3 events to reject the null, and with a significance level of 0.05:

CV.Poisson(SampleSize=10,D=0,M=3,alpha=0.05)

## In the next example, no critical value exist that will give the desired 0.05 alpha level
#  exactly. Instead, a criticical value is provided that makes the alpha as large as
#  possible without exceeding 0.05.

 CV.Poisson(SampleSize=3,D=3,M=1,alpha=0.05)

}


