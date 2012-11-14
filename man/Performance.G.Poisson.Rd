\name{Performance.G.Poisson}
\alias{Performance.G.Poisson}
\title{Power, Expected Signal Time and Expected Length of Serveillance from the Poisson Group MaxSPRT.}
\description{The \code{Performance.G.Poisson} function calculates several performance 
metrics for the continuous Poisson based MaxSPRT. The user specifies the upper limit on the sample 
size (SampleSize), the upper boundary critical value (cv), the maximum number of tests to be perfomed (Looks), the minimum number of events 
required before rejecting the null hypothesis (M) and a relative risk (RR). 
The metrics calculated are the statistical power, the expected time to signal when the null hypothesis is 
rejected and the expected sample size at the end of the analysis whether the null hypothesis was rejected or not. If the user wants to calculate
the statistical significance level alpha, that can be done by setting RR=1, in which case the power output value is the alpha level.}
\usage{
Performance.G.Poisson(SampleSize,cv,Looks=5,M=1,RR=1.5)
}
\arguments{
\item{SampleSize}{the upper limit on the length of surveillance, expressed in terms of the number of events under the null hypothesis. "SampleSize" 
must be greater than 0 and, aiming to avoid unecessary computational effort, it is suggested to do not use values greater than 1000.}
\item{cv}{the critical value for a group sequential MaxSPRT test (see \code{\link[Sequential]{CV.G.Poisson}}).}
\item{Looks}{The maximum number of tests to be performed in the period "SampleSize". It is understood that each test will be performed after each "SampleSize/Looks"
units with respect to the expected number of cases under the null.}
\item{M}{the minimum number of cases for which a signal is allowed to occur(default is 1).}
\item{RR}{the relative risk (>=1) to calculate the power and the expected time to have a signal (default is 1).}
}
\details{
The function \code{Performance.G.Poisson} calculates the power, the expected time to have a signal and the expected
time of surveillance by using the group sequential MaxSPRT for poisson data, given a flat boundary, an upper limit on the number "SampleSize"
of expected events under the null, and a delayed start in the vigilance given by "M".
It is important to point out that large values of "SampleSize", as "SampleSize>1000", for
example, leads to a high computational effort. 

The times when each test will be performed is defined by 
using regular time intervals that are referenced by constant expected number of adverse events between 
looks at the accumulated data. In another words, the MaxSPRT statistic is compared against cv whenever the time is a multiple of \eqn{SampleSize/Looks}. 
}
\value{
\item{Power}{the statistical power for a given flat upper boundary (CV).}
\item{Expected time to signal}{the expected time to have a signal.}
\item{Expected length of surveillance}{the expected length of surveillance.}
}
\author{Ivair Ramos Silva and Martin Kulldorff
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

\section{See also}{
\code{\link[Sequential]{CV.Poisson}}, \code{\link[Sequential]{SampleSize.Poisson}}, \code{\link[Sequential]{CV.G.Poisson}}, 

\code{\link[Sequential]{Performance.Poisson}}, \code{\link[Sequential]{TableCV.PoissonD}}, 

\code{\link[Sequential]{TableCV.PoissonM}}, \code{\link{poisson.test}}
}
\examples{
## calculates the critical value for a significance level of 0.05.

cvt<- CV.G.Poisson(SampleSize=30,alpha=0.05,Looks=5,M=1)

## calculates power, expected signal time and expected surveillance time for a 
#  relative risk equal to 2

#Performance.G.Poisson(SampleSize=30,cvt,Looks=5,M=1,RR=1.5)
}

