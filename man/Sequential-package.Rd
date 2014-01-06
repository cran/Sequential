\name{Sequential-package}
\alias{Sequential-package}
\alias{Sequential}
\docType{package}
\title{Exact Critical Values, Power, Time to Signal and Sample Size for Sequential Analysis with Poisson Data.}
\description{
This package calculates the critical values (CV) to perform continuous or group sequential surveillance by using
the maximized sequential probability ratio test statistic (MaxSPRT) for Poisson data. The critical value provided is based on a Wald type upper boundary, 
which is flat on the scale of the log-likelihood ratio, and on a predetermined maximum sample size. The package \code{Sequential} also calculates the 
statistical power and the expected time to signal when the null hypothesis is rejected. For any desired power and relative risk,
the package can also calculate the required upper limit on the sample size, and the corresponding expected time to signal. All results are exact, 
using iterative numerical calculations, rather than using asymptotic theory or computer simulations.}
\details{
\tabular{ll}{
Package: \tab Sequential\cr
Type: \tab Package\cr
Version: \tab 1.01\cr
Date: \tab 2014-01-06\cr
License: \tab GPL 2\cr
LazyLoad: \tab yes\cr
\tab \cr
Index: \tab \cr
\tab \cr
\code{\link[Sequential]{CV.Poisson}} \tab Critical Values for Continuous Sequential Analysis with\cr
                                     \tab Poisson Data.\cr                                
\code{\link[Sequential]{SampleSize.Poisson}}  \tab Sample Size Calculation for Continuous Sequential Analysis with.\cr
                                              \tab Poisson Data.\cr                                   
\code{\link[Sequential]{Performance.Poisson}}  \tab Power, Expected Signal Time and Sample Size for Continuous Sequen-\cr
                                               \tab tial Analysis with Poisson Data.\cr
\code{\link[Sequential]{CV.G.Poisson}} \tab Critical Values for Group Sequential Analysis with Poisson Data.\cr
\code{\link[Sequential]{Performance.G.Poisson}}  \tab Power, Expected Signal Time and Sample Size for Group Sequential\cr
                                                 \tab  Analysis with Poisson Data.\cr
}
}

\section{Overview}{
Most sequential analysis is based on asymptotic results. This package contains functions for the exact calculation of critical values, statistical power, 
expected time to signal when the null is rejected and the maximum sample size needed when the null is not rejected. This is done for Poisson type date 
with a Wald type upper boundary, which is flat with respect to the likelihood ratio function, and a predetermined upper limit on the sample size. For a 
desired statistical power, it is also possible to calculate the latter. The motivation for this package is post-market near real-time drug and vaccine 
safety surveillance, where the goal is to detect rare but serious safety problems as early as possible, in many cases after only a hand full of 
adverse events. There are also other application areas. 

The basis for this package is the MaxSPRT statistic (Kulldorff et al., 2011), which is a variant of Wald's Sequential Probability Ratio Test (SPRT)
(Wald, 1945,47). MaxSPRT uses a composite alternative hypothesis, and upper bounbdary to reject the null hypothesis when there are more events than expected, 
no lower boundary, and an upper limit on the sample size at which time the sequential analyses end without rejecting the null. 
MaxSPRT was developed for post-market vaccine safety surveillance as part of the Vaccine Safety Datalink project run by the Centers for Disease 
Control and Prevention.

Let \eqn{C_t} be the random variable that counts the number of events up to time t. Suppose that, under
the null hypothesis, \eqn{C_t} has a Poisson distribution with mean \eqn{\mu_t}, where \eqn{\mu_t} is a known function
reflecting the population at risk. Under the alternative hypothesis, suppose that \eqn{C_t} has a Poisson
distribution with mean \eqn{RR C_t}, where "RR" is the unknown increased relative risk due to the vaccine. The MaxSPRT statistic defined
in terms of the log likelihood ratio is given by: 

 \deqn{LLR_t=(\mu_t-c_t)+c_t \log{c_t/\mu_t},}{LLR_t= (\mu_t-c_t)+ c_t log(c_t/\mu_t),}

when \eqn{c_t} is at least \eqn{\mu_t}, and \eqn{LLR_t =0}, otherwise. For continuous sequential analysis, the test statistic \eqn{LLR_t} is monitored at all times t > 0.
The sequential analyses end, and \eqn{H_0} is rejected, 
if and when \eqn{LLR_t \geq CV}. If \eqn{\mu_t=SampleSize=SS}, the sequential analysis ends without rejecting the null hypothesis. "SS" is defined a
priori by the user in order to achieve the desired statistical power. 

If the first event occur sufficiently early, the sequential analyses may end with the null hypothesis rejected after a single events.
There is an optional to require a minimum number of observed events, \eqn{c_t=M}, before the null can be rejected. Setting M in the range [3,6] 
is often a good choice (Kulldorff and Silva, 2012). If there is a delay until the sequential analysis starts, but it continuous continously 
thereafter, there is an option for that as well, requireing a minimum number \eqn{\mu_t=D} of expected events before the null can be rejected. 

With continuous sequential analysis, investigators can repeatedly analyze the data as often as they want, ensuring
that the overall probability of falsely rejecting the null hypothesis at any time during the analysis
is controlled at the desired nominal significance level (Wald, 1945, 1947). Continuous sequential methods are suitable for real-time or near real-time 
monitoring. When data is only analyzed intermittently, group sequential methods are used instead (Jennison and Turnbull, 1999).
The data is then analyzed at regular or irregular discrete time intervals after a certain amount of data is accessible.
Group sequential statistical methods are commonly used in clinical trials, where a trial may be
stopped early due to either efficacy or unexpected adverse events (Jennison and Turnbull, 1999). 

The same test statistic, \eqn{LLR_t}, is used for group sequential analyses (Silva and Kulldorff, 2012
). 
The times when \eqn{LLR_t} is evaluated can be defined in several ways, 
using regular or irregular time intervals that are referenced by calendar period, sample size or some scale involving the distribution of the data. 
In this first version of the package, the group sequential analysis must be conducted with a constant expected number of adverse events between 
looks at the accumulated data. In another words, \eqn{LLR_t} is compared against CV whenever \eqn{\mu_t} is a multiple of \eqn{SS/L}, where L is 
the total number of looks at the data. 

In this package, all critical values, statistical power, epected time to signal and required sample size to achive a 
certain power, are obtained exactly to whatever decimal precision desired, using iterative numerical calculations. None of the results are based on 
asymptotic theory or computer simulations.
}

\author{ Ivair Ramos Silva, Martin Kulldorff.

Maintainer: Ivair Ramos Silva <jamesivair@yahoo.com.br>

}
\section{Acknowledgements}{This package was developed during a one year post-doctoral visit by Dr. Silva at the Department of Population Medicine, 
Harvard Medical School and Harvard Pilgrim Health Care Institute, Boston, MA 02215, USA. The funds supporting this work came from the 
United States Food and Drug Administration, Center for Biologics Evaluation and Research, through
the Mini-Sentinel Post-Rapid Immunization Safety Monitoring program; from the National Council of Scientific and 
Technological Development (CNPq), Minas Gerais, Brazil; and from the Bank for the Development of the Minas Gerais State (BDMG), Minas Gerais, Brazil. 
We are grateful to Claudia Coronel-Moreno for editorial support.}


\references{

Abraham Wald. (1945). Sequential Tests of Statistical Hypotheses, \emph{Annals of Mathematical Statistics},
 \bold{16}, 117--186.

Ivair R. Silva, Martin Kulldorff. (2012). Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Manuscript, to be submitted.

Jennison, V. and Turnbull, B. (1999), Group Sequential Methods with Applications to Clinical
Trials, \emph{no. ISBN 0-8493-0316-8, London: Chapman and Hall/CRC}.

Martin Kulldorff, Ivair R. Silva. (2012). Continuous Sequential Analysis with a Delayed Start. Munuscript, under review.

Martin Kulldorff, Robert L. Davis, Margarette Kolczak, Edwin Lewis, Tracy Lieu, and Richard Platt. (2011). A Maximized Sequential Probability Ratio Test 
for Drug and Safety Surveillance. \emph{Sequential Analysis}, \bold{30}: 58--78.

}

\keyword{Sequential analysis}

\section{See also}{
\code{\link[Sequential]{CV.Poisson}}, \code{\link[Sequential]{SampleSize.Poisson}}, \code{\link[Sequential]{CV.G.Poisson}}, \code{\link[Sequential]{Performance.Poisson}},

\code{\link[Sequential]{Performance.G.Poisson}}, \code{\link[stats]{poisson.test}}, \code{\link[gsDesign]{gsDesign}}.
}


\examples{

## Critical value for a maximum sample size = 10, alpha = 0.05 and minimum number of events = 3:

cvt<- CV.Poisson(SampleSize=10,D=0,M=3,alpha=0.05)

## Statistical power and the expected time for RR=2:

#result<- Performance.Poisson(SampleSize=10,D=0,M=3,cv=cvt,RR=2)


}



