\name{TableCV.PoissonM}
\alias{TableCV.PoissonM}
\docType{data}
\title{
Critical Values Table, Continuous Poisson with Minimum Case Requirement.}
\description{
This table contains exact critical values for the continuous Poisson based MaxSPRT, with a significance level of 0.05. 
Critical values are provided for different options for both the upper limit on the sampel size 'SampleSize' and
on the minimum number of events 'M' needed before rejecting the null hypothesis. 
}
\usage{data(TableCV.PoissonM)}
\format{
  A matrix containing critical values for 30 different values for "SampleSize" and 7 different values for "M".
  }
\details{
The function \code{\link[Sequential]{CV.Poisson}} uses an analitycal procedure to obtain 
exact critical values to be used as flat upper boundaries to perform the MaxSPRT continuous sequential test. When the maximum length of surveillance (SampleSize),
specified by the user, is relatively large, like "SampleSize= 1000", for example, the algorithm can take a long time to calculate the 
critical value. Then, we provide the table \code{TableCV.PoissonD} of critical values to save computing time for a significance level equal to 0.05,
minimum number of expected events to reject the null hypothesis (D) equal to zero and a large set of combinations for "SampleSize" and 
minimum number of observed events to reject the null (M). Each line of \code{TableCV.PoissonD} represents a different "SampleSize" value in the set 
1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10, 12, 15, 20, 25, 30, 40 ,50, 60, 80, 100, 120, 150, 200, 250, 300, 400, 500, 600, 800, 1000, 
and each column represents a "D" value in the set 1, 2, 3, 4, 6, 8, 10. Observe that, for "M=1", the procedure is just the standart continuous sequential test.
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

\section{See also}{
\code{\link[Sequential]{CV.Poisson}}, \code{\link[Sequential]{CV.G.Poisson}}, \code{\link[Sequential]{Performance.G.Poisson}}, 

\code{\link[Sequential]{Performance.Poisson}}, \code{\link[Sequential]{SampleSize.Poisson}},

\code{\link[stats]{poisson.test}}.
}

\examples{

## Taking the critical value to perform the MaxSPRT test for alpha=0.05,
#  minimum number of expected events to reject the null equal to 0,
#  SampleSize=50 and M=10.

data(TableCV.PoissonM)
TableCV.PoissonM[rownames(TableCV.PoissonM)==50,7]
}

