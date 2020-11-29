\name{Sequential-package}
\alias{Sequential-package}
\alias{Sequential}
\docType{package}
\title{Analysis Support, Critical Values, Power, Time to Signal and Sample Size for Sequential Analysis with Poisson and Binomial Data.}
\description{
\code{Sequential} is designed for continuous and group sequential analysis, where statistical hypothesis testing is conducted repeatedly
on accumulating data that gradually increases the sample size. This is different from standard statistical analysis, where a single analysis is performed
using a fixed sample size. It is possible to analyze either Poisson type data or binomial 0/1 type data. For binomial data, it is possible to incorporate an
off-set term to account for variable matching ratios. For Poisson data, the critical value is based on a Wald-type upper boundary, which is flat on the scale
of the log-likelihood ratio, and on a predetermined maximum sample size. Alternatively, it is also possible to apply a user defined alpha spending function
in order to specify non-flat signaling thresholds. For group sequential analyses, there are functions for pre-specified group sizes and for the situation when
the group sizes are not known a priori. It is also possible to perform mixed continuous/group sequential analysis, where, for example, there is at first a big
batch of data that arrives in one group, followed by continuous sequential analysis. All results are exact, based on iterative numerical calculations, rather
than asymptotic theory or computer simulations. In the \code{Sequential} package, there are functions to calculate critical values, statistical power,
expected time to signal, and expected sample size at the end of the sequential analyses whether the null hypothesis was rejected or not.
For example, for any desired power, relative risk and alpha level, the package can calculate the required upper limit on the sample size (maximum length of surveillance),
the critical value needed, and the corresponding expected time to signal when the null hypothesis is rejected.}
\details{
\tabular{ll}{
Package: \tab Sequential\cr
Type: \tab Package\cr
Version: \tab 3.2\cr
Date: \tab 2020-11-26\cr
License: \tab GPL 2\cr
LazyLoad: \tab yes\cr
\tab \cr
Index: \tab \cr
\tab \cr
\code{\link[Sequential]{Performance.Threshold.Binomial}} \tab Power, Expected Time to Signal, Expected\cr
                                              \tab Sample Size and Alpha Spending for User Defined Threshold\cr
                                              \tab with Binomial Data.\cr
\code{\link[Sequential]{Performance.Threshold.CondPoisson}} \tab Power, Expected Time to Signal, Expected\cr
                                              \tab Sample Size and Alpha Spending for User Defined Threshold\cr
                                              \tab with Conditional Poisson Data.\cr
\code{\link[Sequential]{Performance.Threshold.Poisson}} \tab Power, Expected Time to Signal, Expected\cr
                                              \tab Sample Size and Alpha Spending for User Defined Threshold\cr
                                              \tab with Poisson Data.\cr
\code{\link[Sequential]{Analyze.Binomial}} \tab Function to Conduct Group Sequential Analyses for Binomial\cr
                                           \tab Data When the Goup Sizes are not Known a Priori.\cr
\code{\link[Sequential]{AnalyzeSetUp.Binomial}} \tab Function to Set Up the Input Parameters Before Using the\cr
                                                \tab \code{Analyze.Binomial} Function for the First Time.\cr
\code{\link[Sequential]{Analyze.wBinomial}} \tab Function for group sequential analyses of multiple\cr
                                           \tab weighted binomial endpoints.\cr
\code{\link[Sequential]{AnalyzeSetUp.wBinomial}} \tab Function to set up input parameters before using the\cr
                                                \tab \code{Analyze.wBinomial} function for the first time.\cr
\code{\link[Sequential]{Analyze.Poisson}} \tab Function to Conduct Group Sequential Analyses for Poisson\cr
                                           \tab Data When the Goup Sizes are not Known a Priori.\cr
\code{\link[Sequential]{AnalyzeSetUp.Poisson}} \tab Function to Set Up the Input Parameters Before Using the\cr
                                                \tab \code{Analyze.Poisson} Function for the First Time.\cr
\code{\link[Sequential]{Analyze.CondPoisson}} \tab Function to Conduct Group Sequential Analyses for Conditional\cr
                                           \tab Poisson Data When the Goup Sizes are not Known a Priori.\cr
\code{\link[Sequential]{AnalyzeSetUp.CondPoisson}} \tab Function to Set Up the Input Parameters Before Using the\cr
                                                \tab \code{Analyze.CondPoisson} Function for the First Time.\cr
\code{\link[Sequential]{CV.Binomial}} \tab Critical Values for Group and Continuous Sequential\cr
                                     \tab Analysis with Binomial Data.\cr
\code{\link[Sequential]{CV.Poisson}} \tab Critical Values for Group and Continuous Sequential\cr
                                     \tab Analysis with Poisson Data.\cr
\code{\link[Sequential]{CV.CondPoisson}} \tab Critical Values for Continuous Sequential CMaxSPRT for\cr
                                     \tab Limited Information from Historical Cohort Using\cr
                                     \tab Conditional Poisson Data.\cr
\code{\link[Sequential]{Optimal.Binomial}} \tab Optimal alpha spending for continuous and group sequential\cr
                                     \tab analysis with binomial data.\cr
\code{\link[Sequential]{Performance.Binomial}}  \tab Power, Expected Time to Signal and Expected Sample Size\cr
                                               \tab for Group and Continuous Sequential Analysis with\cr
                                               \tab Binomial Data.\cr 
\code{\link[Sequential]{Performance.Poisson}}  \tab Power, Expected Time to Signal and Expected Sample Size\cr
                                               \tab for Continuous Sequential Analysis with\cr 
                                                \tab Poisson Data.\cr
\code{\link[Sequential]{Performance.CondPoisson}}  \tab Power, Expected Time to Signal and Expected Sample\cr
                                               \tab Size for Continuous Sequential CMaxSPRT Using\cr 
                                                \tab Limited Information from Historical Cohort\cr
                                                 \tab with Conditional Poisson Data.\cr
\code{\link[Sequential]{SampleSize.Binomial}}  \tab Sample Size Calculation for Continuous Sequential Analysis with\cr
                                              \tab Binomial Data.\cr                              
\code{\link[Sequential]{SampleSize.Poisson}}  \tab Sample Size Calculation for Continuous Sequential Testing with\cr
                                              \tab Poisson Data.\cr
\code{\link[Sequential]{SampleSize.CondPoisson}}  \tab Sample Size Calculation for Continuous Sequential CMaxSPRT with\cr
                                              \tab Poisson Data.\cr
\code{\link[Sequential]{Performance.AlphaSpend.Binomial}} \tab Power, Expected Time to Signal, Expected\cr
                                              \tab Sample Size and Threshold for User Defined Alpha Spending \cr
                                              \tab with Binomial Data.\cr
\code{\link[Sequential]{Performance.AlphaSpend.CondPoisson}} \tab Power, Expected Time to Signal, Expected\cr
                                              \tab Sample Size and Threshold for User Defined\cr
                                              \tab Alpha Spending with Conditional Poisson Data.\cr
\code{\link[Sequential]{Performance.AlphaSpend.Poisson}} \tab Power, Expected Time to Signal, Expected\cr
                                            \tab Sample Size and Threshold for User Defined \cr
                                              \tab Alpha Spending with Poisson Data.\cr
}
}

\section{Overview}{
Most of the sequential analysis methods found in the literature are based on asymptotic results. In contrast,
this package contains functions for the exact calculation of critical values, statistical power, expected time
to signal when the null is rejected and the maximum sample size needed when the null is not rejected. This is done for
Poisson and binomial type data with a Wald-type upper boundary, which is flat with respect to the likelihood ratio function,
and a predetermined upper limit on the sample size. For a desired statistical power, it is also possible to calculate the latter.
The motivation for this package is post-market near real-time drug and vaccine safety surveillance, where the goal is to detect rare
but serious safety problems as early as possible, in many cases after only a hand full of adverse events. The package can also be used
in other application areas, such as clinical trials.

The basis for this package is the Maximized Sequential Probability Ratio Test (MaxSPRT) statistic (Kulldorff et al., 2011), which is a variant of
Wald's Sequential Probability Ratio Test (SPRT) (Wald, 1945,47). MaxSPRT uses a composite alternative hypothesis, and upper boundary to reject the
null hypothesis when there are more events than expected, no lower boundary, and an upper limit on the sample size at which time the sequential
analyses end without rejecting the null. MaxSPRT was developed for post-market vaccine safety surveillance as part of the Vaccine Safety Datalink
project run by the Centers for Disease Control and Prevention.

In this package, all critical values, alpha spending strategies, statistical power, expected time to
signal and required sample size to achieve a certain power, are obtained exactly to whatever decimal precision desired,
using iterative numerical calculations. None of the results are based on asymptotic theory or computer simulations.

\bold{Poisson Data}

To start, consider continuous sequential analysis for Poisson data. Let \eqn{C_t} be the random variable that counts the number of events up to time t.
Suppose that, under the null hypothesis, \eqn{C_t} has a Poisson distribution with mean \eqn{\mu_t}, where \eqn{\mu_t} is a known function reflecting the population at risk.
Under the alternative hypothesis, suppose that \eqn{C_t} has a Poisson distribution with mean RR\eqn{\mu_t}, where "RR" is the unknown increased relative risk due to the vaccine.
The MaxSPRT statistic defined in terms of the log likelihood ratio is given by: \deqn{LLR_t=(\mu_t-c_t)+c_t \log{c_t/\mu_t},} when \eqn{c_t} is at least \eqn{\mu_t}, and \eqn{LLR_t =0}, otherwise.
For continuous sequential analysis, the test statistic, \eqn{LLR_t}, is monitored at all times \eqn{t \in (0,T]}, where \eqn{T}= SampleSize. SampleSize is defined
a priori by the user in order to achieve the desired statistical power, which can be calculated using the \code{SampleSize.Poisson} function.
The sequential analysis ends, and \eqn{H_0} is rejected if, and when, \eqn{LLR_t \geq CV}, where CV is calculated using the \code{CV.Poisson} function.
If \eqn{\mu_t}= SampleSize, the sequential analysis ends without rejecting the null hypothesis. To calculate other important performance metrics, such as the expected time to signal when
the null hypothesis is rejected, use the \code{Performance.Poisson}, \code{Performance.Threshold.Poisson} and \code{Performance.AlphaSpend.Poisson} functions.
These functions also work for group sequential analysis when the group sizes are fixed a priori.
For with fixed and/or unpredictable group sizes, see the function \code{Analyze.Poisson}. 

If the first event occurs sufficiently early, the sequential analysis may end with the null hypothesis rejected after a single events. There is an option
to require a minimum number of observed events, \eqn{c_t}= M, before the null can be rejected.
Setting M in the range [3,6] is often a good choice (Kulldorff and Silva, 2012). If there is a delay until the sequential analysis starts,
but it continuous continuously thereafter, there is an option for that as well, requiring a minimum number \eqn{\mu_t}= D of expected events before the null can be rejected.

With continuous sequential analysis, investigators can repeatedly analyze the data as often as they want, ensuring that the overall probability of falsely
rejecting the null hypothesis at any time during the analysis is controlled at the desired nominal significance level (Wald, 1945, 1947). Continuous sequential
methods are suitable for real-time or near real-time monitoring. When data is only analyzed intermittently, group sequential methods are used instead
(Chin, 2012; Cook and DeMets, 2007; Xia, 2007; Friedman et al., 2010; Ghosh and Sen, 1991; Jennison and Turnbull, 2000; Mukhopadhyay and Silva, 2002; Whitehead, 1997).
The data is then analyzed at regular or irregular discrete time intervals after a certain amount of data is accessible. Group sequential statistical methods are commonly used in
clinical trials, where a trial may be stopped early due to either efficacy or unexpected adverse events (Jennison and Turnbull, 2000). 

The same test statistic, \eqn{LLR_t}, is used for group sequential analyses (Silva and Kulldorff, 2012). The times when \eqn{LLR_t} is evaluated can be defined in several ways,
using regular or irregular time intervals that are referenced by calendar period, sample size or some scale involving the distribution of the data. For Poisson data,
the group sequential analysis must be conducted with equal size groups, with a constant expected number of adverse events between looks at the accumulating data.
In another words, \eqn{LLR_t} is compared against CV whenever \eqn{\mu_t} is a multiple of SampleSize/Looks, where 'Looks' is the total number of looks at the data.

\bold{Binomial Data}

The MaxSPRT method can also be applied to binomial/Bernoulli data. Let \eqn{n} denote the total number of events that has been observed in a sequential monitoring
up to a certain moment in time. Suppose that these n events are categorized as cases and controls. For example, cases may be adverse events happening
to a person taking drug A, while controls may be the same adverse event happening to someone in a matched set of individuals taking drug B. As another example,
in a self-control sequential analysis, cases may be adverse events happening during the 1-28 days following vaccination, while controls are the same adverse events
happening 29-56 days after vaccination. 

Let \eqn{C_t} to denote the number of cases among the n events, and assume that \eqn{C_t} follows a binomial distribution with success probability equal to \eqn{p}, where \eqn{p = 1=(1 + z)},
and z is the matching ratio between the occurrence of a case and of a control under the null hypothesis. For example, if the probability of having a case (instead of a control)
is \eqn{p = 1=(1 + z) = 0.5}, then z=1 (1:1 matching ratio), or, \eqn{p = 0.25} for z=3 (1:3 matching ratio), etc.

The MaxSPRT statistic (Kulldorff et al., 2011) for a continuous binomial surveillance is:
\deqn{LR_n=\frac{(c_n/n)^{c_n}\left[(n-c_n)/n\right]^{n-c_n}}{\left[1/(1+z)\right]^{c_n}\left[z/(1+z)\right]^{n-c_n}},}
if \eqn{z c_n/(n-c_n)>1}, and \eqn{LR_n= 1} otherwise.

The monitoring is continued until either there is a signal rejecting the null hypothesis \eqn{(LR_n > CV)}
or until \eqn{n=N}, which indicates that the null is not to be rejected. To perform the calculations, use the
\code{CV.Binomial}, \code{SampleSize.Binomial}, \code{Performance.Binomial}, \code{Performance.Threshold.Binomial} and \code{Performance.AlphaSpend.Binomial} functions
for continuous and group sequential fashions. For group sequential analysis, the group sizes are fixed a priori.
For fixed and/or unpredictable group sizes, see the function \code{Analyze.Binomial}.

The main assumptions behind the method above are: (i) the monitoring is truly performed in a continuous fashion; (ii) the matching ratio (z) is constant for all of the n events,
and (iii) it uses a Wald type rejection boundary that is flat in terms of the likelihood function. Relaxing these assumptions, Fireman et al. (2013) developed exact sequential
analysis for group sequential data with varying matching ratios, and for any user specified alpha rejection plan.

\bold{Multiple Weighted Binomial Endpoints}

When there are multiple different outcomes, one can use weights reflecting practical interpretations, such as an outcome
severity in comparison to the others. For example, for two different outcomes, the first with weight w, and the second with weight 1.
This way, a single event of the first outcome is equivalent to w outcomes of the second type. This type of analysis is treated by Silva et al (2019).
To run the analysis with multiple weighted binary endpoints with fixed and/or unpredictable group sizes, see the function \code{Analyze.wBinomial}. 

\bold{Conditional Poisson data with limited information from historical cohort - CMaxSPRT}

In Poisson MaxSPRT, the expected mean \eqn{\mu_t} is assumed to be a known function reflecting the baseline adverse event risk
in the absence of the exposure of interest. In practice, it is estimated with historical data and the uncertainty associated with
the estimated counts may or may not have a non-negligible impact on the performance of the sequential analysis method.
Li and Kulldorff (Li and Kulldorff, 2010) showed in their simulation study that uncertainty in the estimated baseline means can be
ignored when the total number of events in the historical data is at least 5 times the specified upper limit \eqn{T}. Otherwise, it
is recommended to implement a statistical procedure that takes in account the variation from the historical data, i.e. a procedure that
conditionates the likelihood function of the historical Poisson data, here simply called "conditional Poisson data". For this, the
Conditional Maximized Sequential Probabilit Ratio Test (CMaxSPRT) to account for variation in both the historical and surveillance cohorts.

Let \eqn{c} and \eqn{V} denote the total number of events and the cumulative person-time in the historical data, let \eqn{P_k} denote the
cumulative person-time observed in the surveillance population when the \eqn{k}th event occurred. The CMaxSPRT statistic defined in terms
of the log likelihood ratio is given by \deqn{U_k=clog(\frac{c(1+P_k/V)}{c+k})+klog(\frac{k(1+P_k/V)}{P_k/V(c+k)}),} when \eqn{k/c>P_k/V},
and \eqn{U_k=0,} otherwise. In the original publication (Li and Kulldorff, 2010), the method was introduced as a continuous sequential analytic
approach with the upper limit defined in terms of the maximum number of observed events, i.e., \eqn{k \leq K}, and the critical value calculated
via a Monte Carlo approach. A large number of Monte Carlo simulations (e.g., 10 million) might be needed to calculate the critical values with
a reasonable precision. 

In Silva et al. (2016), the method was extended i) with another option of defining the surveillance length in terms of
the maximum cumulative person-time divided by the total cumulative person-time in the historical cohort, i.e., \eqn{P_k/V \le T}, ii) with an exact
calculation of the critical values for both surveillance length definitions, and iii) for group sequential analysis with data updated and analyzed
intermittently instead of continuously. The exact critical values are calculated using the interval havling method to solve for the root of a complex,
non-linear equation such that the overall Type I error rate is preserved at the nominal level. As \eqn{K} increases, the computing time for the exact
critical values increases exponentially. 

Silva et al. (2016) also proposed two approximation methods to calculate the critical values that require
substantially less computing time. One approch may overestimate the critical values and thus is referred to as the conservative approach as it may
yield lower-than-nominal Type I error rates; the other approach may underestmate the critical values and thus is referred to as the liberal approach
as it may yield higher-than-nominal Type I error rates. The recommendation is to use the exact approach when \eqn{K} is small (e.g., 10),
use the conservative approach when \eqn{K} is medium or large but \eqn{c} is small, and use the liberal approach when \eqn{c} is medium (e.g., 50)
or large. Exact calculations for selected tuning parameters show that the three approaches yield very similar results when \eqn{K} and \eqn{c} are reasonably large.

For calculating critical values for a Wald type rejection boundary, use the \code{CV.CondPoisson} function. For statistical power,
expected time to signal, expected time of surveillance, and maximum sample size requirements, use the
\code{Performance.CondPoisson}, \code{Performance.Threshold.CondPoisson}, \code{Performance.AlphaSpend.CondPoisson}, and \code{SampleSize.CondPoisson} functions.
For fixed and/or unpredictable group sizes, see the function \code{Analyze.CondPoisson}.

\bold{Alpha spending function for unpredictable group sizes}

The alpha spending function specifies the cumulative amount, \eqn{F_{\alpha}(t)}, of Type I error probability related to each of the possible values of \eqn{n}.
Thus, at the end of the monitoring the alpha spending corresponds to a value smaller than or equal to the overall amount of Type I error probability defined for
the overall nominal significance level, \eqn{\alpha}.

Denote the single probability of rejecting the null hypothesis at the \eqn{j}-th test by \eqn{\alpha_j}. Then, the alpha spending at test \eqn{i} is given by
\eqn{F_{\alpha}(t_i)=\sum_{j=1}^{i}\alpha_j \leq \alpha}.

There is a vast number of proposals for choosing the shape of the alpha spending function.
Jennison and Turnbull (2000) present a rich discussion about this topic. They dedicated a special attention to the alpha spending of the form:
\eqn{F_{\alpha}(t)=\alpha t^{\rho}}, where \eqn{\rho>0}, and \eqn{t} represents a fraction of the maximum length of surveillance. 

A new aproach for alpha spending selection is the optimal solution, introduced by Silva and Kulldorff (2018).
The optimal solution is based on exact calculations through linear programing, and it is operated through the function \code{Optimal.Binomial}.

To run continuous or group sequential analysis with an user defined alpha spending function, and/or, when the group sizes are not known a prior,
\code{Analyze.Binomial}, \code{Analyze.wBinomial} (for multiple weighted binomial endpoints), \code{Analyze.Poisson}, and \code{Analyze.CondPoisson} should be used for binomial and Poisson, and conditional Poisson data, respectively. 
These functions work differently than the other functions mentioned above.
Those other functions are designed to be used before the start of the sequential analysis, in order to determine what the maximum sample size
and critical value should be. Once the sequential analysis is under way, the test statistic is then calculated using a hand calculator or an
excel spread sheet, and compared with the critical value. The functions \code{Analyze.Binomial}, \code{Analyze.Poisson}, and \code{Analyze.CondPoisson} work very differently, in that they are run at each look at
the accumulating data, whenever a new group of data arrives, and it is meant to perform the test itself, i.e., there is no need to use hand calculators or
excel spread sheets or any other auxiliar code. The results and conclusions, including a descriptive table and illustrative graphics, are automatically
provided after running \code{Analyze.Binomial}, \code{Analyze.wBinomial}, \code{Analyze.Poisson}, or \code{Analyze.CondPoisson}.

Important: before using these functions, though, it is necessary to first run the
functions \code{AnalyzeSetup.Binomial}, \code{AnalyzeSetup.wBinomial}, \code{AnalyzeSetup.Poisson}, or \code{AnalyzeSetup.CondPoisson} once in order to set everything up for the sequential analysis.

}

\section{Comparison with Other R Packages for Sequential Analysis}{
The R  \code{Sequential} package is designed for sequential analysis where statistical hypothesis testing is performed using gradually accumulating data.
It is not designed for quality control problems, where a process is monitored over time to detect an emerging problem due to a sudden increase in the excess risk.
Although the methods for sequential analysis and quality control may seem similar, as they both analyze gradually accumulating data, they are actually very different
in both their purpose and design. Under the sequential hypothesis testing approach, the objective is to quickly determine if there is some intrinsic excess risk,
with the assumption that this risk does not change over time. For example, we may want to know if drug A is better than drug B, and there is no reason to believe
that the behavior of the drugs change over time. In the quality control setting, the objective is instead to detect a possible change in a stochastic process that
may occur in the future, and to detect that change as soon as possible after it occurs. For example, the heart of a hospital patient is beating as it should, but if
there is a sudden deterioration, the alarm should sound as soon as possible without generating a lot of false alarms. This package is only meant for sequential analysis
of the former type, and it should not be used for quality control type problems. For quality control type analyses, there are other R packages available,
such as \code{graphicsQC}, \code{IQCC}, \code{MetaQC}, \code{MSQC}, \code{qcc}, and \code{qcr}.

In a number of ways, the R Sequential package differs from other R packages for sequential analyses. Historically, most sequential analysis has been conducted
using asymptotic statistical theory, and that is also what is used in the \code{gsDesign}, \code{ldbounds}, \code{PwrGSD}, \code{seqDesign}, \code{seqmon}, and \code{sglr} R packages.
In contrast, the R Sequential package is based on exact results, using iterative numerical calculations, rather than using asymptotic theory or computer simulations.

With this package, it is only possible to analyze binomial/Bernoulli, Poisson, or conditional Poisson data. For other probability distributions, such as normal or exponential data,
other R packages should be consulted, such as \code{GroupSeq} or \code{SPRT}. Moreover, all functions in this package uses a one-sided upper bound to reject the null hypothesis,
while the analyses end without rejecting the null when an upper limit on the sample size is reached. For two sided sequential analysis, or other types of rejection
boundaries, other R packages must be used, such as e.g. \code{ldbounds} and \code{Binseqtest}. Finally, in this package, there are functions for both continuous
and group sequential analysis, and it is also possible to analyze situations where some of the data arrives continuously while other parts of the data arrives in groups.
Most other R packages are exclusively designed for group sequential analysis, but there are some that also do continuous sequential analysis, such as \code{Binseqtest} and \code{SPRT},
but \code{Binseqtest} is only for binomial data type, and \code{SPRT} is for simple alternative hypothis, while \code{Sequential} can be used for binomial and Poisson data and is meant to
composite alternative hypothesis. The present package offers the possibility to calculate the expected time to signal through the \code{Performance.Poisson}, \code{Performance.G.Poisson},
\code{Performance.Binomial}, \code{Performance.G.Binomial}, and Performance.CondPoisson functions, which is not offered by the other packages cited above.
Another important advantage of the \code{Sequential} package is the possibility of eliciting, through exact calculations, the minimum sample size needed to
accomplish with target statistical powers through the functions \code{SampleSize.Poisson}, \code{SampleSize.CondPoisson}, and \code{SampleSize.Binomial}.     
}


\author{ Ivair Ramos Silva, Martin Kulldorff.\cr
Maintainer: Ivair Ramos Silva <jamesivair@yahoo.com.br>
}

\references{
Chin R. (2012), Adaptive and Flexible Clinical Trials, Boca Raton, FL: Chapman and Hall/CRC.

Cook TD, DeMets DL. (2007), Introduction to Statistical Methods for Clinical Trials: Chapman and Hall/CRC Texts in Statistical Science. 

Fireman B, et al. (2013) Exact sequential analysis for binomial data with timevarying probabilities. Manuscript in Preparation.

Friedman LM, Furberg CD, DeMets D. (2010), Fundamentals of Clinical Trials, 4th ed.: Springer. 

Ghosh BK, Sen PK. (1991), Handbook of Sequential Analysis, New York: MARCEL DEKKER, Inc. 

Ghosh M, Mukhopadhyay N, Sen PK. (2011), Sequential Estimation: Wiley.

Jennison C, Turnbull B. (2000), Group Sequential Methods with Applications to Clinical
Trials, London: Chapman and Hall/CRC.

Kulldorff M, Davis RL, Kolczak M, Lewis E, Lieu T, Platt R. (2011). A Maximized Sequential Probability Ratio Test for Drug and Safety Surveillance. Sequential Analysis, 30: 58--78.

Kulldorff M, Silva IR. (2015). Continuous Post-market Sequential Safety Surveillance with Minimum Events to Signal. arxiv:1503.01978 [stat.ap]. Forthcoming paper of REVSTAT Statistical Journal.

Li L, Kulldorff M. (2010). A conditional maximized sequential probability ratio test for Pharmacovigilance. Statistics in Medicine, 29 (2), 284--295. 

Mukhopadhyay N, Silva BM. (2002), Sequential Methods and Their Applications, 1th ed.: Chapman and Hall/CRC. 

Silva IR, Kulldorff M. (2015), Continuous versus Group Sequential Analysis for Vaccine and Drug Safety Surveillance. Biometrics, 71 (3), 851--858.

Silva IR, Maro J, Kulldorff M. (2019). Exact Sequential Analysis Using R Sequential. Working Paper.

Silva, IR, and Kulldorff, M. (2018). Optimal Alpha Spending for Sequential Analysis With Binary Data. Working paper.

Silva IR, Gagne J, Najafzadeh M, Kulldorff M. (2019). Exact Sequential Analysis for Multiple Weighted Binomial Endpoints. Working paper.

Silva IR, Li L, Kulldorff M. (2019), Exact conditional maximized sequential probability ratio test adjusted for covariates. Sequential Analysis, 38(1), 115--133.

Silva IR. (2018a). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance with Binomial Data. Statistics in Medicine, 15;37(1), 107--118.

Silva IR. (2018b). Type I Error Probability Spending for Post-Market Drug and Vaccine Safety Surveillance With Poisson Data. Methodol Comput Appl Probab, 20(2), 739--750.

Silva IR., Lopes LM., Dias P., Yih WK. (2019). Alpha Spending for Historical Versus Surveillance Poisson Data With CMaxSPRT. Statistics in Medicine, 38(12), 2126--2138. 

Xia Qi. (2007), A Procedure for Group Sequential Comparative Poisson Trials. Journal of Biopharmaceutical Statistics, 17, 869--881. 

Wald A. (1945), Sequential Tests of Statistical Hypotheses, Annals of Mathematical Statistics, 16, 117--186.

Wald A. (1947), Sequential Analysis. New York: John Wiley and Sons.

Whitehead J. (1997), The Design and Analysis of Sequential Clinical Trials, 2th ed.: Wiley.
}

\section{Acknowledgements}{
Development of the R Sequential package has been funded and supported by:\cr
-	Food and Drug Administration, USA, through the Mini-Sentinel Project (v1.0,1.1,2.0).\cr
-	National Institute of General Medical Sciences, NIH, USA, through grant number R01GM108999 (v2.0 to 3.1).\cr
-	Federal University of Ouro Preto (UFOP), through contract under internal UFOP's resolution CEPE 4600 (v2.0 to 3.1).\cr
-	National Council of Scientific and Technological Development (CNPq), Brazil (v1.0).\cr
-	Bank for Development of the Minas Gerais State (BDMG), Brazil (v1.0).\cr
\cr
Feedback from users is greatly appreciated. Very valuable suggestions concerning the R Sequential package have been received from various individuals, including:\cr
-	Ron Berman, University of California Berkeley.\cr
-	Claudia Coronel-Moreno, Harvard Pilgrim Health Care Institute.\cr
-	Bruce Fireman, Kaiser Permanente Northern California.\cr
-	Josh Gagne, Harvard Medical School and Brigham and Women`s Hospital.\cr
-	Ned Lewis, Kaiser Permanente Northern California.\cr
-	Judith Maro, Harvard Medical School and Harvard Pilgrim Health Care Institute.\cr
-	Azadeh Shoaibi, Food and Drug Administration.\cr
-	Katherine Yih, Harvard Medical School and Harvard Pilgrim Health Care Institute.\cr
-     Jie Tang, Clinical biostatistics, Janssen R and D US, Johnson and Johnson LLC.\cr
-	Tuomo A. Nieminen, The National Institute for Health and Welfare (THL), Finland.\cr
-	Andreia Leite, Department of Infectious Disease Epidemiology, London School of Hygiene and Tropical Medicine.\cr 
-	Laura Hou, Harvard Pilgrim Health Care Institute, Boston, USA.\cr
-	Abdurrahman Abdurrob, Division of Pharmacoepidemiology and Pharmacoeconomics, Department of Medicine, Brigham and Women's Hospital and Harvard Medical School.\cr
-	Kirk Snyder, Information Management Services, Inc.\cr
}

\keyword{Sequential analysis}

\examples{

## Critical value for continuous sequential analyses for Poisson Data.
## Maximum sample size = 10, alpha = 0.05 and minimum number of events = 3:

cvt<- CV.Poisson(SampleSize=10,D=0,M=3,alpha=0.05)

## Statistical power and the expected time to signal for relative risk RR=2:

result<- Performance.Poisson(SampleSize=10,D=0,M=3,cv=cvt,RR=2)

# And if you type:
result

# Then you will see the following:
#          Power ESignalTime ESampleSize
#     [1,] 0.7329625    4.071636    5.654732

}

\section{Version History of the R Sequential Package}{
Version 1.1, February 2013\cr
Exact sequential analysis for Poisson data:\cr
-	Exact continuous sequential analysis.\cr
-	Exact group sequential analysis with pre-defined and constant groups sizes.\cr
-	Wald type rejection boundary.\cr
-	Statistical power, expected time to signal and sample size calculations.\cr
-	User guide.\cr
\cr
Version 1.2,  January 2014\cr
-	Improved code structure and efficiency.\cr
-	More extensive user guide. \cr
\cr
Version 2.0, June 2015\cr
Exact sequential analysis for binomial data:\cr
-	Continuous sequential analysis.\cr
-	Group sequential analysis with pre-defined group sizes.\cr
-	Group sequential analysis with unpredictable group sizes, not specified a priori.\cr
-	Fixed or variable binomial probabilities (matching ratios).\cr
-	User specified alpha spending function.\cr 
-	Statistical power, expected time to signal and sample size calculations.\cr
-	Updated user guide.\cr
\cr
Version 2.0.1, June 2015\cr
-       Correction of bugs in \code{CV.Poisson} function. \cr
-       Updated user guide.\cr
\cr
Version 2.0.2, Octuber 2015\cr
-       Improved user guide.\cr
\cr
Version 2.1, May 2016\cr
Exact sequential analysis for Poisson data:\cr
-	Group sequential analysis with unpredictable group sizes, not specified a priori.\cr
-	User specified alpha spending function.\cr
-	Mixed group-continuous sequential analysis.\cr
-	Statistical power, expected time to signal and sample size calculations for non-constant groups sizes.\cr
Other:\cr
-	Directory address parameter in AnalyzeSetUp functions.\cr
-	Probability parameter in binomial functions.\cr
-	Updated user guide.\cr
\cr
Version 2.1.1, June 2016\cr
-	Correction of bugs in Poisson functions.\cr
-	Updated user guide.\cr
\cr
Version 2.2, July 2016\cr
-	Critical Value, Performance, and SampleSize calculations for CMaxSPRT with Poisson data.\cr
-	Updated user guide.\cr
\cr
Version 2.2.1, September 2016\cr
-     Correction of bugs in \code{CV.Poisson} and \code{CV.G.Poisson} functions. \cr
-     Updated user guide.\cr
\cr
Version 2.3, Dec 2016\cr
-     Correction of bugs in the \code{SampleSize.Binomial} function.\cr
-	Improvement of SampleSize functions for considering vectors for the input parameters R and power.\cr
-	Inclusion of the new functions \code{AnalyzeSetUp.CondPoisson} and \code{Analyze.CondPoisson}.\cr
-     Updated user guide.\cr
\cr
Version 2.3.1, Feb 2017\cr
-     Correction of bugs in \code{Analyze.Binomial} and \code{AnalyzeSetUp.Poisson} functions.\cr
-	Adjustment on the relative risk estimation method for \code{Analyze.Binomial} function.\cr
-     Updated user guide.\cr
\cr
Version 2.3.2, Aug 2017\cr
-     Correction of bugs in \code{Analyze.Binomial}.\cr
-     Updated user guide.\cr
\cr
Version 3.0, Jan 2019\cr
-     Functions \code{Optimal.Binomial} and \code{Analyze.wBinomial}.\cr
-     Updated user guide.\cr
\cr
Version 3.0.1, Feb 2019\cr
-     Updated user guide.\cr
\cr
Version 3.1, Sept 2019\cr
-      New Functions \code{Performance.Threshold.Binomial}, \code{Performance.Threshold.CondPoisson},\cr
      \code{Performance.Threshold.Poisson}, \code{Performance.AlphaSpend.Binomial},\cr
      \code{Performance.AlphaSpend.CondPoisson}, \code{Performance.AlphaSpend.Poisson}.\cr
-     Implemented two-tailed testing for the functions: \code{CV.Poisson}, \code{Performance.Poisson},\cr
      \code{SampleSize.Poisson}, \code{Analyze.wBinomial}, \code{Optimal.Binomial},\cr
      \code{Performance.Threshold.Binomial}, \code{Performance.AlphaSpend.Binomial},\cr
      \code{Performance.Threshold.Poisson}.\cr
-     Updated user guide.\cr
-	Concatenation of functions for continuous and group fashions for critical values and performance.\cr 
\cr
Version 3.2, Nov 2020\cr
-     Function \code{Analyze.Binomial} improved for fixed-width and fixed-accuracy confidence intervals,\cr
-     Function \code{Performance.AlphaSpend.Binomial} improved for new inputs having lower and upper alpha spending in two tailed testing,\cr
-     Correction of bugs in \code{Analyze.Poisson} and \code{CV.Poisson} functions,\cr
-     Updated user guide.\cr
\cr
}
