% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/PseudoBoost.R
\name{PseudoBoost.data.frame}
\alias{PseudoBoost.data.frame}
\title{Perform stagewise pseudo-value regression for the risk 1 in a competing risk setting}
\usage{
\method{PseudoBoost}{data.frame}(data, xmat, times, stepno = 100,
  maxstepno = 100, nu = 0.1, cv = TRUE, multicore = FALSE, ...)
}
\arguments{
\item{data}{A data.frame containing observation times AND status (Coding: 0=Censored, 1=Interesting risk, >1=Other risks)}

\item{xmat}{A numeric matrix containing the covariate values for all patients.}

\item{times}{A numeric vector containing the evaluation times.}

\item{stepno}{A numeric value containing the number of boosting steps to be performed. If you use cross validation (cv=TRUE), this parameter will be ignored.}

\item{maxstepno}{A numeric value containing the maximal number of boosting steps considered during the cross validation (only used if cv=TRUE).}

\item{nu}{A numeric value between 0 and 1, the shrinkage parameter for the stagewise regression algorithm. Setting it to values such as nu=0.1 avoids overfitting in early steps.}

\item{cv}{A boolean value indicating if cross validation should be performed.}

\item{multicore}{A boolean value indication if more than one core should be used for the cross validation (for this the parallel package is needed).}
}
\value{
An object of type PseudoBoost containing the estimates and the performed boosting step number.
}
\description{
This function performs stagewise pseudo-value regression with the help of boosting for the CIF of the interesting risk.
The data frame must contain the observation times and statuses.
Thereby the statuses in the data frame should be natrual numbers and coded as follows: (0=Censored, 1=Interesting risk, >1=Other risks).
Additonally you have to put in the following arguments.
}

