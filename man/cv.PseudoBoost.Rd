% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/PseudoBoost.R
\name{cv.PseudoBoost}
\alias{cv.PseudoBoost}
\title{Determines the optimal number of boosting steps by a K-fold cross-validation.}
\usage{
cv.PseudoBoost(ymat, xmat, subset = 1:nrow(ymat), maxstepno = 100, K = 10,
  multicore = FALSE, folds = NULL, trace = FALSE, ...)
}
\arguments{
\item{ymat}{Matrix of Pseudo Values}

\item{xmat}{n* p matrix of covariates}

\item{subset}{A vector specifying a subset of observations to be used in the fitting process.}

\item{maxstepno}{Maximum number of boosting steps to evaluate, i.e., the returned "optimal" number of boosting steps will be in the range [0,maxstepno]}

\item{K}{Number of folds to be used for cross-validation. If K is larger or equal to the number of non-zero elements in status, leave-one-out-cross-validation is performed.}

\item{multicore}{Indicates whether computations in the cross-validation folds should be performed in parallel, using package parallel.}

\item{folds}{If not NULL, this has to be a list of length K, each element being a vector of indices of fold elements. Useful for employing the same folds for repeated runs.}

\item{trace}{Logical vlue indicating whether progress in estimation should be indicated by printing the number of the cross-validation fold and the index of the covariate updated.}
}
\description{
Determines the optimal number of boosting steps by a K-fold cross-validation.
}

