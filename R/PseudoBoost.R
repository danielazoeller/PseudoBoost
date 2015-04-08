#' Perform stagewise pseudo-value regression for the interesting risk in a competing risk setting
#' 
#' This function performs stagewise pseudo-value regression with the help of boosting for the CIF of the interesting risk.
#' Calls the functions to perform Boosting for pseudo-values for the CIF of the first risk in a competing risk setting depending on the type of the object.
#' The object can either be a data frame containing the observation times and statuses or a numeric matrix containing the pseudo-values for the CIF for the interesting risk obtained with the help of the package "prodlim" (jackknife(prodlim(Hist(observation-times-vector,status-vector) ~ 1), times=evaluation-times-vector, cause=interesting-cause-number)).
#' Thereby the statuses in the data frame should be natrual numbers and coded as follows: (0=Censored, 1=Interesting risk, >1=Other risks). 
#' Additonally you have to put in the following arguments.
#' @param object A data.frame containing observation times AND statuses or a numeric matrix containing the pseudo-values for the CIF for the risk 1
#' @param xmat A numeric matrix containing the covariate values for all patients.
#' @param times A numeric vector containing the evaluation times.
#' @param stepno A numeric value containing the number of boosting steps to be performed. If you use cross validation (cv=TRUE), this parameter will be ignored.
#' @param maxstepno A numeric value containing the maximal number of boosting steps considered during the cross validation (only used if cv=TRUE).
#' @param nu A numeric value between 0 and 1, the shrinkage parameter for the stagewise regression algorithm. Setting it to values such as nu=0.1 avoids overfitting in early steps.
#' @param cv A boolean value indicating if cross validation should be performed.
#' @param multicore A boolean value indication if more than one core should be used for the cross validation (for this the parallel package is needed).
#' @return An object of type PseudoBoost containing the estimates and the performed boosting step number.
#' @import parallel
#' @import prodlim
#' @export 
PseudoBoost <- function(object,...){
  UseMethod("PseudoBoost",object)
}

#' Perform stagewise pseudo-value regression for the risk 1 in a competing risk setting
#' 
#' This function performs stagewise pseudo-value regression with the help of boosting for the CIF of the interesting risk.
#' The data frame must contain the observation times and statuses.
#' Thereby the statuses in the data frame should be natrual numbers and coded as follows: (0=Censored, 1=Interesting risk, >1=Other risks). 
#' Additonally you have to put in the following arguments.
#' @param data A data.frame containing observation times AND status (Coding: 0=Censored, 1=Interesting risk, >1=Other risks)
#' @param xmat A numeric matrix containing the covariate values for all patients.
#' @param times A numeric vector containing the evaluation times.
#' @param stepno A numeric value containing the number of boosting steps to be performed. If you use cross validation (cv=TRUE), this parameter will be ignored.
#' @param maxstepno A numeric value containing the maximal number of boosting steps considered during the cross validation (only used if cv=TRUE).
#' @param nu A numeric value between 0 and 1, the shrinkage parameter for the stagewise regression algorithm. Setting it to values such as nu=0.1 avoids overfitting in early steps.
#' @param cv A boolean value indicating if cross validation should be performed.
#' @param multicore A boolean value indication if more than one core should be used for the cross validation (for this the parallel package is needed).
#' @return An object of type PseudoBoost containing the estimates and the performed boosting step number.
#' @import parallel
#' @import prodlim
#' @export 
PseudoBoost.data.frame <- function(data,xmat,times,stepno=100,maxstepno=100,nu=0.1,cv=TRUE,multicore=FALSE,...){
  require(prodlim)
  if(ncol(data)==2){
    data <- as.data.frame(cbind(rep(0,nrow(data)),data))
  } else if(ncol(data)!=3){
    stop("The data object must be a data.frame containing either the start times, stop times and statuses or the observed times and statuses!")
  } 
  #ymat <- PseudoCIF(data,times)
  ymat <- jackknife(prodlim(Hist(data[,2],data[,3]) ~ 1), times=times, cause=1)
  if (cv) cv.res <- cv.PseudoBoost(ymat,xmat,maxstepno=maxstepno,nu=nu,multicore=multicore,...)
  stepno <- cv.res$optimal.step
  res <- PseudoBoost(ymat,xmat,stepno=stepno,nu=nu,...)
  res$times <- times
  
  return(res) 
}

#' Perform stagewise pseudo-value regression for the first risk in a competing risk setting
#' 
#' This function performs stagewise pseudo-value regression with the help of boosting for the CIF of the interesting risk.
#' The ymat a numeric matrix containing the pseudo-values for the CIF for the interesting risk obtained with the help of the package "prodlim" (jackknife(prodlim(Hist(observation-times-vector,status-vector) ~ 1), times=evaluation-times-vector, cause=interesting-risk-number)).
#' Thereby the statuses in the data frame should be natrual numbers and coded as follows: (0=Censored, 1=Interesting risk, >1=Other risks). 
#' Additonally you have to put in the following arguments.
#' @param ymat A numeric matrix containing the pseudo-values for the CIF for the risk 1
#' @param xmat A numeric matrix containing the covariate values for all patients.
#' @param times A numeric vector containing the evaluation times.
#' @param stepno A numeric value containing the number of boosting steps to be performed. If you use cross validation (cv=TRUE), this parameter will be ignored.
#' @param maxstepno A numeric value containing the maximal number of boosting steps considered during the cross validation (only used if cv=TRUE).
#' @param nu A numeric value between 0 and 1, the shrinkage parameter for the stagewise regression algorithm. Setting it to values such as nu=0.1 avoids overfitting in early steps.
#' @param cv A boolean value indicating if cross validation should be performed.
#' @param multicore A boolean value indication if more than one core should be used for the cross validation (for this the parallel package is needed).
#' @return An object of type PseudoBoost containing the estimates and the performed boosting step number.
#' @import parallel
#' @import prodlim
#' @export 
PseudoBoost.default <- function(ymat,xmat,stepno=10,nu=0.1,trace=FALSE) {
  n <- nrow(xmat)
  toffset <- matrix(0,n,ncol(ymat))
  actual.offset <- rep(0,n)
  beta.inter <- rep(0,ncol(ymat))
  x.inter <- t(rep(1,n))
  beta.mat <- matrix(0,ncol(ymat),ncol(xmat))
  
  mean.y <- colMeans(ymat)
  penalty <- mean.y*(1-mean.y)*n*(1/nu-1)
  
  beta.list <- list()
  
  #   log link
  #actual.score <- sum((exp(actual.offset)*actual.x*(actual.y - exp(actual.offset))))
  #actual.fisher <- -sum(exp(actual.offset)*actual.x^2*(actual.y - exp(actual.offset)) + exp(actual.offset)*actual.x*-exp(actual.offset)*actual.x)
  #   logit link
  get.score <- function(x,y,offset) sum(exp(-offset)*x/(1+exp(-offset))^2 * (y - 1/(1+exp(-offset))))
  get.fisher <- function(x,y,offset) (exp(-offset)*x^2*(y*exp(-offset)^2 - y - 2*exp(-offset) + 1))/
    (exp(-offset)^4+4*exp(-offset)^3+6*exp(-offset)^2+4*exp(-offset)+1)
  
  #   boosting steps
  for (actual.step in 1:stepno) {    
    #   update the intercept terms            
    beta.update <- rep(1,ncol(ymat))
    beta.dirty <- rep(TRUE,ncol(ymat))
    
    iter.max <- 2000
    actual.iter <- 1
    while (any(beta.dirty) && actual.iter < iter.max) {    
      exp.toffset <- exp(toffset[,beta.dirty,drop=FALSE])
      
      #   log
      
      actual.score <- colSums(exp.toffset * (ymat[,beta.dirty,drop=FALSE] - exp.toffset))
      actual.fisher <- colSums(exp.toffset^2)
      #actual.fisher <- -colSums(ymat[,beta.dirty,drop=FALSE]*exp.toffset-2*exp.toffset^2)
      
      
      #   logit
      #div.exp.m.toffset <- exp.m.toffset/(1+exp.m.toffset)^2
      #exp.m.toffset <- exp(-toffset[,beta.dirty,drop=FALSE])
      #y.minus <- ymat[,beta.dirty,drop=FALSE] - 1/(1+exp.m.toffset)
      ##actual.score <- get.score.vec(drop(x.inter),ymat,toffset)
      #actual.score <- colSums(div.exp.m.toffset * y.minus)
      ##actual.fisher <- get.fisher.vec(drop(x.inter),ymat,toffset)
      
      ##actual.fisher <- -colSums((-1 + exp.m.toffset)*div.exp.m.toffset/(1+exp.m.toffset) * y.minus - div.exp.m.toffset^2)
      
      #actual.fisher <- -colSums((exp.m.toffset*(ymat[,beta.dirty,drop=FALSE]*(exp.m.toffset^2 - 1)-2*exp.m.toffset+1))/
      #                          (exp.m.toffset^4+4*exp.m.toffset^3+6*exp.m.toffset^2+4*exp.m.toffset+1))
      
      
      beta.update[] <- 0
      beta.update[beta.dirty] <- ifelse(is.finite(actual.fisher) & (actual.fisher > 0),actual.score/actual.fisher,0)
      
      #beta.update[beta.dirty] <- actual.score/actual.fisher
      
      beta.inter <- beta.inter + beta.update
      toffset <- t(t(toffset) + beta.update)
      
      beta.dirty <- abs(beta.update) > 0.001
      actual.iter <- actual.iter + 1
    }
    if (actual.iter == iter.max) cat("While-loop stopped \n")
    
    if (actual.step == 1) beta.list[[length(beta.list)+1]] <- cbind(beta.inter,beta.mat)
    
    exp.toffset <- exp(toffset)
    
    #   log
    scoremat <- t(crossprod(xmat,exp.toffset * (ymat - exp.toffset)))
    fishermat <- t(crossprod(xmat^2,exp.toffset^2))
    #fishermat <- -t(crossprod(scale(xmat,center=FALSE)^2,ymat*exp.toffset-2*exp.toffset^2))
    
    #   logit
    #exp.m.toffset <- exp(-toffset)
    #div.exp.m.toffset <- exp.m.toffset/(1+exp.m.toffset)^2
    #y.minus <- ymat - 1/(1+exp.m.toffset)    
    #scoremat <- t(crossprod(xmat,div.exp.m.toffset * y.minus))
    #fishermat <- -t(crossprod(xmat^2,exp.m.toffset*(ymat*(exp.m.toffset^2 - 1) - 2*exp.m.toffset + 1)/
    #                                 (exp.m.toffset^4+4*exp.m.toffset^3+6*exp.m.toffset^2+4*exp.m.toffset+1)))
    
    #fishermat <- -t(crossprod(xmat^2,(-1 + exp.m.toffset)*div.exp.m.toffset/
    #                 (1+exp.m.toffset) * y.minus) - 
    #               crossprod(xmat^2,div.exp.m.toffset^2))
    
    #fishermat <- ifelse(fishermat > 0,fishermat,matrix(apply(fishermat,2,function(arg) min(arg[arg > 0])),nrow=nrow(fishermat),ncol=ncol(fishermat),byrow=TRUE))
    
    #fishermat < t(ifelse(t(fishermat) > 0,t(fishermat),apply(fishermat,2,function(arg) min(arg[arg > 0]))))
    
    #cat("******* score:")
    #print(scoremat)
    #cat("******* fisher:")
    #print(fishermat)
    
    #cat("******* T:")
    #print(scoremat^2/(fishermat))
    
    #   score statistic for each time point
    scoret <- scoremat^2/(fishermat)# + penalty)
    #   combine by Fisher's method; degrees of freedom not important, as only ranking needed
    comb.t <- -2*colSums(log(ifelse(is.finite(scoret),1-pchisq(scoret,df=1),1)))
    
    #   update coefficient of the covariate that has the largest score statistic
    selected.index <- which.max(comb.t)
    
    if (trace) cat(selected.index," ",sep="")
    
    #cat("******* selected:")
    #print(selected.index)
    
    #   regularized update; only fraction \nu of the unregularized update 
    #beta.update <- scoremat[,selected.index]/(fishermat[,selected.index]+penalty)
    
    #cat("******* update:")
    #print(scoremat/fishermat)
    
    #if(any(is.na(fishermat[,selected.index])))  print(fishermat[,selected.index])
    #if(any(is.na(scoremat[,selected.index])))  print(scoremat[,selected.index])
    
    beta.update <- nu*ifelse(fishermat[,selected.index] > 0,scoremat[,selected.index]/fishermat[,selected.index],0)
    #beta.update <- ifelse(fishermat[,selected.index] > 0,nu*scoremat[,selected.index]/fishermat[,selected.index],0)
    #beta.update <- scoremat[,selected.index]/(fishermat[,selected.index]+penalty)
    #print(beta.update)
    #beta.update <- scoremat[,selected.index]/(ifelse(is.finite(fishermat[,selected.index]),fishermat[,selected.index],0)+penalty)
    #print(beta.update)
    
    #EINGEFUEGT! DANIELA ZOELLER 03.03.
    # beta.update[is.na(beta.update)]=0
    # beta.update[is.nan(beta.update)]=0
    beta.update = ifelse(is.finite(beta.update),beta.update,0)
    
    beta.mat[,selected.index] <- beta.mat[,selected.index] + beta.update
    if(any(is.na(beta.mat)))  print(beta.mat)
    if(any(is.nan(beta.mat)))  print(beta.mat)
    #print(beta.mat[,5:9])
    for (i in 1:ncol(ymat)) toffset[,i] <- toffset[,i] + xmat[,selected.index]*beta.update[i]
    
    beta.list[[length(beta.list)+1]] <- cbind(beta.inter,beta.mat)
  }
  
  if (trace) cat("\n")
  
  res <- list(coefficients=beta.list,stepno=stepno)
  class(res) <- "PseudoBoost"
  return(res)
}

#' Calculates predicted linear predictor or response (CIF) based on PseudoBoost-Object
#' 
#' @param object A PseudoBoost-object
#' @param newdata A n.new * p matrix with new covariate values.
#' @param subset An optional vecotr specifying a subset of observations to be used for evaluation
#' @param at.step A scalar of boosting step at which prediction is wanted.
#' @param times A vector with T time points where prediction is wanted.
#' @param type Type of prediction to be returned: "lp" gives the linear predictor, "response" the resulting CIF.
#' @export 
predict.PseudoBoost <- function(object,newdata,subset=NULL,at.step=NULL,times=NULL,type = c("lp","response")) {
  type <- match.arg(type)
  
  if (is.null(at.step)) {
    at.step <- object$stepno
  }
  
  if (!is.null(times) && is.null(object$times)) stop("'times' missing, and no default available")
  
  if (is.null(subset)) {
    subset.index <- 1:nrow(newdata)
  } else {
    subset.index <- (1:nrow(newdata))[subset]
  }
  
  #   fit obtained from formula interface
  if (!is.null(object$terms)) {    
    newdata <- model.matrix(object$terms, data = model.frame(object$formula, data = newdata))
  } else {
    newdata <- cbind(rep(1,nrow(newdata)),newdata)
  }
  
  
  linear.predictor <- lapply(object$coefficients[at.step],function(arg) tcrossprod(newdata,arg))
  
  
  #   interpolation
  if (!is.null(times)) {
    eval.index <- unlist(lapply(ifelse(times < min(object$times),min(object$times),times),
                                function(arg) rev(which(object$times <= arg))[1]))    
    linear.predictor <- lapply(linear.predictor,function(arg) arg[,eval.index,drop=FALSE])
  }
  
  if (type == "response") {
    if (length(linear.predictor) == 1) {
      #return(1/(1+exp(-linear.predictor[[1]])))
      return(exp(linear.predictor[[1]]))
    } else {
      return(lapply(linear.predictor,function(arg) 1/(1+exp(-arg))))
    }
  }
  
  if (length(linear.predictor) == 1) {
    return(linear.predictor[[1]])
  }
  
  linear.predictor
}

#' Determines the optimal number of boosting steps by a K-fold cross-validation.
#' 
#' @param ymat Matrix of Pseudo Values
#' @param xmat n* p matrix of covariates
#' @param subset A vector specifying a subset of observations to be used in the fitting process.
#' @param maxstepno Maximum number of boosting steps to evaluate, i.e., the returned "optimal" number of boosting steps will be in the range [0,maxstepno]
#' @param K Number of folds to be used for cross-validation. If K is larger or equal to the number of non-zero elements in status, leave-one-out-cross-validation is performed.
#' @param multicore Indicates whether computations in the cross-validation folds should be performed in parallel, using package parallel.
#' @param folds If not NULL, this has to be a list of length K, each element being a vector of indices of fold elements. Useful for employing the same folds for repeated runs.
#' @param trace Logical vlue indicating whether progress in estimation should be indicated by printing the number of the cross-validation fold and the index of the covariate updated.
#' @export 
cv.PseudoBoost <- function (ymat,xmat,subset=1:nrow(ymat),maxstepno=100,K = 10,multicore=FALSE,folds=NULL,trace=FALSE,...) {
  if (!is.null(folds) && length(folds) != K) stop("'folds' has to be of length 'K'")
  
  subset <- (1:length(time))[subset]
  
  if (is.null(folds)) folds <- split(sample(1:length(subset)), rep(1:K,length=length(subset)))
  
  eval.fold <- function(actual.fold, ...) {
    if (trace) cat("cv fold ", actual.fold, ": ", sep = "")
    
    cv.fit <- PseudoBoost.default(ymat[-folds[[actual.fold]],],xmat[-folds[[actual.fold]],],
                                  stepno=maxstepno,trace=trace,...)
    
    unlist(lapply(predict(cv.fit,newdata=xmat[folds[[actual.fold]],],at.step=0:maxstepno,type="response"),
                  function(arg) mean((ymat[folds[[actual.fold]],] - arg)^2)))
  }
  
  eval.success <- FALSE
  
  if (multicore) {
    if (!require(parallel)) {
      warning("package 'parallel' not found, i.e., parallelization cannot be performed using this package")
    }
    else {
      if (multicore > 1) {
        criterion <- matrix(unlist(mclapply(1:length(folds),eval.fold, mc.preschedule=FALSE,mc.cores = multicore,...)),
                            nrow=length(folds),byrow=TRUE)
      }
      else {
        criterion <- matrix(unlist(mclapply(1:length(folds),eval.fold, mc.preschedule=FALSE,...)),
                            nrow=length(folds),byrow=TRUE)
      }
      eval.success <- TRUE
    }
  }
  if (!eval.success) {
    criterion <- matrix(unlist(lapply(1:length(folds), eval.fold,...)),nrow=length(folds),byrow=TRUE)
  }
  
  mean.criterion <- apply(criterion,2, mean)
  
  list(mean.ipec=mean.criterion,optimal.step = which.min(mean.criterion)-1,folds=folds)
}

#' Calculates predicted response (CIF) based on PseudoBoost-Object
#' 
#' @param object A PseudoBoost-object
#' @param newdata A n.new * p matrix with new covariate values.
#' @param subset An optional vecotr specifying a subset of observations to be used for evaluation
#' @param at.step A scalar of boosting step at which prediction is wanted.
#' @param times A vector with T time points where prediction is wanted.
#' @export 
predictProb.PseudoBoost <- function (object,newdata,subset=NULL,at.step=NULL,times=NULL) {
  
  if (is.null(at.step)) {
    at.step <- object$stepno
  }
  
  if (!is.null(times) && is.null(object$times)) stop("'times' missing, and no default available")
  
  if (is.null(subset)) {
    subset.index <- 1:nrow(newdata)
  } else {
    subset.index <- (1:nrow(newdata))[subset]
  }
  
  predict(object,newdata=newdata,at.step=at.step,type="response")
}