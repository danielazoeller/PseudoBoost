#' Perform smoothed stagewise pseudo-value regression for the first risk in a competing risk setting and calculate confidence bands based on pointwise confidence intervals.
#' 
#' Creates RepMean modified datasets with the help of pairwise switches of the observation times with on observation time smaller than switch.t. For these dataset the function PseudoBoost will be called seperately and the mean value is taken for each effect estimate at each considered time point.
#' @param data A data.frame containing observation times AND statuses
#' @param xmat A numeric matrix containing the covariate values for all patients.
#' @param times A numeric vector containing the evaluation times.
#' @param stepno A numeric value containing the number of boosting steps to be performed. If you use cross validation (cv=TRUE), this parameter will be ignored.
#' @param maxstepno A numeric value containing the maximal number of boosting steps considered during the cross validation (only used if cv=TRUE).
#' @param nu A numeric value between 0 and 1, the shrinkage parameter for the stagewise regression algorithm. Setting it to values such as nu=0.1 avoids overfitting in early steps.
#' @param cv A boolean value indicating if cross validation should be performed.
#' @param multicore A boolean value indication if more than one core should be used for the cross validation (for this the parallel package is needed).
#' @param RepMean A numeric value indicating the number of modified datasets used the estimate.
#' @param switch.to A numeric value indicating the number of the first risk up to which the observations can be switched. It sets the value switch.t. If switch.t is specified, switch.to will be ignored.
#' @param switch.t A numeric value indicating the time point up to which the observations can be switched. If this value is not NULL, switch.to will be ignored.
#' @param seed.start A numeric value indicating the first seed value. 
#' @param trace A boolean value indicating if additonal information should be printed out during the process.
#' @param switches A numeric value indicating the number of switches to be made.
#' @param trace A boolean value indicating if additonal information should be printed to the console (=TRUE)
#' @param part A numeric value between 0 and 1 indicating the fraction which should be used in the resampling approach (=0.632=
#' @param RepCb A numeric value indicating the number of Resampling data sets (=100)
#' @param alpha A numeric value between 0 and 1 indicating the type 1 error (=0.05)
#' @param cv_est A boolean value indicating if cross validation should be performed for the estimates seperately (=TRUE)
#' @return An object of type ciMeanPseudoBoost containing the estimates and the performed boosting step number.
#' @export 

cbMeanPseudoBoost <- function(object,...){
  UseMethod("cbMeanPseudoBoost",object)
}

#' Perform smoothed stagewise pseudo-value regression for the first risk in a competing risk setting and calculate confidence bands based on pointwise confidence intervals.
#' 
#' Creates RepMean modified datasets with the help of pairwise switches of the observation times with on observation time smaller than switch.t. For these dataset the function PseudoBoost will be called seperately and the mean value is taken for each effect estimate at each considered time point.
#' @param data A data.frame containing observation times AND statuses
#' @param xmat A numeric matrix containing the covariate values for all patients.
#' @param times A numeric vector containing the evaluation times.
#' @param stepno A numeric value containing the number of boosting steps to be performed. If you use cross validation (cv=TRUE), this parameter will be ignored.
#' @param maxstepno A numeric value containing the maximal number of boosting steps considered during the cross validation (only used if cv=TRUE).
#' @param nu A numeric value between 0 and 1, the shrinkage parameter for the stagewise regression algorithm. Setting it to values such as nu=0.1 avoids overfitting in early steps.
#' @param cv A boolean value indicating if cross validation should be performed.
#' @param multicore A boolean value indication if more than one core should be used for the cross validation (for this the parallel package is needed).
#' @param RepMean A numeric value indicating the number of modified datasets used the estimate.
#' @param switch.to A numeric value indicating the number of the first risk up to which the observations can be switched. It sets the value switch.t. If switch.t is specified, switch.to will be ignored.
#' @param switch.t A numeric value indicating the time point up to which the observations can be switched. If this value is not NULL, switch.to will be ignored.
#' @param seed.start A numeric value indicating the first seed value. 
#' @param trace A boolean value indicating if additonal information should be printed out during the process.
#' @param switches A numeric value indicating the number of switches to be made.
#' @param trace A boolean value indicating if additonal information should be printed to the console (=TRUE)
#' @param part A numeric value between 0 and 1 indicating the fraction which should be used in the resampling approach (=0.632=
#' @param RepCb A numeric value indicating the number of Resampling data sets (=100)
#' @param alpha A numeric value between 0 and 1 indicating the type 1 error (=0.05)
#' @param cv_est A boolean value indicating if cross validation should be performed for the estimates seperately (=TRUE)
#' @return An object of type ciMeanPseudoBoost containing the estimates and the performed boosting step number.
#' @export 
cbMeanPseudoBoost.default <- function(data,xmat,times,stepno=100,maxstepno=100,nu=0.1,cv=TRUE,multicore=FALSE,RepMean=50,switch.to=150,switch.t=NULL,seed.start=NULL,trace=TRUE,switches=NULL,part=0.632,RepCb=100,alpha=0.05,cv_est=TRUE,step.point=0.001,...){
  
  obs.time <- data[[1]]
  status <- data[[2]]
  
  if(RepMean==0){
    RepMean <- 1
  }
  
  if(!is.null(switch.t)){
    if(trace) cat("The argument switch.t is used instead of the argument switch.to.\n")
    if(switch.t==0){
      switch.to <- 0
    } else{
      switch.to <- rev(which(sort(obs.time[status==1])<=switch.t))[1]
    }
  } else{
    if(switch.to == 0){
      switch.t <- 0
    } else{
      switch.t <- sort(obs.time[status==1])[switch.to]
    }
  }
  
  if(trace) cat("For smoothing the first observation times with a observation time <= ",switch.t," will be switched.\n This corresponds to switching up to the ",switch.to," observation time of the risk 1.\n")
  
  
  if(is.null(seed.start)){
    seed.start <- round(runif(1,0,12374914))
  }
  
  if(length(stepno) != RepMean){
    stepno=rep(stepno[1],RepMean)
  }
  
  if(is.null(switches)){
    switches <- 1000
  }
  
  res.mean <- meanPseudoBoost(data,xmat=xmat,times=times,stepno=stepno,maxstepno=maxstepno,nu=nu,cv=cv_est,multicore=multicore,RepMean=RepMean,switch.to=switch.to,seed.start=seed.start,trace=FALSE,switches=switches)
  
  res.mean.ci <- list()
  n <- nrow(xmat)
  nsub <- ceiling(part*n)
  sub <- seq(n)
  
  sample_quantile <- function(x,alpha){
    x <- sort(x)
    B <- length(x)
    qL <- alpha/2
    qU <- 1-qL
    
    posL <- B*qL
    
    if(round(posL)!=posL) { 
      B <- B+1
      posL <- floor(B*qL)
    }
    
    if(posL==0) posL=1
    
    posU <- B -posL
    
    if(alpha<=0){
      posL <- 1
      posU <- B
    }
    
    ci <- c(x[posL],x[posU])
    names(ci) <- c(as.character(qL),as.character(qU))
    
    return(ci)
  }
  
  if (!cv){
    stepno <- res.mean$stepno
  }
  
  for(j in 1:RepCb){
    cat("\nResample: ",j,"\n")
    set.seed(seed.start + (j+1)*RepMean*2005)
    seed.start.Mean <- round(runif(1,0,54648721))
    
    choice <- sample(sub,nsub,replace=FALSE)
    
    xmatsub <- xmat[choice,]
    
    obs.timesub <- obs.time[choice]
    statussub <- status[choice]
    
    
    
    zwischen <- meanPseudoBoost(as.data.frame(cbind(obs.timesub,statussub)),xmat=xmatsub,times=times,stepno=stepno,maxstepno=maxstepno,nu=nu,cv=cv,multicore=multicore,RepMean=RepMean,switch.to=switch.to,switch.t=switch.t,seed.start=seed.start.Mean,trace=FALSE,switches=switches)
    
    
    if(j == 1){
      res.zwischen <- zwischen
    } else {    
      for (xindex in 1:(ncol(xmat)+1)) {
        res.zwischen[[xindex]] <- cbind(res.zwischen[[xindex]],zwischen[[xindex]])
      }
    }       
  }
  for (xindex in 1:(ncol(xmat)+1)) { 
    res.mean.ci[[xindex]] <- cbind(est_mean=res.mean[[xindex]],t(apply(res.zwischen[[xindex]],MARGIN=1,FUN=sample_quantile,alpha=alpha)))
  }
  
  #Reduce Pointwise alpha till overall alpha fits
  alpha.point <- rep(alpha,ncol(xmat)+1)
  alpha.res <- rep(alpha,ncol(xmat+1))
  for (xindex in 1:(ncol(xmat)+1)) { 
    low <- res.zwischen[[xindex]]<res.mean.ci[[xindex]][,2] #Resample estimate smaller than lower ci
    up <- res.zwischen[[xindex]]>res.mean.ci[[xindex]][,3] #Resample estimate smaller than upper ci
    interim <- sum(apply(low | up,MARGIN=2,FUN=any))/RepCb #Fraction where resample estiamte is not in ci at at least one timepoint
    if(interim <= alpha){
      unvalid <- FALSE
    }else{
      unvalid <- TRUE
    }
    while(unvalid & alpha.point[xindex]>0){
      alpha.point[xindex] <- alpha.point[xindex] - step.point
      res.mean.ci[[xindex]] <- cbind(est_mean=res.mean[[xindex]],t(apply(res.zwischen[[xindex]],MARGIN=1,FUN=sample_quantile,alpha=alpha.point[xindex])))
      
      low <- res.zwischen[[xindex]]<res.mean.ci[[xindex]][,2]
      up <- res.zwischen[[xindex]]>res.mean.ci[[xindex]][,3]
      interim <- sum(apply(low | up,MARGIN=2,FUN=any))/RepCb
      if(interim <= alpha){
        unvalid <- FALSE
      }
    }
    
    if(alpha.point[xindex]<=0 | interim<=0){
      cat("Confidence band for Beta ",xindex-1," consists of all resample estimates! Try reducing step.point.\nAlpha: ",interim,"\nAlpha.point: ",alpha.point[xindex],"\n")
    } else{
      if(trace) cat("Beta ",xindex-1,"\n Alpha: ",interim,"\n Alpha.point: ",alpha.point[xindex],"\n")
    }
    alpha.res[xindex] <- interim 
  }
  
  res.mean.ci[[xindex+1]] <- stepno
  res.mean.ci[[xindex+2]] <- alpha.point
  res.mean.ci[[xindex+3]] <- alpha.res
  names(res.mean.ci) <- c(names(res.mean),"alpha.point","alpha.res")
  
  
  class(res.mean.ci) <- "cbMeanPseudoBoost"
  return(res.mean.ci)
  
  #Idee: Verringere alpha_punktweise solange, bis 1/B * Summe(1[ein Schätzer außerhalb der punktweisen Intervalle]) etwa alpha_simultan
}


cbMeanPseudoBoost.simulation <- function(n,times, stepno=100,maxstepno=100,nu=0.1,cv=TRUE,multicore=FALSE,RepMean=50,switch.to=150,seed.start=NULL,trace=TRUE,RepCb=100,alpha=0.05,cv_est=TRUE,switches=NULL,step.point=0.001,...){
  
  if(RepMean==0){
    RepMean <- 1
  }
  
  
  
  
  if(is.null(seed.start)){
    seed.start <- round(runif(1,0,12374914))
  }
  
  if(length(stepno) != RepMean){
    stepno=rep(stepno[1],RepMean)
  }
  
  
  
  #res.mean <- meanPseudoBoost(data,xmat=xmat,times=times,stepno=stepno,maxstepno=maxstepno,nu=nu,cv=cv_est,multicore=multicore,RepMean=RepMean,switch.to=switch.to,seed.start=seed.start,trace=FALSE,switches=switches)
  
  res.mean.ci <- list()
  
  sample_quantile3 <- function(x,alpha){
    x <- sort(x)
    B <- length(x)
    qL <- alpha/2
    qU <- 1-qL
    
    posL <- B*qL
    
    if(round(posL)!=posL) { 
      B <- B+1
      posL <- floor(B*qL)
    }
    
    if(posL==0) posL=1
    
    posU <- B -posL
    
    if(alpha<=0){
      posL <- 1
      posU <- B
    }
    
    ci <- c(median(x),x[posL],x[posU])
    names(ci) <- c("Median",as.character(qL),as.character(qU))
    
    return(ci)
  }
  
  
  for(j in 1:RepCb){
    cat("\nSimulation: ",j,"\n")
    
    set.seed(seed.start + j*RepMean*2005)
    
    p <- 10
    beta <- c(c(1,-1),rep(0,p-2))
    beta2 <- c(c(0,0,1,-1),rep(0,p-4))
    
    xmat <- matrix(rnorm(n*p),n,p)
    real.time <- -(log(runif(n)))/(1*exp(drop(xmat %*% beta)))
    real.time2 <- -(log(runif(n)))/(1*exp(drop(xmat %*% beta2)))
    real.time <- ifelse(real.time > 0.1,0.1+real.time2,real.time)
    cens.time <- rexp(n,rate=1/1)
    second.time <- rexp(n,rate=1/1)
    status <- ifelse(real.time <= cens.time,1,0)
    obs.time <- ifelse(real.time <= cens.time,real.time,cens.time)
    
    compete.status <- ifelse(obs.time <= second.time,status,2)
    compete.time <- ifelse(obs.time <= second.time,obs.time,second.time)
    
    
    obs.time <- compete.time
    status <- compete.status
    
    if(switch.to == 0){
      switch.t <- 0
    } else{
      switch.t <- sort(obs.time[status==1])[switch.to]
    }
    
    if(is.null(switches)){
      switches <- 100*max(which(sort(obs.time[status==1])<=switch.t))
    }
    
    if(trace) cat("For smoothing the first observation times with a observation time <= ",switch.t," will be switched.\n This corresponds to switching up to the ",switch.to," observation time of the risk 1.\n")
    
    
    seed.start.Mean <- round(runif(1,0,54648721))
    
    
    zwischen <- meanPseudoBoost(as.data.frame(cbind(obs.time,status)),xmat=xmat,times=times,stepno=stepno,maxstepno=maxstepno,nu=nu,cv=cv,multicore=multicore,RepMean=RepMean,switch.to=switch.to,switch.t=switch.t,seed.start=seed.start.Mean,trace=FALSE,switches=switches)
    
    
    if(j == 1){
      res.zwischen <- zwischen
    } else {    
      for (xindex in 1:(ncol(xmat)+1)) {
        res.zwischen[[xindex]] <- cbind(res.zwischen[[xindex]],zwischen[[xindex]])
      }
    }       
  }
  for (xindex in 1:(ncol(xmat)+1)) { 
    res.mean.ci[[xindex]] <- cbind(t(apply(res.zwischen[[xindex]],MARGIN=1,FUN=sample_quantile3,alpha=alpha)))
  }
  
  #Reduce Pointwise alpha till overall alpha fits
  alpha.point <- rep(alpha,ncol(xmat)+1)
  alpha.res <- rep(alpha,ncol(xmat+1))
  for (xindex in 1:(ncol(xmat)+1)) { 
    low <- res.zwischen[[xindex]]<res.mean.ci[[xindex]][,2] #Resample estimate smaller than lower ci
    up <- res.zwischen[[xindex]]>res.mean.ci[[xindex]][,3] #Resample estimate smaller than upper ci
    interim <- sum(apply(low | up,MARGIN=2,FUN=any))/RepCb #Fraction where resample estiamte is not in ci at at least one timepoint
    if(interim <= alpha){
      unvalid <- FALSE
    }else{
      unvalid <- TRUE
    }
    while(unvalid & alpha.point[xindex]>0){
      alpha.point[xindex] <- alpha.point[xindex] - step.point
      res.mean.ci[[xindex]] <- cbind(t(apply(res.zwischen[[xindex]],MARGIN=1,FUN=sample_quantile3,alpha=alpha.point[xindex])))
      
      low <- res.zwischen[[xindex]]<res.mean.ci[[xindex]][,2]
      up <- res.zwischen[[xindex]]>res.mean.ci[[xindex]][,3]
      interim <- sum(apply(low | up,MARGIN=2,FUN=any))/RepCb
      if(interim <= alpha){
        unvalid <- FALSE
      }
    }
    
    if(alpha.point[xindex]<=0 | interim<=0){
      cat("Confidence band for Beta ",xindex-1," consists of all resample estimates! Try reducing step.point.\nAlpha: ",interim,"\nAlpha.point: ",alpha.point[xindex],"\n")
    } else{
      if(trace) cat("Beta ",xindex-1,"\n Alpha: ",interim,"\n Alpha.point: ",alpha.point[xindex],"\n")
    }
    alpha.res[xindex] <- interim 
  }
  
  res.mean.ci[[xindex+1]] <- stepno
  res.mean.ci[[xindex+2]] <- alpha.point
  res.mean.ci[[xindex+3]] <- alpha.res
  res.mean.ci[[xindex+4]] <- times
  names(res.mean.ci) <- c(names(res.zwischen),"alpha.point","alpha.res","evaluation.times")
  
  
  class(res.mean.ci) <- "cbMeanPseudoBoost"
  return(res.mean.ci)
  
  #Idee: Verringere alpha_punktweise solange, bis 1/B * Summe(1[ein Schätzer außerhalb der punktweisen Intervalle]) etwa alpha_simultan
}

#' Plot confidence bands and estimates.
#' 
#' @param object Object of type cbMeanPseudoBoost
#' @param est A numeric vector indicating the index of the estimates to be plotted.
#' @param ci A numeric vector indicating the index of the estimates to be plotted with confidence bands.
#' @param trans A boolean vector indicating if the results should be transformed (trans=TRUE => Plot exponential of the estimates).
#' @param name A string value indicating the name of the resulting PDF. E.g. name="results.pdf"
#' @return A PDF document with the name "name".
#' @export 
plot.cbMeanPseudoBoost <- function(object,est, ci,alpha=0.05,trans=TRUE,name="results.pdf"){
  pdf(name)
  if (length(est) != 0) {
    if(trans==TRUE){plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(eval.times)),ylim=c(0,2))}
    else{plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(eval.times)),ylim=c(-2,2))}
    
    if(trans==TRUE){abline(h=1, col = "lightgray")}
    else {abline(h=0, col = "lightgray")}
    
    linetyp <- 1
    
    for (i in 1:length(est)) {
      if(trans==TRUE){lines(eval.times,exp(object[[est[i]]][,1]),lty=linetyp)}
      else {lines(eval.times,object[[est[i]]][,1],lty=linetyp)}
      linetyp <- linetyp+1
    }
    legend("topright",legend=names(as.list(object[est])),lty=seq(linetyp),bty="n")
    title("Estimater for Parameters")
    
  }
  
  if (length(ci) != 0){
    for (i in 1:length(ci)) {
      if(trans==TRUE){plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(eval.times)),ylim=c(0,2))}
      else{plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(eval.times)),ylim=c(-2,2))}
      
      if(trans==TRUE){abline(h=1, col = "lightgray")}
      else {abline(h=0, col = "lightgray")}
      
      linetyp=2
      
      for (j in c(1,2,3)){
        if(trans==TRUE){lines(eval.times,exp(object[[ci[i]]][,j]),lty=min(linetyp,j))}
        else{lines(eval.times,object[[ci[i]]][,j],lty=min(linetyp,j))}
      }
      
      if (ci[i] == 1) {
        legend("topright",legend=c("Estimator of Intercept","confidencebounds"),lty=seq(linetyp),bty="n")
      } else {
        legend("topright",legend=c(paste("Estimator of beta",ci[i]-1,sep="_"),"confidencebounds"),lty=seq(linetyp),bty="n")
      }
      title(paste("Estimater and ",(1-alpha)*100,"%-confidence band for ",names(as.list(object[ci[i]])),sep=""))
      
    }  
  }  
  dev.off()
  
}