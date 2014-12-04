#' Perform smoothed stagewise pseudo-value regression for the first risk in a competing risk setting
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
#' @param switch.to A numeric value indicating the number of the first risk up to which the observations can be switched. It sets the value switch.t. If switch.t is specified, switch.to will be ignored. If switch.to=0, switch.t will be set to zero.
#' @param switch.t A numeric value indicating the time point up to which the observations can be switched. If this value is not NULL, switch.to will be ignored. If switch.t=0, there will be no switching.
#' @param seed.start A numeric value indicating the first seed value. 
#' @param trace A boolean value indicating if additonal information should be printed out during the process.
#' @param switches A numeric value indicating the number of switches to be made. Setting switches=0 will suppress smoothing. If switches is NULL, 1000 switches will be performed.
#' @return An object of type meanPseudoBoost containing the estimates, the performed boosting step number and the evaluation times.
#' @import prodlim
#' @export 
meanPseudoBoost <- function(object,...){
  UseMethod("meanPseudoBoost",object)
}

#' Perform smoothed stagewise pseudo-value regression for the first risk in a competing risk setting
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
#' @param switch.to A numeric value indicating the number of the first risk up to which the observations can be switched. It sets the value switch.t. If switch.t is specified, switch.to will be ignored. If switch.to=0, switch.t will be set to zero.
#' @param switch.t A numeric value indicating the time point up to which the observations can be switched. If this value is not NULL, switch.to will be ignored. If switch.t=0, there will be no switching.
#' @param seed.start A numeric value indicating the first seed value. 
#' @param trace A boolean value indicating if additonal information should be printed out during the process.
#' @param switches A numeric value indicating the number of switches to be made. Setting switches=0 will suppress smoothing. If switches is NULL, 1000 switches will be performed.
#' @return An object of type meanPseudoBoost containing the estimates, the performed boosting step number and the evaluation times.
#' @import prodlim
#' @export 
meanPseudoBoost.default <- function(data,xmat,times,stepno=100,maxstepno=100,nu=0.1,cv=TRUE,multicore=FALSE,RepMean=50,switch.to=150,switch.t=NULL,seed.start=NULL,trace=TRUE,switches = NULL,...){
  require(prodlim)
  
  obs.time <- data[[1]]
  status <- data[[2]]
  
  
  status <- status[order(obs.time)]
  xmat <- xmat[order(obs.time),]
  obs.time <- obs.time[order(obs.time)]
  
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
  
  if(is.null(switches)){
    switches <- 1000
  }
  
  if(trace) cat("For smoothing the first observation times with a observation time <= ",switch.t," will be switched.\n This corresponds to switching up to the ",switch.to," observation time of the risk 1.\n")
  
  if(is.null(seed.start)){
    seed.start <- round(runif(1,0,749853798))
  }
  
  
  
  ymat <- jackknife(prodlim(Hist(obs.time,status) ~ 1), times=times, cause=1)
  
  if(RepMean==0){
    RepMean <- 1
  }
  
  if(length(stepno) != RepMean){
    stepno=rep(stepno[1],RepMean)
  }
  
  
  n <- nrow(xmat)
  sub <- seq(n)
  ids <- seq(1:n)
  ids.sub <- ids[1:length(which(obs.time<=switch.t))]
  ids.rest <- ids[-ids.sub]
  
  if(switch.to==0){
    switches <- 0
  }
  
  reslist <- list() #One element will contain all resulting values for one beta_i for every time point for every sample
  res.mean <- list() #One element will contain estimate and confidencebounds for one beta_i for every time point
  
  ymatsub.list = vector("list",n)
  
  if(switch.t != 0){
    for (i in 1:RepMean) {
      cat("Replication:",i,"\n")
      set.seed(seed.start+i*100)
      
      
      pairs <- sample(seq(length(ids.sub)-1),switches,replace=TRUE)
      obs.timesub <- obs.time
      
      if(length(pairs)>0){
        for(pairs.id in 1:length(pairs)){
          zwischen <- obs.timesub[pairs[pairs.id]]
          obs.timesub[pairs[pairs.id]] <- obs.timesub[pairs[pairs.id]+1]
          obs.timesub[pairs[pairs.id]+1] <- zwischen
        }
      }
      
      statussub <- status
      
      xmatsub <- xmat
      
      
      
      ymatsub <- jackknife(prodlim(Hist(obs.timesub,statussub) ~ 1), times=times, cause=1)
      if (cv) {
        cv.res <- cv.PseudoBoost(ymatsub,xmatsub,maxstepno=maxstepno,nu=nu,trace=FALSE,...)
        stepno[i] <- cv.res$optimal.step
      }
      
      
      ressub <- PseudoBoost(ymatsub,xmatsub,stepno=stepno[i],nu=nu,trace=FALSE)
      
      l = length(ressub$coefficients)
      ressub <- ressub$coefficients[[l]]
      
      if(cv) cat("Number of Boosting steps: ",l,"\n")
      
      
      if (i ==1) {
        for (xindex in 1:(ncol(xmat)+1)) {
          reslist[[xindex]] <- ressub[,xindex]
        }
      } else {
        for (xindex in 1:(ncol(xmat)+1)) {
          reslist[[xindex]] <- cbind(reslist[[xindex]],ressub[,xindex])
        }
      }
    }
    
    #reslist[[xindex+1]] <- stepno
    
    
    
    #   if(switch.to==0){ 
    #     for (xindex in 1:(ncol(xmat)+1)) {
    #       res.mean[[xindex]] <- reslist[[xindex]]
    #       
    #       
    #       
    #       if (xindex == 1){
    #         rescinames[xindex] <- "Intercept"
    #       } else {
    #         rescinames[xindex] <- paste("beta",(xindex-1),sep="_")
    #       }
    #     }
    #   }else{
    if(RepMean>1){
      for (xindex in 1:(ncol(xmat)+1)) {
        res.mean[[xindex]] <- apply(reslist[[xindex]],MARGIN=1,FUN=mean)
        
      }
    } else {
      for (xindex in 1:(ncol(xmat)+1)) {
        res.mean[[xindex]] <- reslist[[xindex]]
        
      }
    }
    
  } else {
    set.seed(seed.start)
    
    if (cv) {
      cv.res <- cv.PseudoBoost(ymat,xmat,maxstepno=maxstepno,nu=nu,trace=FALSE,...)
      stepno[1] <- cv.res$optimal.step
    }
    
    ressub <- PseudoBoost(ymat,xmat,stepno=stepno[1],nu=nu,trace=FALSE)
    
    l = length(ressub$coefficients)
    ressub <- ressub$coefficients[[l]]
    
    if(cv) cat("Number of Boosting steps: ",l,"\n")
    
    
    for (xindex in 1:(ncol(xmat)+1)) {
      res.mean[[xindex]] <- ressub[,xindex]
    }
    
  }
  
  rescinames <- rep("stepno",ncol(xmat) + 2)  
  
  for (xindex in 1:(ncol(xmat)+1)){
    if (xindex == 1){
      rescinames[xindex] <- "Intercept"
    } else {
      rescinames[xindex] <- paste("beta",(xindex-1),sep="_")
    }
  }
  #  }
  res.mean[[xindex+1]] <- stepno
  
  names(res.mean) <- rescinames
  
  res.mean[[xindex+2]] <- times
  names(res.mean)[xindex+2] <- "evaluationTimes"
  
  class(res.mean) <- "meanPseudoBoost"
  return(res.mean)
}

#' Plot an object of class meanPseudoBoost. Output as PDF.
#' 
#' @param object An object of class meanPseudoBoost created by meanPseudoBoost().
#' @param est A numeric vector containing the indices of the estimated to be plotted. est=c(1) will e.g. plot the estimated intercept.
#' @param comb A boolean vector indicating if the results should be combined in one plot (comb=TRUE) or if each estimate should be plotted seperately (comb=FALSE).
#' @param trans A boolean vecotr indicating if the results should be transformed (trans=TRUE => Plot exponential of the estimates).
#' @param name A string value indicating the name of the resulting PDF. E.g. name="results.pdf"
#' @return A PDF document with the name "name".
#' @export 
plot.meanPseudoBoost <- function(object,est,comb=TRUE,trans=TRUE,name="results.pdf"){
  pdf(name)
  times <- object$evaluationTimes
  if (comb) {
    if(trans==TRUE){plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(times)),ylim=c(0,2))}
    else{plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(times)),ylim=c(-2,2))}
    
    if(trans==TRUE){abline(h=1, col = "lightgray")}
    else {abline(h=0, col = "lightgray")}
    
    linetyp <- 1
    
    for (i in 1:length(est)) {
      if(trans==TRUE){lines(times,exp(object[[est[i]]]),lty=linetyp)}
      else {lines(times,object[[est[i]]],lty=linetyp)}
      linetyp <- linetyp+1
    }
    legend("topright",legend=names(as.list(object[est])),lty=seq(linetyp),bty="n")
    title("Estimater for Parameters")
    
  }
  
  else {
    for (i in 1:length(est)) {
      if(trans==TRUE){plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(times)),ylim=c(0,2))}
      else{plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(times)),ylim=c(-2,2))}
      
      if(trans==TRUE){abline(h=1, col = "lightgray")}
      else {abline(h=0, col = "lightgray")}
      
      linetyp=2
      
      for (j in c(1)){
        if(trans==TRUE){lines(times,exp(object[[est[i]]]),lty=min(linetyp,j))}
        else{lines(times,object[[est[i]]],lty=min(linetyp,j))}
      }
      
      if (est[i] == 1) {
        legend("topright",legend=c("Estimator of Intercept","confidencebounds"),lty=seq(linetyp),bty="n")
      } else {
        legend("topright",legend=c(paste("Estimator of beta",est[i]-1,sep="_"),"confidencebounds"),lty=seq(linetyp),bty="n")
      }
      title(paste("Estimater for ",names(as.list(object[est[i]])),sep=""))
      
    }  
  }  
  dev.off()
  
}