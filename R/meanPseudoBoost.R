meanPseudoBoost <- function(object,...){
  UseMethod("meanPseudoBoost",object)
}

meanPseudoBoost.default <- function(data,xmat,times,stepno=100,maxstepno=100,nu=0.1,cv=TRUE,multicore=FALSE,RepMean=50,switch.to=150,switch.t=NULL,seed.start=NULL,trace=TRUE,switches = 1000,...){
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
    for (xindex in 1:(ncol(xmat)+1)) {
      res.mean[[xindex]] <- apply(reslist[[xindex]],MARGIN=1,FUN=mean)
      
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
  
  class(res.mean) <- "meanPseudoBoost"
  return(res.mean)
}

plot.meanPseudoBoost <- function(object,eval.times,est,comb=TRUE,trans=TRUE,name="results.pdf"){
  pdf(name)
  if (comb) {
    if(trans==TRUE){plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(eval.times)),ylim=c(0,2))}
    else{plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(eval.times)),ylim=c(-2,2))}
    
    if(trans==TRUE){abline(h=1, col = "lightgray")}
    else {abline(h=0, col = "lightgray")}
    
    linetyp <- 1
    
    for (i in 1:length(est)) {
      if(trans==TRUE){lines(eval.times,exp(object[[est[i]]]),lty=linetyp)}
      else {lines(eval.times,object[[est[i]]],lty=linetyp)}
      linetyp <- linetyp+1
    }
    legend("topright",legend=names(as.list(object[est])),lty=seq(linetyp),bty="n")
    title("Estimater for Parameters")
    
  }
  
  else {
    for (i in 1:length(est)) {
      if(trans==TRUE){plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(eval.times)),ylim=c(0,2))}
      else{plot(0,type="n",xlab="time",ylab="coefficient",xlim=c(0,1.3*max(eval.times)),ylim=c(-2,2))}
      
      if(trans==TRUE){abline(h=1, col = "lightgray")}
      else {abline(h=0, col = "lightgray")}
      
      linetyp=2
      
      for (j in c(1)){
        if(trans==TRUE){lines(eval.times,exp(object[[est[i]]]),lty=min(linetyp,j))}
        else{lines(eval.times,object[[est[i]]],lty=min(linetyp,j))}
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