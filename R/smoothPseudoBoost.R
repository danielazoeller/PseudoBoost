#' Perform smoothed stagewise pseudo-value regression for the first risk in a competing risk setting
#' 
#' Creates RepSmooth modified datasets with the help of randomly changed CIF-values and corresponding observation times. For these dataset the function PseudoBoost will be called seperately and the mean value is taken for each effect estimate at each considered time point.
#' @param data A data.frame containing observation times AND statuses
#' @param xmat A numeric matrix containing the covariate values for all patients.
#' @param times A numeric vector containing the evaluation times.
#' @param stepno A numeric value containing the number of boosting steps to be performed. If you use cross validation (cv=TRUE), this parameter will be ignored.
#' @param maxstepno A numeric value containing the maximal number of boosting steps considered during the cross validation (only used if cv=TRUE).
#' @param nu A numeric value between 0 and 1, the shrinkage parameter for the stagewise regression algorithm. Setting it to values such as nu=0.1 avoids overfitting in early steps.
#' @param cv A boolean value indicating if cross validation should be performed.
#' @param multicore A boolean value indication if more than one core should be used for the cross validation (for this the parallel package is needed).
#' @param RepSmooth A numeric value indicating the number of modified datasets used the estimate.
#' @param smooth_para A numeric value used for smoothing with the beta-distribution (Variance of beta-distribution).
#' @param seed.start A numeric value indicating the first seed value. 
#' @param trace A boolean value indicating if additonal information should be printed out during the process.
#' @return An object of type meanPseudoBoost containing the estimates, the performed boosting step number and the evaluation times.
#' @import prodlim
#' @import etm
#' @export 
smoothPseudoBoost <- function(object,...){
  UseMethod("smoothPseudoBoost",object)
}

#' Perform smoothed stagewise pseudo-value regression for the first risk in a competing risk setting
#' 
#' Creates RepSmooth modified datasets with the help of randomly changed CIF-values and corresponding observation times. For these dataset the function PseudoBoost will be called seperately and the mean value is taken for each effect estimate at each considered time point.
#' @param data A data.frame containing observation times AND statuses
#' @param xmat A numeric matrix containing the covariate values for all patients.
#' @param times A numeric vector containing the evaluation times.
#' @param stepno A numeric value containing the number of boosting steps to be performed. If you use cross validation (cv=TRUE), this parameter will be ignored.
#' @param maxstepno A numeric value containing the maximal number of boosting steps considered during the cross validation (only used if cv=TRUE).
#' @param nu A numeric value between 0 and 1, the shrinkage parameter for the stagewise regression algorithm. Setting it to values such as nu=0.1 avoids overfitting in early steps.
#' @param cv A boolean value indicating if cross validation should be performed.
#' @param multicore A boolean value indication if more than one core should be used for the cross validation (for this the parallel package is needed).
#' @param RepSmooth A numeric value indicating the number of modified datasets used the estimate.
#' @param smooth_para A numeric value used for smoothing with the beta-distribution (Variance of beta-distribution).
#' @param seed.start A numeric value indicating the first seed value. 
#' @param trace A boolean value indicating if additonal information should be printed out during the process.
#' @return An object of type meanPseudoBoost containing the estimates, the performed boosting step number and the evaluation times.
#' @import prodlim
#' @import etm
#' @export 
smoothPseudoBoost.default <- function(data,xmat,times,stepno=100,maxstepno=100,nu=0.1,cv=TRUE,multicore=FALSE,RepSmooth=20,smooth_para = 2,seed.start=NULL,trace=TRUE,...){
  require(prodlim)
  set.seed(seed.start)
  
  obs.time <- data[[1]]
  status <- data[[2]]
    
    
  if(length(stepno) != RepSmooth){
    if(RepSmooth>0){
      stepno=rep(stepno[1],RepSmooth)
    } else{
      stepno=stepno[1]
    }
  }
  
  
  n <- nrow(xmat)
    
  reslist <- list() #One element will contain all resulting values for one beta_i for every time point for every sample
  res.mean <- list() #One element will contain estimate and confidencebounds for one beta_i for every time point
  
  ymatsub.list = vector("list",n)
  
  if(RepSmooth>0){
    etmCIF_kurz <- function(data,n_data){
      require(etm)
      states <- data$status
      states[states==0] <- "cens" 
      states <- factor(states)
      state.names<-levels(factor(data$status))
      tra <- matrix(FALSE, length(state.names), length(state.names))
      tra[1, 2:length(state.names)] <- TRUE
      CIF <- list()
      CIF[[1]] <- etm(data.frame(id = seq(1,n_data),from=rep(0,n_data),to=states,entry=data$entry,exit=data$exit,cov=rep(1,n_data)), state.names=state.names,tra=tra, cens.name="cens",s=0,covariance=FALSE)
      return(CIF)
    }
    
     time_smooth <- matrix(rep(obs.time,n*RepSmooth),n,RepSmooth)
    value <- matrix(rep(0,n*RepSmooth),n,RepSmooth)
    
    cause_1 <- which(status==1)
    
    CIF <- trprob(etmCIF_kurz(data.frame(entry=rep(0,n),exit=obs.time,status=status),n)[[1]], "0 1", obs.time)
    CIF_sort <- data.frame(time = sort(obs.time), CIF=CIF[order(obs.time)])
    
    #on [0,1]:
    #cif_delta_min <- (1-(min(CIF[cause_1])))/(min(CIF[cause_1]))
    #mini <- cif_delta_min/((1+cif_delta_min)^2)
    #on [0,cif_max]
    cif_max <- min(c(1,max(CIF[cause_1]+0.01)))
    cif_delta_min <- (1-(min(CIF[cause_1])/cif_max))/(min(CIF[cause_1])/cif_max)
    mini <- cif_delta_min/((1+cif_delta_min)^2)
    
    #auf[0,cif_max]:
    smooth_scale <- smooth_para/((cif_max)^2)
    if(mini<smooth_scale){
      smooth_scale <- mini-(mini/100)
      smooth_para_new <- smooth_scale*((cif_max)^2)
      cat("ERROR: Variance of beta-distribution (smooth_para) must be smaller. Value has been set to be a little bit smaller: ",smooth_para_new)
    }
    for(i in cause_1){
      #beta distribution on [0,1]
#       cif_delta_i <- (1-CIF[i])/CIF[i]
#       
#       a <- cif_delta_i - (((1+cif_delta_i)^2)*smooth_para)
#       b <- (smooth_para) * ((1+cif_delta_i)^3)
#       p <- a/b
#       
#       q=cif_delta_i * p
#             
#       value[i,] <- rbeta(20,p,q)
      
      #beta distribution on [0,cif_max]
      cif_i <- CIF[i]/cif_max
      cif_delta_i <- ((1-cif_i)/cif_i)
      
      a <- cif_delta_i - (((1+cif_delta_i)^2)*smooth_scale)
      b <- (smooth_scale) * ((1+cif_delta_i)^3)
      p <- a/b
      
      q=cif_delta_i * p
      
      value[i,] <- cif_max * rbeta(20,p,q)
      #mehrere Varianten zur Bestimmung der Zeit mit dieser CIF: 
      #1. Bereich mit letzte CIF<berechnete oder erste CIF>=berechnete
      #2. Start oder Stop-Zeit dieses Bereichs
      #3. Intrapoliert zwischen Sprungzeiten
      for(j in 1:RepSmooth){#Variante 3
        if(any(CIF_sort$CIF>value[i,j]) && any(CIF_sort$CIF<value[i,j])){
          arg <- min(which(CIF_sort$CIF>value[i,j]))
          time_smooth[i,j] <- (CIF_sort$time[arg-1]*((CIF_sort$CIF[arg]-value[i,j])/(CIF_sort$CIF[arg]-CIF_sort$CIF[arg-1]))) + (CIF_sort$time[arg]*((value[i,j]-CIF_sort$CIF[arg-1])/(CIF_sort$CIF[arg]-CIF_sort$CIF[arg-1])))
        } else if (any(CIF_sort$CIF>value[i,j])){
          time_smooth[i,j] <- min(CIF_sort$time[cause_1])
        } else {
          time_smooth[i,j] <- max(CIF_sort$time[cause_1])
        }
        
      }
      
      #Versuch 1
#       for(j in 1:RepSmooth){#Variante 3
#         index <- which(CIF_sort[,2]==value[i,j])
#         if(length(index)==0 && any(CIF_sort$CIF>value[i,j]) && any(CIF_sort$CIF<value[i,j])){
#           arg <- min(which(CIF_sort$CIF>value[i,j]))
#           time_smooth[i,j] <- (CIF_sort$time[arg-1]*((CIF_sort$CIF[arg]-value[i,j])/(CIF_sort$CIF[arg]-CIF_sort$CIF[arg-1]))) + (CIF_sort$time[arg]*((value[i,j]-CIF_sort$CIF[arg-1])/(CIF_sort$CIF[arg]-CIF_sort$CIF[arg-1])))
#         } else if (any(CIF_sort$CIF>value[i,j]) && any(CIF_sort$CIF<value[i,j])){
#           arg <- min(which(CIF_sort$CIF>value[i,j]))
#           time_smooth[i,j] <- (CIF_sort$time[arg-1]*((CIF_sort$CIF[arg]-value[i,j])/(CIF_sort$CIF[arg]-CIF_sort$CIF[arg-2])))
#         } else {
#           time_smooth[i,j] <- max(CIF_sort$time)
#         }
#         
#       }
    }
  
  
    for (i in 1:RepSmooth) {
      set.seed(seed.start+i*100)
      cat("Replication:",i,"\n")
      
      obs.timesub <- time_smooth[,i]
      ymat <- jackknife(prodlim(Hist(obs.timesub,status) ~ 1), times=times, cause=1)
      if (cv) {
        cv.res <- cv.PseudoBoost(ymat,xmat,maxstepno=maxstepno,nu=nu,trace=FALSE,...)
        stepno[i] <- cv.res$optimal.step
      }
      
      
      ressub <- PseudoBoost(ymat,xmat,stepno=stepno[i],nu=nu,trace=FALSE)
      
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
    if(RepSmooth>1){
      for (xindex in 1:(ncol(xmat)+1)) {
        res.mean[[xindex]] <- apply(reslist[[xindex]],MARGIN=1,FUN=mean)
        
      }
    } else {
      for (xindex in 1:(ncol(xmat)+1)) {
        res.mean[[xindex]] <- reslist[[xindex]]
        
      }
    }
    
  } else {
    set.seed(seed.start+100)
    ymat <- jackknife(prodlim(Hist(obs.time,status) ~ 1), times=times, cause=1)
    
    
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
  
  class(res.mean) <- "smoothPseudoBoost"
  return(res.mean)
}

#' Plot an object of class smoothPseudoBoost. The estimated effects will be plotted against the evaluation times saved in the object. Output as PDF.
#' 
#' @param object An object of class meanPseudoBoost created by smoothPseudoBoost().
#' @param est A numeric vector containing the indices of the estimated to be plotted. est=c(1) will e.g. plot the estimated intercept.
#' @param comb A boolean vector indicating if the results should be combined in one plot (comb=TRUE) or if each estimate should be plotted seperately (comb=FALSE).
#' @param trans A boolean vector indicating if the results should be transformed (trans=TRUE => Plot exponential of the estimates).
#' @param name A string value indicating the name of the resulting PDF. E.g. name="results.pdf"
#' @return A PDF document with the name "name".
#' @export 
plot.smoothPseudoBoost <- function(object,est,comb=TRUE,trans=TRUE,name="results.pdf"){
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