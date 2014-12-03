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
  class(CIF) <- "etmCIF"
  return(CIF)
}

PseudoCIF <- function(data,times){
  require(etm)
  data <- as.data.frame(cbind(entry=data[[1]],exit=data[[2]],status=data[[3]]))
  n_data <- nrow(data)
  ids <- as.list(seq(1:n_data))
  
  #Versuch mit etm:
  #states <- data$status
  #states[states==0] <- "cens" 
  #states <- factor(states)
  #state.names<-levels(factor(data$status))
  #tra <- matrix(FALSE, length(state.names), length(state.names))
  #tra[1, 2:length(state.names)] <- TRUE
  #all <- n_data * trprob(etmCIF_kurz(data,n_data)[[1]], "0 1", times)
  #pseudoCIF <- t(sapply(ids,FUN=function(arg) {all - (n_data-1) * trprob(etmCIF_kurz(data[-arg,],n_data-1)[[1]], "0 1", times)}))
  
  #mit etmCIF:
  #all <- n_data * trprob(etmCIF(Surv(entry, exit, status != 0) ~ 1, data, etype = status, failcode = 1)[[1]], "0 1", times)
  #pseudoCIF <- t(sapply(ids,FUN=function(arg) {all - (n_data-1) * trprob(etmCIF(Surv(entry, exit, status != 0) ~ 1, data[-arg,], etype = status, failcode = 1)[[1]], "0 1", times)}))
  
  #mit etmCIF_kurz:
  all <- n_data * trprob(etmCIF_kurz(data, n_data)[[1]], "0 1", times)
  pseudoCIF <- t(sapply(ids,FUN=function(arg) {all - (n_data-1) * trprob(etmCIF_kurz(data[-arg,], n_data-1)[[1]], "0 1", times)}))
  
  return(pseudoCIF)
}