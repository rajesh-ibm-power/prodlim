jackknife <- function(object,times,...){
  if (object$model=="survival")
    jackknife.survival(object=object,times=times,...)
  else if (object$model=="competing.risks")
    jackknife.competing.risks(object=object,times=times,...)
  else stop("No method for jackknifing this object.")
}

leaveOneOut <- function(object,times,...){
  if (object$model=="survival")
    leaveOneOut.survival(object=object,times=times,...)
  else if (object$model=="competing.risks")
    leaveOneOut.competing.risks(object=object,times=times,...)
  else stop("No method for jackknifing this object.")
}

jackknife.survival <- function(object,times,...){
  S <- predict(object,times=times,newdata=object$model.response)
  Sk <- leaveOneOut.survival(object,times,...)
  N <- NROW(Sk)
  Jk <- t(N*S-t((N-1)*Sk))
  ##   Jk <- N*S-(N-1)*Sk
  Jk
}

jackknife.competing.risks <- function(object,times,cause,...){
  F <- predict(object,times=times,newdata=object$model.response,cause=cause)
  Fk <- leaveOneOut.competing.risks(object,times,cause,...)
  N <- NROW(Fk)
  ##   Jk <- N*F-(N-1)*Fk
  Jk <- t(N*F-t((N-1)*Fk))
  Jk
}


leaveOneOut.survival <- function(object,times,lag=0,...){
  stopifnot(object$covariate.type==1)
  ##
  time <- object$time
  Y <- object$n.risk
  D <- object$n.event
  Y <- Y[D>0]
  time <- time[D>0]
  D <- D[D>0]
  NU <- length(time)
  obstimes <- object$model.response[,"time"]
  status <- object$model.response[,"status"]
  N <- length(obstimes)
  ##
  S <- predict(object,times=time,newdata=object$model.response)
  ## idea: find the at-risk set for pseudo-value k by
  ##       substracting 1 in the period where subj k is
  ##       at risk. need the position of obstime.k in time ...
  ##   pos <- match(obstimes,time)
  pos <- sindex(jump.times=time,eval.times=obstimes)
  ##
  loo <- do.call("rbind",lapply(1:N,function(k){
    Dk <- D
    if (status[k]==1) Dk[pos[k]] <- Dk[pos[k]]-1
    Yk <- Y-c(rep(1,pos[k]),rep(0,NU-pos[k]))
    cumprod(1-Dk/Yk)
  }))
  out <- loo
  if (!missing(times)){
    found <- sindex(jump.times=time,eval.times=times)+1
    if (lag==0)
      out <- cbind(1,out)[,found,drop=TRUE]
    else
      out <- cbind(1,cbind(1,out))[,found,drop=TRUE]
  }
  out
}

leaveOneOut.competing.risks <- function(object,times,cause,...){
  stopifnot(object$covariate.type==1)
  mr <- object$model.response
  states <- states(object)
  if (missing(cause)) {
    C <- 1
    cause <- states[1]
  }
  else{
    C <- match(cause,states,nomatch=0)
    if (length(C)>1 || C==0) stop("Cause must match exactly one of the names of object$n.event.")
  }
  D <- object$n.event[[C]]
  #  it is enough to consider time points where an event occurs
  time <- object$time[D>0]
  Y <- object$n.risk[D>0]
  sFit <- prodlim(Hist(time,status)~1,data=data.frame(unclass(mr)))
  S <- sFit$surv[D>0]
  D <- D[D>0]
  lagSk <- leaveOneOut.survival(sFit,times=time,lag=1)
  NU <- length(time)
  obstimes <- mr[,"time"]
  status <- mr[,"status"]
  E <- events(object)
  N <- length(obstimes)
  ## idea: see leaveOneOut.survival
  pos <- sindex(jump.times=time,eval.times=obstimes)
  loo <- do.call("rbind",lapply(1:N,function(k){
    Dk <- D
    if (status[k]==1 && E[k]==cause) Dk[pos[k]] <- Dk[pos[k]]-1
    Yk <- Y-c(rep(1,pos[k]),rep(0,NU-pos[k]))
    Sk <- as.numeric(lagSk[k,,drop=TRUE])
    Hk <- Dk/Yk
    Fk <- cumsum(Sk*Hk)
    Fk
  }))
  out <- loo
  if (!missing(times)){
    found <- sindex(jump.times=time,eval.times=times)+1
    out <- cbind(1,out)[,found,drop=TRUE]
  }
#  out[NROW(out)]
  out
}


