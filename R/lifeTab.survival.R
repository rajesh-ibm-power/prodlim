lifeTab.survival <- function(object,times,newdata,stats,intervals=FALSE,percent=TRUE,showTime=TRUE){
  # {{{ get the indices
  IndeX <- predict(object,newdata=newdata,level.chaos=0,times=times,type="list")
  # }}}
  # {{{ times
  times <- IndeX$times
  Ntimes <- IndeX$dimensions$time
  pindex <- IndeX$indices$time
  # }}}
  # {{{ covariate strata
  Nstrata <- IndeX$dimensions$strata
  findex <- IndeX$indices$strata
  # }}}
  # {{{ stats
  if (missing(stats) || ((!missing(stats)) && is.null(stats)))
    stats <- list(c("n.event",0),c("n.lost",0))
  else{
    stats <- c(list(c("n.event",0),c("n.lost",0)),stats)
  }
  # }}}
  # {{{ no. at atrisk, events, and censored
  if (intervals==FALSE){
    if (is.null(object$clustervar)){
      ## only one column for n.risk
      xxx <- .C("life_table",pred.nrisk=integer(Ntimes*Nstrata),pred.nevent=integer(Ntimes*Nstrata),pred.nlost=integer(Ntimes*Nstrata),nrisk=as.integer(object$n.risk),nevent=as.integer(object$n.event),nlost=as.integer(object$n.lost),as.double(times),as.double(object$time),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(Nstrata),as.integer(Ntimes),intervals=as.integer(FALSE),NAOK=FALSE,PACKAGE="prodlim")
      out <- data.frame(n.risk=xxx$pred.nrisk)
    }
    else{
      xxx <- .C("life_table",pred.nrisk=integer(Ntimes*Nstrata),pred.nevent=integer(Ntimes*Nstrata),pred.nlost=integer(Ntimes*Nstrata),nrisk=as.integer(object$n.risk[,1]),nevent=as.integer(object$n.event[,1]),nlost=as.integer(object$n.lost[,1]),as.double(times),as.double(object$time),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(Nstrata),as.integer(Ntimes),intervals=as.integer(FALSE),NAOK=FALSE,PACKAGE="prodlim")
      out <- data.frame(n.risk=xxx$pred.nrisk)
      for (cv in 1:length(object$clustervar)){
        yyy <- .C("life_table",pred.nrisk=integer(Ntimes*Nstrata),pred.nevent=integer(Ntimes*Nstrata),pred.nlost=integer(Ntimes*Nstrata),nrisk=as.integer(object$n.risk[,1+cv]),nevent=as.integer(object$n.event[,1+cv]),nlost=as.integer(object$n.lost[,1+cv]),as.double(times),as.double(object$time),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(Nstrata),as.integer(Ntimes),intervals=as.integer(FALSE),NAOK=FALSE,PACKAGE="prodlim")}
      outCV <- data.frame(n.risk=yyy$pred.nrisk)
      names(outCV) <- paste(object$clustervar,names(outCV))
      out <- cbind(out,outCV)
    }
  }
  # }}}
  # ------------------------Intervals---------------------------
  else{
    #,----
    #|      get no. at risk at the left limit of the interval
    #|      and count events and censored excluding the left limit
    #`----
    lagTimes <- c(min(min(object$time),0)-.1 , times[-length(times)])
    if (is.null(object$clustervar)){
      ## only one column in n.event and n.risk
      xxx <- .C("life_table",pred.nrisk=integer(Ntimes*Nstrata),pred.nevent=integer(Ntimes*Nstrata),pred.nlost=integer(Ntimes*Nstrata),nrisk=as.integer(object$n.risk),nevent=as.integer(object$n.event),nlost=as.integer(object$n.lost),as.double(times),as.double(object$time),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(Nstrata),as.integer(Ntimes),intervals=as.integer(TRUE),NAOK=FALSE,PACKAGE="prodlim")
      out <- data.frame(n.risk=xxx$pred.nrisk,n.event=xxx$pred.nevent,n.lost=xxx$pred.nlost)
      lagxxx <- .C("life_table",pred.nrisk=integer(Ntimes*Nstrata),pred.nevent=integer(Ntimes*Nstrata),pred.nlost=integer(Ntimes*Nstrata),nrisk=as.integer(object$n.risk),nevent=as.integer(object$n.event),nlost=as.integer(object$n.lost),as.double(lagTimes),as.double(object$time),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(Nstrata),as.integer(Ntimes),intervals=as.integer(TRUE),NAOK=FALSE,PACKAGE="prodlim")
      out$n.risk <- lagxxx$pred.nrisk
    }
    else{
      xxx <- .C("life_table",pred.nrisk=integer(Ntimes*Nstrata),pred.nevent=integer(Ntimes*Nstrata),pred.nlost=integer(Ntimes*Nstrata),nrisk=as.integer(object$n.risk[,1]),nevent=as.integer(object$n.event[,1]),nlost=as.integer(object$n.lost[,1]),as.double(times),as.double(object$time),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(Nstrata),as.integer(Ntimes),intervals=as.integer(FALSE),NAOK=FALSE,PACKAGE="prodlim")
      out <- data.frame(n.risk=xxx$pred.nrisk,n.event=xxx$pred.nevent,n.lost=xxx$pred.nlost)
      lagxxx <- .C("life_table",pred.nrisk=integer(Ntimes*Nstrata),pred.nevent=integer(Ntimes*Nstrata),pred.nlost=integer(Ntimes*Nstrata),nrisk=as.integer(object$n.risk[,1]),nevent=as.integer(object$n.event[,1]),nlost=as.integer(object$n.lost[,1]),as.double(lagTimes),as.double(object$time),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(Nstrata),as.integer(Ntimes),intervals=as.integer(TRUE),NAOK=FALSE,PACKAGE="prodlim")
      out$n.risk <- lagxxx$pred.nrisk
      for (cv in 1:length(object$clustervar)){
        yyy <- .C("life_table",pred.nrisk=integer(Ntimes*Nstrata),pred.nevent=integer(Ntimes*Nstrata),pred.nlost=integer(Ntimes*Nstrata),nrisk=as.integer(object$n.risk[,1+cv]),nevent=as.integer(object$n.event[,1+cv]),nlost=as.integer(object$n.lost[,1+cv]),as.double(times),as.double(object$time),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(Nstrata),as.integer(Ntimes),intervals=as.integer(FALSE),NAOK=FALSE,PACKAGE="prodlim")
        outCV <- data.frame(n.risk=yyy$pred.nrisk,n.event=yyy$pred.nevent,n.lost=yyy$pred.nlost)
        lagyyy <- .C("life_table",pred.nrisk=integer(Ntimes*Nstrata),pred.nevent=integer(Ntimes*Nstrata),pred.nlost=integer(Ntimes*Nstrata),nrisk=as.integer(object$n.risk[,1+cv]),nevent=as.integer(object$n.event[,1+cv]),nlost=as.integer(object$n.lost[,1+cv]),as.double(lagTimes),as.double(object$time),as.integer(object$first.strata[findex]),as.integer(object$size.strata[findex]),as.integer(Nstrata),as.integer(Ntimes),intervals=as.integer(TRUE),NAOK=FALSE,PACKAGE="prodlim")
        outCV$n.risk <- lagyyy$pred.nrisk
        names(outCV) <- paste(object$clustervar,names(outCV))
        out <- cbind(out,outCV)
      }
    }
  }
  if (!is.null(stats)){
    statsList <- lapply(stats,function(x){
      if (percent==TRUE && length(grep(x[1],c("n.event","n.lost","n.risk"),val=FALSE))==0)
        100*as.numeric(c(x[2],object[[x[1]]])[pindex+1])
      else
        as.numeric(c(x[2],object[[x[1]]])[pindex+1])
    })
    names(statsList) <- sapply(stats,function(x)x[[1]])
    add <- do.call("cbind",statsList)
    add <- add[,match(colnames(add),colnames(out),nomatch=FALSE)==0]
    if (NROW(out)==1)
      out <- data.frame(c(out,add))
    else
      out <- cbind(out,add)
  }
  
  #
  # ----------------split according to covariate strata----------------
  #
  if (Nstrata > 1) {
    split.cova <- rep(1:Nstrata,rep(Ntimes,Nstrata))
    out <- split(out,split.cova)
    names(out) <- IndeX$names.strata
    out <- lapply(out,function(x){
      x <- as.matrix(x)
      if (showTime==TRUE){
        if (intervals==TRUE)
          x <- cbind(time0=c(0,round(times[-length(times)],2)),time1=times,x)
        else
          x <- cbind(time=times,x)
        rownames(x) <- 1:NROW(x)
      }
      else{ # times are rownames
        if (intervals==TRUE)
          rownames(x) <- paste("(",paste(c(0,round(times[-length(times)],2)),round(times,2),sep="-"),"]",sep="")
        else
          rownames(x) <- round(times,2)
      }
      x
    })
  }
  else{
    out <- as.matrix(out)
    if (showTime==TRUE){
      if (intervals==TRUE)
        out <- cbind(time0=c(0,round(times[-length(times)],2)),time1=times,out)
      else
        out <- cbind(time=times,out)
      rownames(out) <- 1:NROW(out)
    }
    else{ # times are rownames
      if (intervals==TRUE)
        rownames(out) <- paste("(",paste(c(0,round(times[-length(times)],2)),round(times,2),sep="-"),"]",sep="")
      else
        rownames(out) <- round(times,2)
    }
  }
  out
}