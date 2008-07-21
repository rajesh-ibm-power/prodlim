predictNevent <- function(object,times,X,intervals=FALSE,stat.alist){
  model <- object$model
  cotype <- object$covariate.type
  p <- predict(object,newdata=X,level.chaos=0,times=times,type="list")
  times <- p$times
  NT <- p$dimensions$time
  NR <- p$dimensions$strata
  pindex <- p$indices$time
  findex <- p$indices$strata
  N <- NROW(object$model.response)
  if(N<=0) stop("Null response")
  Nevent <- as.list(data.frame(object$n.event))
  Nrisk <- as.list(data.frame(object$n.risk))
  dimE <- length(Nevent)
  dimR <- length(Nrisk)
  if (dimR>dimE) stop("Somethings wrong: dim(object$nrisk) > dim(object$nevent))")
  tempE <- lapply(1:dimE,function(e){
    .C("life_table",
       pred.nrisk=integer(NT*NR),
       pred.nevent=integer(NT*NR),
       pred.nlost=integer(NT*NR),
       nrisk=as.integer(Nrisk[[min(e,dimR)]]),
       nevent=as.integer(Nevent[[e]]),
       nlost=as.integer(object$n.lost),
       as.double(times),
       as.double(object$time),
       as.integer(object$first.strata[findex]),
       as.integer(object$size.strata[findex]),
       as.integer(NR),
       as.integer(NT),
       intervals=FALSE,
       NAOK=FALSE,
       PACKAGE="prodlim")
  })
  pred.nevent <- do.call("cbind",lapply(tempE,function(x){x$pred.nevent}))
  pred.nlost <- tempE[[1]]$pred.nlost
  if (intervals!=TRUE)
    pred.nrisk <- do.call("cbind",lapply(tempE[1:dimR],function(x){x$pred.nrisk}))
  else{
    tempR <- lapply(1:dimR,function(r){
      inc <- min(diff(object$time))/10
      Rtimes <- c(min(min(object$time),0)-.1 , times[-length(times)])
      .C("life_table",
         pred.nrisk=integer(NT*NR),
         pred.nevent=integer(NT*NR),
         pred.nlost=integer(NT*NR),
         nrisk=as.integer(Nrisk[[r]]),
         nevent=as.integer(Nevent[[1]]),
         nlost=as.integer(object$n.lost),
         as.double(Rtimes),
         as.double(object$time),
         as.integer(object$first.strata[findex]),
         as.integer(object$size.strata[findex]),
         as.integer(NR),
         as.integer(NT),
         intervals=intervals,
         NAOK=FALSE,
         PACKAGE="prodlim")
    })
    pred.nrisk <- do.call("cbind",lapply(tempR,function(x){x$pred.nrisk}))
  }
  out <- data.frame(n.risk=pred.nrisk,n.lost=pred.nlost,n.event=pred.nevent)
  if (!missing(stat.alist)){
    names <- sapply(stat.alist,function(x)x[1])
    initial <- sapply(stat.alist,function(x)x[2])
    cause <- as.numeric(sapply(stat.alist,function(x)x[3]))
    add <- lapply(1:length(names),function(i){
      oo <- object[[names[i]]]
      if(all(sapply(cause,is.na)==FALSE) && !all(cause <= NCOL(oo))) stop(paste("This competing risk model knows only",NCOL(oo),"causes."),call.=FALSE)
      if (is.null(dim(oo)))
        as.numeric(c(initial[i],oo)[pindex+1])
      else
        as.numeric(c(initial[i],oo[,cause[i]])[pindex+1])
    })
    if (all(!is.na(cause)))
      names(add) <- paste(names,cause,sep=".")
    else
      names(add) <- names
    add <- do.call("cbind",add)
    #    colnames(add) <- names
    out <- cbind(out,add)
  }
  if (NR > 1) {
    split.cova <- rep(1:NR,rep(NT,NR))
    out <- split(out,split.cova)
    names(out) <- p$names.strata
    out <- lapply(out,function(x){
      if (intervals==TRUE)
        rownames(x) <- paste("(",paste(c(0,times[-length(times)]),times,sep="-"),"]",sep="")
      else
        rownames(x) <- round(times,2)
      x
    })
  }
  else{ 
    if (intervals==TRUE)
      rownames(out) <- paste("(",paste(c(0,times[-length(times)]),times,sep="-"),"]",sep="")
    else rownames(out) <- times
  }
  out
}
  
