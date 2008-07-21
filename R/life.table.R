life.table <- function(object,
                       times,
                       newdata,
                       max.tables=20,
                       verbose=1){
  if (object$covariate.type>1){
    if (missing(newdata) || length(newdata)==0){
      X <- object$X
      if (NROW(X)>max.tables){
        if (verbose>0) warning(call.=TRUE,immediate.=TRUE,paste("\nPredicted survival probabilities for",NROW(X),"covariate values available.\nShown are only the first, the median and the last table ...\nto see other tables use arguments `X' and `max.tables'\n"))
        X <- X[c(1,round(median(1:NROW(X))),NROW(X)),,drop=FALSE]
      }
    }
    else{
      X <- newdata
    }
  }
  else X <-  NULL
  times <- times[times<=max(object$time)]
  x <- predict(object,newdata=X,level.chaos=0,times=times,type="list")
  times <- x$times
  NT <- x$dimensions$time
  NR <- x$dimensions$strata
  pindex <- x$indices$time
  findex <- x$indices$strata
  ltab <- .C("life_table",
             pred.nrisk=integer(NT*NR),
             pred.nevent=integer(NT*NR),
             pred.nlost=integer(NT*NR),
             nrisk=as.integer(object$n.risk),
             nevent=as.integer(object$n.event),
             nlost=as.integer(object$n.lost),
             as.double(times),
             as.double(object$time),
             as.integer(object$first.strata[findex]),
             as.integer(object$size.strata[findex]),
             as.integer(NR),
             as.integer(NT),
             NAOK=FALSE,
             PACKAGE="prodlim")
  iid <- is.null(object$clustervar)
  if (!iid){  
    ltab.cluster <- .C("life_table",
                       pred.nrisk=integer(NT*NR),
                       pred.nevent=integer(NT*NR),
                       pred.nlost=integer(NT*NR),
                       nrisk=as.integer(object$n.risk[,2]),
                       nevent=as.integer(object$n.event[,2]),
                       nlost=as.integer(object$n.lost),
                       as.double(times),
                       as.double(object$time),
                       as.integer(object$first.strata[findex]),
                       as.integer(object$size.strata[findex]),
                       as.integer(NR),
                       as.integer(NT),
                       NAOK=FALSE,
                       PACKAGE="prodlim")
    if (length(object$conf.int)>0)
      out <- split(data.frame(time=rep(times,NR),
                              nrisk=ltab$pred.nrisk,
                              nrisk.cluster=ltab.cluster$pred.nrisk,
                              nevent=ltab$pred.nevent,
                              ncluster.with.event=ltab.cluster$pred.nevent,
                              nlost=ltab$pred.nlost,
                              hazard=c(0,object$hazard)[pindex+1],
                              surv=c(1,object$surv)[pindex+1],
                              se.surv=c(0,object$se.surv)[pindex+1],
                              lower=c(0,object$lower)[pindex+1],
                              upper=c(0,object$upper)[pindex+1]),rep(1:NR,rep(NT,NR)))
    else
      out <- split(data.frame(time=rep(times,NR),
                              nrisk=ltab$pred.nrisk,
                              nrisk.cluster=ltab.cluster$pred.nrisk,
                              nevent=ltab$pred.nevent,
                              ncluster.with.event=ltab.cluster$pred.nevent,
                              nlost=ltab$pred.nlost,
                              hazard=c(0,object$hazard)[pindex+1],
                              surv=c(1,object$surv)[pindex+1],
                              se.surv=c(0,object$se.surv)[pindex+1]),rep(1:NR,rep(NT,NR)))
  }
  else{
    if (object$model=="competing.risks"){
      NE <- NCOL(object$cuminc)
      pred.nevent <- ltab$pred.nevent
      for (j in 2:NE){
        ltab.cr <- .C("life_table",
                      pred.nrisk=integer(NT*NR),
                      pred.nevent=integer(NT*NR),
                      pred.nlost=integer(NT*NR),
                      nrisk=as.integer(object$n.risk),
                      nevent=as.integer(object$n.event[,j]),
                      nlost=as.integer(object$n.lost),
                      as.double(times),
                      as.double(object$time),
                      as.integer(object$first.strata[findex]),
                      as.integer(object$size.strata[findex]),
                      as.integer(NR),
                      as.integer(NT),
                      NAOK=FALSE,
                      PACKAGE="prodlim")
        pred.nevent <- cbind(pred.nevent,ltab.cr$pred.nevent)
      }
      colnames(pred.nevent) <- as.character(1:NE)
      if (object$conf.int>0){
        out <- lapply(1:NE,function(e){
          inner <- data.frame(time=rep(times,NR),
                              nrisk=ltab$pred.nrisk,
                              nevent=pred.nevent[,e],
                              surv=c(1,object$surv)[pindex+1],
                              cuminc=rbind(rep(0,NE),object$cuminc)[pindex+1,e],
                              se.cuminc=rbind(rep(0,NE),object$se.cuminc)[pindex+1,e],
                              lower=rbind(rep(0,NE),object$lower)[pindex+1,e],
                              upper=rbind(rep(0,NE),object$upper)[pindex+1,e])
          inner})
        names(out) <- paste("event",1:NE)
      }
      else{
        out <- data.frame(time=rep(times,NR),
                          nrisk=ltab$pred.nrisk,
                          nevent=pred.nevent,
                          surv=c(1,object$surv)[pindex+1],
                          cuminc=rbind(rep(0,NE),object$cuminc)[pindex+1,],
                          se.cuminc=rbind(rep(0,NE),object$se.cuminc)[pindex+1,])
      }
    }
    else{
      if(NR==1){
        out <- data.frame(time=times,
                          nrisk=ltab$pred.nrisk,
                          nevent=ltab$pred.nevent,
                          nlost=ltab$pred.nlost,
                          hazard=c(0,object$hazard)[pindex+1],
                          surv=c(1,object$surv)[pindex+1],
                          se.surv=c(0,object$se.surv)[pindex+1])
        if (length(object$conf.int)>0)
          out <- cbind(out,data.frame(lower=c(0,object$lower)[pindex+1],upper=c(0,object$upper)[pindex+1]))
      }
      else{
        x.hazard <- matrix(c(0,object$hazard)[pindex+1],nrow=NR,ncol=NT,byrow=TRUE)
        x.surv <- matrix(c(1,object$surv)[pindex+1],nrow=NR,ncol=NT,byrow=TRUE)
        x.se.surv <- matrix(c(0,object$se.surv)[pindex+1],nrow=NR,ncol=NT,byrow=TRUE)
        x.nrisk <- matrix(ltab$pred.nrisk,nrow=NR,ncol=NT,byrow=TRUE)
        x.nlost <- matrix(ltab$pred.nlost,nrow=NR,ncol=NT,byrow=TRUE)
        x.nevent <- matrix(ltab$pred.nevent,nrow=NR,ncol=NT,byrow=TRUE)
        if (length(object$conf.int)>0){
          x.lower <- matrix(c(0,object$lower)[pindex+1],nrow=NR,ncol=NT,byrow=TRUE)
          x.upper <- matrix(c(0,object$upper)[pindex+1],nrow=NR,ncol=NT,byrow=TRUE)
        }
        out <- lapply(1:NR,function(r){
          inner <- data.frame(time=times,
                              nrisk=x.nrisk[r,],
                              nevent=x.nevent[r,],
                              nlost=x.nlost[r,],
                              hazard=x.hazard[r,],
                              surv=x.surv[r,],
                              se.surv=x.se.surv[r,])
          if (length(object$conf.int)>0){
            inner <- cbind(inner,data.frame(lower=x.lower[r,],upper=x.upper[r,]))
          }
          inner
        })
        names(out) <- x$names.strata
      }
    }
  }
  out
}
  
