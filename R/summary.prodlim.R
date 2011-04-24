# {{{ header
"summary.prodlim" <- function(object,
                              times,
                              newdata,
                              verbose=TRUE,
                              max.tables=20,
                              surv=TRUE,
                              cause,
                              intervals=FALSE,
                              percent=FALSE,
                              showTime=TRUE,
                              ...) {
  # }}}
  # {{{  classify the situation
  cens.type <- object$cens.type         # uncensored, right or interval censored
  model <- object$model                 # survival, competing risks or multi-state
  ## cluster <- object$clustervar          # clustered data?
  cotype <- object$covariate.type       # no, discrete, continuous or both
  # }}}
  # {{{  times
  jump.times <- object$time
  if (missing(times) && (length(times <- jump.times) > 50)) 
    times <- quantile(sort(unique(jump.times)))
  times <- sort(unique(times))
  if (any(times>max(jump.times)))
    warning(call.=TRUE,immediate.=TRUE,paste("\n","Time(s) ",paste(times[times>max(jump.times)],collapse=", ")," are beyond the maximal follow-up time ",max(jump.times),"\n"))
  ntimes <- length(times)
  if (cens.type=="intervalCensored"){
    temp <- data.frame(time=paste("(",paste(signif(object$time[1,],2),
                         signif(object$time[2,],2),
                         sep="-"),"]",sep=""),
                       n.risk=signif(object$n.risk,2),
                       n.event=signif(object$n.event,2),
                       ##    n.lost=object$n.lost,
                       surv=object$surv)
  }
  else{
    # }}}
# {{{  covariates
    if (cotype>1){
      if (missing(newdata) || length(newdata)==0){
        X <- object$X
        if (NROW(X)>max.tables){
          if (verbose>0) warning(call.=TRUE,immediate.=TRUE,paste("\nLife tables are available for",NROW(X),"different covariate constellations.\nShown are only the first, the median and the last table ...\nto see other tables use arguments `X' and `max.tables'\n"))
          X <- X[c(1,round(median(1:NROW(X))),NROW(X)),,drop=FALSE]
        }
      }
      else X <- newdata
    }
    else
      X <- NULL
    
    if (model=="survival") {
      stats <- list(c("surv",1),c("se.surv",0))
      if (!is.null(object$conf.int))
        stats <- c(stats,list(c("lower",0),c("upper",1)))
      if (surv==FALSE){
        object$cuminc <- 1-object$surv
        object$se.cuminc <- object$se.surv
        cuminc.upper <- 1-object$lower
        cuminc.lower <- 1-object$upper
        object$lower <- cuminc.lower
        object$upper <- cuminc.upper
        stats <- list(c("cuminc",0),c("se.cuminc",0))
        if (!is.null(object$conf.int))
          stats <- c(stats,list(c("lower",0),c("upper",1)))
      }
    }
    if (model=="competing.risks"){
      stats <- list(c("cuminc",0),c("se.cuminc",0))
      if (!is.null(object$conf.int))
        stats <- c(stats,list(c("lower",0),c("upper",0)))
    }
    temp <- lifeTab(object=object,times=times,newdata=X,stats=stats,intervals=intervals,percent=percent,showTime=showTime)
    if (model=="competing.risks" & !missing(cause)){
      if (is.numeric(cause) && !is.numeric(names(temp)))
        cause <- attr(object$model.response,"states")[cause]
        # found <- match(cause,attr(object$model.response,"states"),nomatch=FALSE))
        Found <- match(cause,names(temp))
      if (all(Found)) temp <- temp[Found]
      else warning("could not find requested causes in attributes of object$mode.response")
    }
  }
  if (verbose==TRUE){
    if (model=="survival")
      if (cotype==1)
        print(temp,quote=FALSE,...)
      else
        print.listof(temp,quote=FALSE,...)
    else
      if (model=="competing.risks")
        for (cc in 1:length(temp)){
          cat("\n\n----------> Cause: ",names(temp)[cc],"\n\n")
          if (cotype==1)
            print(temp[[cc]],quote=FALSE,...)
          else
            print.listof(temp[[cc]],quote=FALSE,...)
        }
  }
  # }}}
invisible(temp)
}
