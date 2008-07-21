"summary.prodlim" <- function(object,
                              times,
                              newdata,
                              verbose=TRUE,
                              max.tables=20,
                              incidence=FALSE,
                              cause=1,
                              intervals=FALSE,
                              proz=FALSE,
                              ...) {
  # model characteristics
  # --------------------------------------------------------------------
  cens.type <- object$cens.type
  model <- object$model
  cluster <- object$clustervar
  cotype <- object$covariate.type
  
  # times
  # --------------------------------------------------------------------
  jump.times <- object$time
  if (missing(times) && (length(times <- jump.times) > 50)) 
    times <- quantile(sort(unique(jump.times)))
  times <- sort(unique(times))
  ntimes <- length(times)
  
  if (cens.type=="intervalCensored"){
    print(data.frame(time=paste("(",paste(object$time[1,],object$time[2,],sep="-"),"]",sep=""),
                     n.risk=round(object$n.risk,2),
                     n.event=round(object$n.event,2),
                     ##    n.lost=object$n.lost,
                     surv=object$surv))
  }
  else{
    
    # covariates
    # --------------------------------------------------------------------
    if (cotype>1){
      if (missing(newdata) || length(newdata)==0){
        X <- object$X
        ##       names(X) <- sapply(strsplit(names(X),"strata.\|NN."),function(x)x[2])
        if (NROW(X)>max.tables){
          if (verbose>0) warning(call.=TRUE,immediate.=TRUE,paste("\nPredicted survival probabilities for",NROW(X),"covariate values available.\nShown are only the first, the median and the last table ...\nto see other tables use arguments `X' and `max.tables'\n"))
          X <- X[c(1,round(median(1:NROW(X))),NROW(X)),,drop=FALSE]
        }
      }
      else X <- newdata
    }
    else
      X <- NULL

    if (model=="survival") {
      stat.alist <- list(c("surv",1),c("se.surv",0))
      if (!is.null(object$conf.int))
        stat.alist <- c(stat.alist,list(c("lower",0),c("upper",1)))
      if (incidence==TRUE){
        object$incidence <- 1-object$surv
        object$se.inc <- object$se.surv
        object$inc.upper <- 1-object$lower
        object$inc.lower <- 1-object$upper
        stat.alist <- list(c("incidence",0),c("se.inc",0))
        if (!is.null(object$conf.int))
          stat.alist <- c(stat.alist,list(c("inc.lower",0),c("inc.upper",1)))
      }
    }
    if (model=="competing.risks"){
      stat.alist <- list(c("cuminc",1,cause),c("se.cuminc",0,cause))
      if (!is.null(object$conf.int))
        stat.alist <- c(stat.alist,list(c("lower",0,cause),c("upper",1,cause)))
    }
    temp <- predictNevent(object,times,X,stat.alist=stat.alist,intervals=intervals)
    
    if (verbose==TRUE)
      if (cotype>1){
        if (proz==TRUE){
          for (i in 1:length(temp)){
            for (n in sapply(stat.alist,function(x)x[1]))
              temp[[i]][,n] <- 100* temp[[i]][,n]
          }
        }
        print.listof(temp,quote=FALSE,...)
      }
      else{
        if(proz==TRUE){
          for (n in sapply(stat.alist,function(x)x[1]))
            temp[,n] <- 100* temp[,n]
        }
        print(temp,quote=FALSE,...)
      }
    invisible(temp)
  }
}
