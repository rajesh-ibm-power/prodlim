SimSurv <- function(N,
                    surv,
                    cens,
                    cova,
                    verbose=1,
                    ...){
  
  default.surv.args <- list(dist="rweibull",args=list(shape=1),baseline=1/100,link="exp",coef=c(1,-1),transform=NULL)
  default.cens.args <- list(dist="rexp",args=NULL,baseline=1/100,link="exp",max=NULL,type="right",coef=NULL,transform=NULL)
  if (missing(cova)) cova <- list(X1=list("rnorm",mean=0,sd=2),X2=list("rbinom",size=1,prob=.5))
  smartA <- resolveSmartArgs(call=match.call(),
                             keys=c("surv","cens","cova"),
                             ignore=c("n","N","surv","cens","cova","verbose"),
                             defaults=list("surv"=default.surv.args,"cens"=default.cens.args),
                             verbose=TRUE)
  surv <- smartA$surv
  
  censnotwanted <- (!missing(cens) && is.logical(cens) && cens==FALSE)
  
  if (!censnotwanted) cens <- smartA$cens

  # ------------------------resolving covariates------------------------
  
  X.matrix <- do.call(resolveX,c(list(N=N),object=list(cova)))
  NP <- NCOL(X.matrix)

  # -----------special links between survival and covariates-----------
  
  if (length(surv$transform)>0){
    survSpecials <- TRUE
    surv.X <- transformX(X=X.matrix,transform=surv$transform,transName="f")
  }
  else{
    survSpecials <- FALSE
    surv.X <- X.matrix
  }
  
  # ------------------resolving the linear predictors------------------

  linpred.surv <- resolveLinPred(X=surv.X,coef=surv$coef,verbose=verbose)

  # ---------------------------survival time---------------------------
  
  surv.time <- SimSurvInternal(N=N,dist=surv$dist,args=surv$args,link=surv$link,baseline=surv$baseline,linpred=linpred.surv)
  
  # ---------------------------censoring----------------------------

  #  if (!censnotwanted && length(cens$type)>0)
  #    censType <- match(cens$type,c("interval","right"),nomatch=0)

  if (censnotwanted==TRUE)
    cens.time <- rep(Inf,N)
  else{

    # --------special links between right censoring and covariates-----------
    
    if (length(cens$transform)>0){
      censSpecials <- TRUE
      cens.X <- transformX(X=X.matrix,transform=cens$transform,transName="f")
    }
    else{
      censSpecials <- FALSE
      cens.X <- X.matrix
    }

    linpred.cens <- resolveLinPred(X=cens.X,coef=cens$coef,verbose=verbose)
    
    cens.time <- SimSurvInternal(N=N,dist=cens$dist,args=cens$args,link=cens$link,baseline=cens$baseline,linpred=linpred.cens)
    
    if (is.numeric(cens$max)) cens.time <- pmin(cens.time, cens$max)

    
    if (!censnotwanted && cens$type=="interval"){
      if (is.null(cens$compliance)) cens$compliance <- .95
      if (is.null(cens$unit)) cens$unit <- round(median(cens.time)/3,2)
      if (is.null(cens$lateness)) cens$lateness <- round(median(cens.time)/30,2)
      icens <- SimSurvInternalIntervalCensored(N=N,
                                               unit=cens$unit,
                                               lateness=cens$lateness,
                                               compliance=cens$compliance,
                                               withdraw.time=cens.time,
                                               event.time=surv.time)
    }
  }

  # ------the censoring status: 0= right, 1= observed, 2= interval------

  status <- as.numeric(surv.time <= cens.time)

  # ------------------------preparing the output------------------------
  
  out <- data.frame(cbind(time = pmin(surv.time, cens.time),
                          status = status,
                          uncensored.time=surv.time,
                          cens.time=cens.time))
  if (!is.null(X.matrix)) out <- cbind(out,X.matrix)
  if (survSpecials) out <- cbind(out,surv.X[0<match(names(surv.X),names(out),nomatch=0)])
  if (!censnotwanted){
    if (censSpecials) out <- cbind(out,cens.X[0<match(names(cens.X),names(out),nomatch=0)])
    if (cens$type=="interval"){
      out <- cbind(icens,out)
      if (cens$type=="interval"){
        out$status[is.infinite(out$R)] <- 0
        out$cens.time <- pmin(out$cens.time,out$L)
        out$status[!(is.infinite(out$R))] <- 2
        out$status[out$L==out$R] <- 1
      }
    }}

  # --------------------------sorting the data--------------------------
  
  out <- out[order(out$time,-out$status),]
  

  
  # -----------------reporting the censoring percentage-----------------

  if (verbose>0){
    cat(paste(round(100 * sum((1 - out$status>0))/N), "% right censoring"), "\n")
  }
  
  # -----------------------adding some attributes-----------------------
  
  attr(out,"formula") <- formula("Hist(time,status)~1")
  attr(out,"call") <- match.call()
  class(out) <- c("SimSurv",class(out))
  row.names(out) <- rep(1:N)
  out
}


SimSurvInternal <- function(N, dist, args, link, baseline, linpred){
  if (length(linpred)==0) linpred <- 0
  Stime <- switch(dist,
                  "runif"= do.call("runif",c(list(n=N),args)),
                  "rexp"= (1/baseline) * (-log(runif(N)) * exp(-linpred)),
                  "rweibull"= (- (log(runif(N)) * (1 / baseline) * exp(-linpred)))^(1/args$shape),
                  "rgompertz"=(1/args$alpha) * log(1 - (args$alpha/baseline) * (log(runif(N)) * exp(-linpred))))
  Stime
}

SimSurvInternalTimeVarying <- function(N,dist, args, coef, baseline, x, fun){
  #  idea from yanqing sun:
  # lambda = 1/sqrt(t) * exp(sqrt(t)) 
  #
  #  switch(fun,{
  #    Stime <- -log(1-runif(N))/baseline
  #    Stime[x!=0] <-  (-baseline - sqrt(baseline^2-4*x[x!=0]*log(1-runif(sum(x!=0)))))/(2*x[x!=0])
  #    Stime[Stime<0] <-  (-baseline + sqrt(baseline^2-4*x[Stime<0]*log(1-runif(sum(Stime<0)))))/(2*x[Stime<0])
  #    Stime
  #  },{
  #    Stime <- -log(1-runif(N))/baseline
  #    Stime[x!=0] <-  (exp(-coef*log(1-runif(N[x!=0]))/x[x!=0])-1)/coef
  #    Stime
  #  })
}

SimSurvInternalIntervalCensored <- function(N,
                                            unit,
                                            lateness,
                                            compliance,
                                            withdraw.time,
                                            event.time){
  Intervals <- do.call("rbind",lapply(1:N,function(i){
    schedule <- seq(0,withdraw.time[i],unit)
    M <- length(schedule)
    g <- c(0,rep(unit,M))
    # introduce normal variation of the visit times
    g <- g+c(abs(rnorm(1,0,lateness)),rnorm(M,0,lateness))
    grid <- c(0,cumsum(g))
    # remove visits after the end of follow-up time
    grid <- grid[grid<withdraw.time[i]]
    # remove intermediate visits
    if (compliance<1){
      stopifnot(compliance>0)
      missed <- rbinom(length(grid),1,compliance)==0
      grid <- grid[missed==FALSE]
    }
    if (length(grid)==0){
      L <- 0
      R <- Inf
    }
    else{
      posTime <- sindex(jump.times=grid,
                        eval.times=event.time[i])
      L <- grid[posTime]
      R <- grid[posTime+1]
      if (is.na(R)){
        R <- Inf
      }
    }
    c(L=L,R=R)
  }))
  out <- data.frame(Intervals)
  out
}
  
find.baseline <- function(x=.5,
                          setting,
                          verbose=FALSE){
  N <- setting$N
  f <- function(y){
    setting$cens.baseline <- y
    ncens <- sum(do.call("SimSurv",replace(setting,"verbose",verbose))$status==0)
    x-ncens/N
  }
  base.cens <- uniroot(f,c(exp(-50),1000000),tol=.0000001,maxiter=100)$root
  new.setting <- setting
  new.setting$cens.baseline <- base.cens
  do.call("SimSurv",replace(new.setting,"verbose",TRUE))
  new.setting
}

quantile.SimSurv <- function(x,B=10,na.rm=FALSE,probs=.9){
  callx <- attr(x,"call")
  nix <- do.call("rbind",lapply(1:B,function(b){
    quantile(eval(callx)$time,probs)
  }))
  nix <- colMeans(nix)
  nix
}

