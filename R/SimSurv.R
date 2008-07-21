"SimSurv" <- function (N,
                       surv = list(dist="rweibull",
                         args=list(shape=1),
                         baseline=1,
                         link="exp",
                         coef=1,
                         transform=NULL,
                         method="simulation"),
                       cens = list(dist="rexp",
                         args=NULL,
                         baseline=1/100,
                         link="exp",
                         max=NULL,
                         type="right",
                         coef=0,
                         transform=NULL,
                         method="transform"),
                       cova = list(X1=list("rnorm",mean=0,sd=2),
                         X2=list("rbinom",size=1,prob=.5)),
                       verbose=1,
                       ...){
  extra <- match.call()
  extra[[1]] <- as.name("list")
  extra <- extra[match(names(extra),c("N","surv","cens","cova","verbose"),nomatch=0)==0]
  for (e in names(extra)[-1]){# the first `name' is ""
    ename <- unlist(strsplit(e,"\\."))
    switch(ename[1],
           "surv"={
             e2 <- match(ename[2],names(surv),nomatch=0)
             if (e2>0){
               surv[[e2]] <- eval(extra[[e]])
             }
             else
               warning("extra argument ", e , " does not match an argument of surv and is therefore ignored")},
           "cens"={
             e2 <- match(ename[2],names(cens),nomatch=0)
             if (e2>0)
               cens[[e2]] <- eval(extra[[e]])
             else
               warning("extra argument ", e , " does not match an argument of cens and is therefore ignored")},
           warning("extra argument ", e , " is misspecified and therefore ignored"))
  }
  if (NROW(cova)!=N)
    X.matrix <- data.frame(sapply(cova, function(x) {
      do.call(x[[1]], c(n = N, x[-1]))
    }))
  else X.matrix <- cova
  
  NP <- NCOL(X.matrix)

  if (NP != length(surv$coef))
    if (length(surv$coef)==1){
      surv$coef <- rep(surv$coef,NP)
      if (verbose) warning("The survival regression coefficient is used for all covariates.",call.=FALSE)
    }
    else
    stop("Number of covariates and number of survival regression coefficients differ.")
  
  if (!length(surv$baseline)) surv$baseline <- 1

  
  # functional form for linking to cova 
  # --------------------------------------------------------------------

  transform <- list(surv=surv$transform,cens=cens$transform)
  method <- list(surv=surv$method,cens=cens$method)
  
  special.link.surv <- length(transform$surv)>0
  
  if (special.link.surv){
    found.link.surv <- match(names(transform$surv),names(X.matrix),nomatch=FALSE)
    which.link.surv <- match(names(X.matrix),names(transform$surv),nomatch=FALSE)>0
    if (any(!(found.link.surv)))
      stop("Argument `transform$surv' must be a list whose names match the covariate names")
    if (verbose>1){
      print(paste("Survival times linked to covariate(s)",names(transform$surv),"via"))
      print(transform$surv)
    }
    surv.X <- do.call("cbind",lapply(names(X.matrix),function(p){
      if (match(p,names(transform$surv),nomatch=FALSE)){
        sapply(X.matrix[,p],transform$surv[[p]])
      }
      else
        X.matrix[,p]
    }))
    colnames(surv.X) <- names(X.matrix)
    colnames(surv.X)[which.link.surv] <- paste("f",colnames(surv.X)[which.link.surv],sep=".")
  }
  
  # functional form linking cova to censoring  time
  # --------------------------------------------------------------------
  
  special.link.cens <- length(transform$cens)>0
  
  if (special.link.cens){
    if (!length(transform$cens)==NP)
      stop("Wrong number of link functions for censoring model.")
    if (verbose>1){ print("Censoring linked to cova via");print(transform$cens)}
    cens.X <- do.call("cbind",lapply(1:NP,function(p){sapply(X.matrix[,p],transform$cens[[p]])}))
    colnames(cens.X) <- paste("g",names(X.matrix),sep=".")
  }

  # linear predictors
  # --------------------------------------------------------------------
  if (special.link.surv)
    linpred.surv <- colSums(surv$coef * t(surv.X))
  else
    linpred.surv <- colSums(surv$coef * t(X.matrix))

  if (length(cens$coef)>0)
    if (special.link.cens)
      linpred.cens <- colSums(cens$coef * t(cens.X))
    else
      linpred.cens <- colSums(cens$coef * t(X.matrix))
  else linpred.cens <- NULL
  
  surv.method <- match.arg(method$surv,c("simulation","transform"))
  cens.method <- match.arg(method$cens,c("simulation","transform"))
  ##   if (time.varying==TRUE)
  ##     surv.time <- SimSurvInternalTimeVarying(N,surv$dist,surv$args,surv$coef,surv$baseline,X.matrix[,1,drop=TRUE],method,timevar.fun)
  ##   else
  surv.time <- SimSurvInternal(N,surv$dist,surv$args,surv$link,surv$baseline,linpred.surv,surv.method)
  
  if (length(cens$type) > 0 && match(cens$type,c("interval","right"),nomatch=0)) {
    if (length(cens$coef)>0){
      if (NP != length(cens$coef))
        if (length(cens$coef)==1){
          cens$coef <- rep(cens$coef,NP)
          if (verbose) warning("The censoring regression coefficient is used for all covariates.",call.=FALSE)
        }
        else
          stop("Number of covariates and number of censoring regression coefficients differ.")
    }
    if (!length(cens$baseline)) cens$baseline <- 1
    
    ##     print(list(N,cens$dist,cens$args,cens$link,cens$baseline,linpred.cens,cens.method))
    cens.time <- SimSurvInternal(N,cens$dist,cens$args,cens$link,cens$baseline,linpred.cens,cens.method)
    if (is.numeric(cens$max)) cens.time <- pmin(cens.time, cens$max)
    out <- data.frame(cbind(time = pmin(surv.time, cens.time), status = as.numeric(surv.time <= cens.time)),
                      X.matrix)
    if (cens$type=="interval"){
      max <- max(cens.time)
      if (is.null(cens$compliance)) cens$compliance <- .8
      if (is.null(cens$unit)) cens$unit <- round(max/10,2)
      if (is.null(cens$lateness)) cens$lateness <- round(max/30,2)
      icens <- SimSurvInternalIntervalCensored(N=N,unit=cens$unit,lateness=cens$lateness,compliance=cens$compliance,withdraw.time=cens.time,event.time=surv.time)
      out <- cbind(icens,out)
      out$status[out$L<out$R & !is.na(out$R)] <- 2
      
    }
  ##     if (!missing(cluster))
  ##       out$ID <- ID
    
    if (special.link.cens) out <- cbind(out,cens.X)
    if (special.link.surv) out <- cbind(out,surv.X[,which.link.surv,drop=FALSE])
    
    #    if (keep.uncensored)
    out <- cbind(out,uncensored.time=surv.time)
    if (verbose>0) cat(paste(round(100 * sum((1 - out$status>0))/N), "%Censoring"), "\n")
    out <- out[order(out$time,-out$status),]
    attr(out,"formula") <- formula(paste("Hist(time,status)~", paste(names(X.matrix), collapse = "+")))
  }
  else{
    out <- data.frame(cbind(time = surv.time, status=rep(1,N), X.matrix))
    if (special.link.cens) out <- cbind(out,cens.X)
    if (special.link.surv) out <- cbind(out,surv.X[,which.link.surv,drop=FALSE])
  }
#  if (special.link.cens) attr(out,"cova.link.cens") <- transform$cens
#  if (special.link.surv) attr(out,"cova.link.surv") <- transform$surv
  attr(out,"call") <- match.call()
  class(out) <- c("SimSurv",class(out))
  row.names(out) <- rep(1:N)
  out
}


SimSurvInternal <- function(N,dist, args, link, baseline, linpred, method){
  if (method=="transform"){
    if (length(linpred)>0)
      Stime <- sapply(do.call(link, list(x = linpred)), function(lp) {
        do.call(dist, c(list(n = 1, lp), args)) })
    else{
      Stime <- do.call(dist, c(list(n=N),args))
    }
  }
  else{  ## via simulation lemma
    Stime <- switch(dist,
                    "rexp"= (1/baseline) * (-log(runif(N)) * exp(-linpred)),
                    "rweibull"= - (log(runif(N)) * (1 / baseline) * exp(-linpred))^args$shape,
                    "rgompertz"=(1/args$alpha) * log(-(args$alpha/baseline) * (log(runif(N)) * exp(-linpred)) + 1))
  }
}

SimSurvInternalTimeVarying <- function(N,dist, args, coef, baseline, x, method,fun){
  switch(fun,{
    Stime <- -log(1-runif(N))/baseline
    Stime[x!=0] <-  (-baseline - sqrt(baseline^2-4*x[x!=0]*log(1-runif(sum(x!=0)))))/(2*x[x!=0])
    Stime[Stime<0] <-  (-baseline + sqrt(baseline^2-4*x[Stime<0]*log(1-runif(sum(Stime<0)))))/(2*x[Stime<0])
    Stime
  },{
    Stime <- -log(1-runif(N))/baseline
    Stime[x!=0] <-  (exp(-coef*log(1-runif(N[x!=0]))/x[x!=0])-1)/coef
    Stime
  })
  
}

SimSurvInternalIntervalCensored <- function(N,unit,lateness,compliance,withdraw.time,event.time){
  interval <- do.call("rbind",lapply(1:N,function(i){
    fixed <- seq(0,withdraw.time[i],unit)
    stopifnot(unit>lateness)
    v <- rnorm(length(fixed),0,lateness)
    prob <- rbinom(length(fixed),1,compliance)
    grid <- sort((fixed+v)[prob==1])
    found <- sum(grid<event.time[i])
    if (is.na(found)||found==0){
      L <- 0
      R <- grid[1]
    }
    else{
      L <- grid[found]
      R <- grid[found+1]
    }
    if (is.na(R)) L <- max(L,withdraw.time[i])
    ##     print(c("L"=L,"R"=R,event.time[i]))
    if ((!is.na(L+R))&& (L>event.time[i]||R<event.time[i]))
      print(cbind(g=grid,L,R))
    c("L"=L,"R"=R)
  }))
  interval <- data.frame(interval)
  interval$L[interval$L<0] <- 0
  interval  
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
