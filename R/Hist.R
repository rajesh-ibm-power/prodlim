"Hist" <- function(time,
                   event,
                   id=NULL,
                   cens.code="0") {
  
  # resolving the `time' argument
  # --------------------------------------------------------------------
  
  NT <- NROW(time)
  
  if (is.list(time) && NT==2){
    cens.type <- "intervalCensored"
    L <- time[[1]]
    R <- time[[2]]
    stopifnot(is.numeric(L))
    stopifnot(is.numeric(R))
  }
  else{
    if (NCOL(time)==2){
      cens.type <- "intervalCensored"
      L <- time[,1,drop=TRUE]
      R <- time[,2,drop=TRUE]
      stopifnot(is.numeric(L))
      stopifnot(is.numeric(R))
    }

    else{
      N <- NT
      stopifnot(is.numeric(time))
      cens.type <- "rightCensored"
      status <- rep(1,NT) ## only temporary
    }
  }
  
  if (cens.type=="intervalCensored"){
    N <- length(L)
    stopifnot(length(L)==length(R))
    stopifnot(all(L<=R || is.na(R)))
    status <- rep(2,N)
    status[L==R] <- 1
    status[is.infinite(R) | is.na(R) | as.character(R)==cens.code] <- 0
    R[status==0] <- Inf
  }
  
  # resolving argument `event'
  # --------------------------------------------------------------------

  if (missing(event)){
    model <- "onejump"
    if (cens.type=="intervalCensored"){
      event <-  rep(1,N)
    }
    else {
      event <-  rep(1,N)
      warning("argument `event' is missing:\nassume uncensored observations of a survival model\nwith only one event per subject")
    }
  }
  else{
    NE <- NROW(event)
    if (!is.null(id)){
      if (class(event) %in% c("list","data.frame"))
        stop("When id is given event must be a vector\nshowing the transitions of a multi states model in longitudinal form.")
      model <- "multi.states"
      if (cens.type=="intervalCensored"){
        sorted <- order(id,L)
        L <- L[sorted]
        R <- R[sorted]
        L <- L[duplicated(id)] ## remove the resp. first time
        R <- R[duplicated(id)] ## remove the resp. first time 
        status <- status[sorted]
      }
      else{
        sorted <- order(id,time)
        time <- time[sorted]
        time <- time[duplicated(id)] ## remove the resp. first time
        status <- status[duplicated(id)]
        ## status <- status[sorted] consists only of 1's
      }
      last.id <- c(diff(id[sorted]) != 0, 1)
      first.id <- c(1, diff(id[sorted]) != 0)
      from <- event[sorted][last.id!=1]
      to <- event[sorted][first.id!=1]
      status[is.na(to) | is.infinite(to) | as.character(to)==cens.code] <- 0
    }
    else{
      if (is.list(event) && NE==2){
        model <- "multi.states"
        from <- event[[1]]
        to <- event[[2]]
        status[is.na(to) | is.infinite(to) | as.character(to)==cens.code] <- 0
      }
      else{
        if (NCOL(event)==2){
          model <- "multi.states"
          from <- event[,1,drop=TRUE]
          to <- event[,2,drop=TRUE]
          status[is.na(to) | is.infinite(to) | as.character(to)==cens.code] <- 0
        }
        else{
          if (is.logical(event)) event <- as.numeric(event)
          status[is.na(event) | is.infinite(event) | as.character(event)==cens.code] <- 0
          model <- "onejump"
        }
      }
    }
  }
  
  if (all(status==0)) stop("All observations are censored")
  if (all(status==1)) cens.type <- "uncensored"
  
  if(model=="onejump"){
    
    # 2-state and competing.risks models
    # --------------------------------------------------------------------

    if (is.factor(event)){
      states <- levels(event)
      ## states <- states[match(state.order,states)]
    }
    else
      states <- sort(as.character(unique(event)))
    
    states <- as.character(states[states!=cens.code])
    
    if (length(states)>1)
      model <- "competing.risks"
    else
      model <- "survival"
    
    if (cens.type=="intervalCensored"){
      if (model=="survival")
        history <- cbind(L=L,
                         R=R,
                         status=status)
      else
        history <- cbind(L=L,
                         R=R,
                         status=status,
                         event=as.integer(factor(event,levels=c(cens.code,states))))
    }
    else{
      if (model=="survival")
        history <- cbind(time=time,status=status)
      else{
        history <- cbind(time=time,
                         status=status,
                         event=as.integer(factor(event,levels=c(cens.code,states))))
      }
    }
  }
  else{
    
    # multi.state models
    # --------------------------------------------------------------------
    
    if (any(as.character(from)==as.character(to))) stop("Data contain transitions from state x to state x")
    
    eventISfactor <- as.numeric(is.factor(from)) + as.numeric(is.factor(to))

    if (eventISfactor==1) stop("Components of `event' have different classes")
    
    if (eventISfactor==2)
      states <- unique(c(levels(from),levels(to)))
    else
      states <- as.character(unique(c(from,to)))

    states <- as.character(states[states!=cens.code])
    ## states <- states[match(state.order,states)]
    
    if (cens.code %in% levels(from)) stop("Code for censored data used wrongly in argument `event'")

    if (cens.code %in% levels(from)) stop("Code for censored data used wrongly in argument `event'")
    
    if (cens.type=="intervalCensored"){
      history <- cbind(L=L,
                       R=R,
                       status=status,
                       from=as.integer(factor(from,levels=c(cens.code,states))),
                       to=as.integer(factor(to,levels=c(cens.code,states))))
    }
    else{
      # print(levels(factor(from,levels=c(cens.code,states))))
      history <- cbind(time=time,
                       status=status,
                       from=as.integer(factor(from,levels=c(cens.code,states))),
                       to=as.integer(factor(to,levels=c(cens.code,states))))
    }
  }
  class(history) <- "Hist"
  attr(history,"states") <- states
  attr(history,"cens.type") <- cens.type
  attr(history,"cens.code") <- as.character(cens.code)
  attr(history,"model") <- model
  history
}


"[.Hist" <- function(x,i,j,drop=FALSE) {
  # If only 1 subscript is given, the result will still be a Hist object
  # If the second is given extract the relevant columns as a matrix
  if (missing(j)) {
    y <- x
    y <- unclass(y)
    y <- y[i,, drop=FALSE]
    attr(y,"class") <- attr(x,"class")
    ## attr(y,"dimnames") <- attr(x,"dimnames")
    attr(y,"states") <- attr(x,"states")
    attr(y,"model") <- attr(x,"model")
    attr(y,"cens.type") <- attr(x,"cens.type")
    attr(y,"cens.code") <- attr(x,"cens.code")
    y
  }
  else {
    class(x) <- NULL
    NextMethod("[")
  }
}


is.na.Hist <- function(x) {
  as.vector( (1* is.na(unclass(x)))%*% rep(1, ncol(x)) >0)
}

str.Hist <- function(x){
  class(x) <- "matrix"
  str(x)
}

