"Hist" <- function(time,
                   event,
                   id=NULL,
                   cens.code="0") {
  
  # resolving the `time' argument
  # --------------------------------------------------------------------
  if (is.matrix(time)) time <- data.frame(time)
  if (class(time)=="list"){
    if (length(time) !=2 || length(time[[1]])!=length(time[[2]]))
      stop("Argument time has a wrong format")
    time <- data.frame(time)
  }
  if (is.data.frame(time)){
    cens.type <- "intervalCensored"
    L <- time[[1]]
    R <- time[[2]]
    N <- length(L)
    stopifnot(is.numeric(L))
    stopifnot(is.numeric(R))
    stopifnot(all(L<=R || is.na(R)))
    status <- rep(2,N)
    status[L==R] <- 1
    status[is.infinite(R) | is.na(R) | (L!=R & as.character(R)==cens.code)] <- 0
    ## the last part of the condition achieves to things:
    ##     1. for multi-state models allow transitions to a censored state
    ##     2. to ignore this, if an event occured exactly at time 0 and 0 is the cens.code
    R[status==0] <- Inf
  }
  else{
    stopifnot(is.numeric(time))
    cens.type <- "rightCensored"
    N <- length(time)
    status <- rep(1,N) ## temporary dummy
  }
  
  # resolving the argument `event' 
  # --------------------------------------------------------------------

  if (missing(event)){
    model <- "onejump"
    event <-  rep(1,N)
    warning("Argument event is missing:\nassume observations of a survival model\nand only one event per subject")
  }
  else{
    if (is.matrix(event)) event <- data.frame(event)
    if (class(event)=="list"){
      if (length(event) !=2 || length(event[[1]])!=length(event[[2]]))
        stop("Argument event has a wrong format")
      event <- data.frame(event)
    }
    if (!is.data.frame(event)){
      if (is.null(id)){
        model <- "onejump"
        if (is.logical(event)) event <- as.numeric(event)
        status[is.na(event) | is.infinite(event) | as.character(event)==cens.code] <- 0
      }
      else{
        ## inFormat <- "longitudinal"
        stopifnot(is.numeric(id) || is.factor(id))
        model <- "multi.states"
        if (cens.type=="intervalCensored"){
          stop("Dont know the order of transitions for interval censored observations.")
        }
        else{
          # 1. sort the observations by id and time
          sorted <- order(id,time)
          time <- time[sorted]
          ## status <- status[sorted] consists only of 1's
          id <- id[sorted]
          event <- event[sorted]
          # time <- time[duplicated(id)] ## remove the resp. first time
          # status <- status[duplicated(id)]
          last.id <- c(diff(id[sorted]) != 0, 1)
          first.id <- c(1, diff(id[sorted]) != 0)
          from <- event[last.id!=1]
          to <- event[first.id!=1]
          # 2. get back to the original order
          time <- time[sorted]
          id <- id[sorted]
          event <- event[sorted]          
          status[is.na(to) | is.infinite(to) | as.character(to)==cens.code] <- 0
        }
      }
    }
    else{
      ## inFormat <- "from2to"
      model <- "multi.states"
      from <- event[[1]]
      to <- event[[2]]
      status[is.na(to) | is.infinite(to) | as.character(to)==cens.code] <- 0
    }
  }
  
  if (all(status==0)) stop("All observations are censored")
  if (all(status==1)) cens.type <- "uncensored"
  
  if(model=="onejump"){
    
    # 2-state and competing.risks models
    # --------------------------------------------------------------------

    if (is.factor(event)){
      event <- factor(event) # drop unused levels
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

    if (eventISfactor==1) stop("Components of event have different classes")
    
    if (eventISfactor==2)
      states <- unique(c(levels(from),levels(to)))
    else
      states <- as.character(unique(c(from,to)))

    states <- as.character(states[states!=cens.code])
    ## states <- states[match(state.order,states)]
    
    if (cens.code %in% levels(from)) stop("Code for censored data used wrongly in argument event")

    if (cens.code %in% levels(from)) stop("Code for censored data used wrongly in argument event")
    
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
  if (!is.null(id)) history <- cbind(history,id)
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

states <- function(x){
  UseMethod("states")
}

states.Hist <- function(x){
  attr(x,"states")
}

states.prodlim <- function(x){
  attr(x$model.response,"states")
}
