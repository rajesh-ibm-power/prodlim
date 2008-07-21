summary.Hist <- function(object, verbose=TRUE,...){
  
  D <- object[,"status",drop=TRUE]
  states <- attr(object,"states")
  cens.code <- attr(object,"cens.code")
  
  # resolving events and model states 
  # --------------------------------------------------------------------
  model <- attr(object,"model")
  model.string <- paste("response of a", model,"model")
  if (model=="survival")
    E <- object[,"status",drop=TRUE]
  else
    if (model=="competing.risks")
      E <- object[,"event",drop=TRUE]

  if (model=="multi.states"){
    from <- object[,"from"]
    to <- object[,"to"]
    stupid.stupid.factor.levels <- as.integer(factor(c(cens.code,states),levels=c(cens.code,states)))
    sorted.stupid.stupid.factor.levels <- c(stupid.stupid.factor.levels[-1],stupid.stupid.factor.levels[1])
    code.to <- factor(to,levels=sorted.stupid.stupid.factor.levels,
                      labels=c(states,"unknown"))
    code.from <- factor(from,levels=sorted.stupid.stupid.factor.levels,
                        labels=c(states,"unknown"))
    trans.frame <- data.frame(from=code.from,to=code.to)
    Transitions <- apply(cbind(as.character(code.from),as.character(code.to)),1,paste,collapse=" -> ")
    obnoxious.factor.levels <- unique(Transitions)
    Transitions <- factor(Transitions,obnoxious.factor.levels)
    transitions <- table(Transitions)
    state.types <- factor(as.numeric(match(states,unique(code.from),nomatch=0)!=0) + 2*as.numeric(match(states,unique(code.to),nomatch=0)!=0),
                          levels=c(1,2,3))
    names(state.types) <- states
    levels(state.types) <- c("initial","absorbing","transient")
    summary.out <- list(states=table(state.types),transitions=transitions,trans.frame=trans.frame)
    if (verbose==TRUE){
      state.table <- as.matrix(transitions)
      colnames(state.table) <- c("Freq")
    }
  }
  else{
    summary.out <- NULL
  }
  # resolving the censoring mechanism
  # --------------------------------------------------------------------
  if (verbose==TRUE){
    cens.type <- attr(object,"cens.type")
    capitalize <- function(x) {
      s <- strsplit(x, " ")[[1]]
      paste(toupper(substring(s, 1,1)), substring(s, 2), sep="", collapse=" ")
    }
    
    cens.string <- capitalize(cens.type)
    Observations <- switch(cens.type,
                           "intervalCensored"=factor(D,levels=c(1,2,0),labels=c("exact.time","intervalCensored","rightCensored")),
                           "rightCensored"=factor(D,levels=c(1,0),labels=c("event","right.censored")),
                           "uncensored"=factor(D,labels=c("event")))

    Freq <- table(Observations)

    cat(paste("\n",cens.string," ",model.string,"\n",sep=""))

    cat("\nNo.Observations:",NROW(object),"\n\nPattern:\n")
    
    switch(model,"survival"={
      prmatrix(cbind(names(Freq),Freq),quote=FALSE,rowlab=rep("",NROW(Freq)))},
           "competing.risks"={
             if (cens.type=="rightCensored"){
               stupid.stupid.factor.levels <- as.integer(factor(c(cens.code,states),levels=c(cens.code,states)))
               sorted.stupid.stupid.factor.levels <- c(stupid.stupid.factor.levels[-1],stupid.stupid.factor.levels[1])
               events <- factor(object[,"event"],
                                levels=sorted.stupid.stupid.factor.levels,
                                labels=c(states,"unknown"))
             }
             else{
               stupid.stupid.factor.levels <- as.integer(factor(c(cens.code,states),levels=c(cens.code,states)))
               sorted.stupid.stupid.factor.levels <- c(stupid.stupid.factor.levels[-1],stupid.stupid.factor.levels[1])
               events <- factor(object[,"event"],
                                levels=sorted.stupid.stupid.factor.levels,
                                labels=c(states,"unknown"))
             }
             prout <- table("Cause"=events,as.character(Observations))
             print(prout)
           },
           "multi.states"={
             ##              cat("Observations:\n")
             print(table(Transitions,Observations))
             ##              print(cbind(Transition=transitions))
             ##              print(table(transitions,Observations))
           })
    ## prmatrix(transitions,quote=FALSE,rowlab=rep("",NROW(transitions)))})
  }
  if (!is.null(summary.out))
    invisible(summary.out)
  else
    invisible(NULL)
}
