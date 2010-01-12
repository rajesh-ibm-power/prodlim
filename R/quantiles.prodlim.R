"quantile.prodlim" <- function(object,
                               q,
                               newdata,
                               level.chaos=1,
                               mode="list",
                               ...){
  stopifnot(object$model=="survival")
  if (missing(q)) q <- c(1,.75,0.5,.25,0)
  q <- 1-q ## ts says correctly that this is a survival function
  sumx <- summary(object,newdata=object$X,times=object$time,showTime=TRUE,verbose=FALSE)
  getQ <- function(sum){
    out <- do.call("cbind",lapply(c("surv","lower","upper"),function(w){
      nanana=is.na(sum[,w])
      xxx=sum[,w][!nanana]
      ttt=sum[,"time"][!nanana]
      found <- 2+sindex(jump.times=xxx,eval.times=q,comp="greater",strict=FALSE)
      inner <- c(as.vector(c(0,ttt)[found]))
      inner
    }))
    out <- data.frame(out)
    out <- cbind(q,out)
    names(out) <- c("q","quantile","lower","upper")
    out}
  if (is.null(object$X)) getQ(sumx)
  else lapply(sumx,getQ)
}
  
