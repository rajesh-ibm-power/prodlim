jackknife.prodlim <- function(x,times){
  ##   xdata <- model.frame(formula=eval(x$call$formula),data=eval(x$call$data))
  xdata <- data.frame(unclass(x$model.response))
  stopifnot(x$type==1)
  stopifnot(x$model=="survival")
  N <- NROW(xdata)
  if (missing(times)) {
    times <- x$time
    S <- x$surv
  }
  else{
    S <- predict(x,times=times,type="surv")
  }
  jackframe <- do.call("rbind",lapply(1:N,function(k){
    ##     print(k)
    xdata.k <- xdata[-k,]
    call.k <- x$call
    call.k$data <- xdata.k
    fitk <- eval(call.k)
    Sk <- predict(fitk,times=times,type="surv")
    Jk <- N * S - (N-1) * Sk
    Jk
  }))
  out <- jackframe
  class(jackframe) <- c("jackknife","prodlim")
  jackframe
}
