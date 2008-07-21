"plot.prodlim" <- function(x,
                           what,
                           cause=1,
                           newdata,
                           add = FALSE,
                           col,
                           lty,
                           lwd,
                           ylim,
                           xlim,
                           xlab="Time",
                           ylab,
                           legend=TRUE,
                           legend.args=NULL,
                           mark.time=FALSE,
                           cex.mark.time=1.5,
                           conf.int=TRUE,
                           conf.int.args=list(col="gray"),
                           atrisk=ifelse(add,FALSE,TRUE),
                           atrisk.args,
                           timeOrigin=0,
                           ...){

  # 1. extracting the values 
  # --------------------------------------------------------------------
  if (x$cens.type=="intervalCensored") conf.int <- FALSE
  if (x$cens.type=="intervalCensored") atrisk <- FALSE
  
  model <- x$model
  
  if (missing(what))
    what <- switch(model,
                   "survival"="survival",
                   "competing.risks"="incidence",
                   "multi.states"="hazard")
  else
    what <- match.arg(what,
                      c("survival","incidence","hazard"))

  clusterp <- !is.null(x$clustervar)
  
  if (x$cens.type=="intervalCensored")
    plot.times <- sort(unique(x$time[2,]))
  else{
    plot.times <- sort(unique(x$time))
    if (plot.times[1]>timeOrigin) plot.times <- c(timeOrigin,plot.times)
  }
  if (x$covariate.type>1){
    if (missing(newdata) || length(newdata)==0){
      newdata <- x$X
      ## names(newdata) <- sapply(strsplit(names(newdata),"strata.\|NN."),function(x)x[2])
      if (NROW(newdata)>10){
        ## warning(paste("\nModel is fitted at ",NROW(newdata)," different covariate values ... \nonly the fits at the first, the median and the last entry of fit$X  are shown. see argument `newdata' for other behaviors.\n\n"))
        newdata <- newdata[c(1,round(median(1:NROW(newdata))),NROW(newdata)),,drop=FALSE]
      }
    }
    if (what=="survival" || model=="survival"){
      Y <- predict(x,times = plot.times,newdata = newdata,level.chaos=1,type="surv")
    }
    else{
      if (length(cause)!=1) stop("With covariates, only cumulative incidence of one cause can be plotted. make use of the  argument `add'.")
      temp <- predictCuminc(x, newdata = newdata, times = plot.times, mode="matrix",cause=cause,level.chaos=1)
      Y <- split(temp[[1]],1:NROW(temp[[1]]))
      ## print(Y)
      names(Y) <- rownames(temp[[1]])
      # if (!("title" %in% names(legend.args))) legend.args <- c(title=paste("Cause",cause),legend.args)
      # print(Y)
      # Y <- unlist(lapply(temp,function(x){split(x,1:NROW(x))}),recursive=FALSE)
    }
  }
  else{
    if (model=="survival" || what=="survival")
      Y <- list(predict(x,times=plot.times,level.chaos=0,type="surv"))
    else{
      Y <- predict(x,times=plot.times,cause=cause,type="cuminc")
    }
  }

  if (model=="survival" && what=="incidence") Y <- lapply(Y,function(x){1-x})
  if (!is.null(x$reverse) && x$reverse==TRUE) Y <- lapply(Y,function(x){1-x})
  
  # 2. setting graphical parameters and generating empty plot
  # --------------------------------------------------------------------

  nlines <- length(Y)
  
  if (missing(ylab)) ylab <- ifelse(what=="survival",ifelse(x$reverse==TRUE,"Censoring probability","Survival probability"),"Cumulative incidence")
  if (missing(xlab)) xlab <- "Time"
  if (missing(lwd)) lwd <- 3
  
  if (atrisk==TRUE && !add){
    oldmar <- par()$mar
    if (missing(atrisk.args)) atrisk.args <- list()
    if (!("interspace" %in% names(atrisk.args))) atrisk.args$interspace <- 1
    if (!("dist" %in% names(atrisk.args))) atrisk.args$dist <- .3
    if (!length(atrisk.args$automar) || atrisk.args$automar==TRUE){
      newmar <- par()$mar + ## user level mar
        c(par()$mgp[2] + ## distance of xlab from xaxis
          atrisk.args$dist+ ## distance of atrisk numbers from xlab
          ifelse(clusterp,2,1)*nlines ## number of atrisk lines
          + 1, ## one extra line below the bottom number atrisk line
          0,
          0,
          0)
      par(mar=newmar)
    }
  }
  if (!add) {
    if (missing(xlim)) xlim <- c(0, max(plot.times))
    if (missing(ylim)) ylim <- c(0, 1)
    plot(0, 0, type = "n", ylim = ylim, xlim = xlim, xlab = xlab, ylab = ylab, ...)
  }
  if (missing(col)) col <- 1:nlines
  if (missing(lty)) lty <- rep(1, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  if (atrisk==TRUE) par(mar=oldmar) ## reset
  
  # 3. adding the lines 
  # --------------------------------------------------------------------
  
  nix <- lapply(1:nlines, function(s) {
    lines(x = plot.times, y = Y[[s]], type = "s", col = col[s], lty = lty[s], lwd = lwd)
  })
  
  if (mark.time==TRUE && model=="survival"){
    mark.times <- sort(unique(x$time[x$n.lost>0]))
    mark.pos <- match(mark.times,plot.times)
    nix <- lapply(1:nlines, function(s) {
      points(x = mark.times,
             y = Y[[s]][mark.pos],
             pch="I",
             type = "p",
             col = col[s],
             lty = lty[s],
             lwd = lwd,
             cex=cex.mark.time)
    })}
  
    if (conf.int==TRUE){
      if (!("col" %in% names(conf.int.args))) conf.int.args$col <- "gray55"
      if (!("times" %in% names(conf.int.args)))
        at <- ceiling(quantile(1:length(x$time),c(.25,.5,.75,1)))
      else
        at <- sindex(x$time,conf.int.args$times)
      ##     print(x$time)
      ##     print(conf.int.args$times)
      ##     print(cbind(conf.int.args$times[at>0],x$lower[at],x$upper[at]))
      apply(cbind(conf.int.args$times[at>0],x$lower[at],x$upper[at]),1,function(ci){
        segments(x0=ci[1],x1=ci[1],y0=ci[2],y1=ci[3],lwd=2,col=conf.int.args$col,lty=1)
        #      points(ci[1],ci[2],pch="-",col=col,lty=1)
        #      points(ci[1],ci[3],pch="-",col=col,lty=1)
      })
    }

  # adding the no. of individuals at risk
  # --------------------------------------------------------------------
  
  if (missing(newdata)) newdata <- NULL

  if (atrisk==TRUE && !add){
    do.call("atRisk",c(list(x=x),list(newdata=newdata),atrisk.args))
  }

  # legend
  # --------------------------------------------------------------------
  if(legend==TRUE && !add && !is.null(names(Y))){
    save.xpd <- par()$xpd
    par(xpd=TRUE)
    if (!("legend" %in% names(legend.args))) legend.args <- c(list(legend=names(Y)),legend.args)
    if (!("lwd" %in% names(legend.args))) legend.args <- c(list(lwd=lwd),legend.args)
    if (!("col" %in% names(legend.args))) legend.args <- c(list(col=col[1:nlines]),legend.args)
    if (!("lty" %in% names(legend.args))) legend.args <- c(list(lty=lty[1:nlines]),legend.args)
    ##     if (!("cex" %in% names(legend.args))) legend.args <- c(list(cex=1.5),legend.args)
    if (!("bty" %in% names(legend.args))) legend.args <- c(list(bty="n"),legend.args)
    if (!("y.intersp" %in% names(legend.args))) legend.args <- c(list(y.intersp=1.3),legend.args)
    if (!("x" %in% names(legend.args))) legend.args <- c(list(x=0),legend.args)
    if (!("y" %in% names(legend.args))) legend.args <- c(list(y=ifelse(what=="survival",ylim[1]+0.3*diff(ylim),ylim[2])),legend.args)
    do.call("legend",legend.args)
    par(xpd=save.xpd)
  }
  invisible(x)
}


atRisk <- function(x,
                   formula,
                   data,
                   curves,
                   newdata,
                   times,
                   line,
                   interspace,
                   cex,
                   labels,
                   pos,
                   adj,
                   dist,
                   adjust.labels=TRUE){

  if (missing(times)) times <- seq(0,x$maxtime,x$maxtime/10)

  if (!missing(formula) & !missing(data))
    x <- prodlim(formula,data)
  else
    if (class(x)!="prodlim")
      stop("Either formula and data or x must be given")

  px <- predictNevent(object=x,times=times,X=newdata)

  if (is.matrix(px) || is.data.frame(px))
    sumx <- lapply(px[,grep("n.risk",colnames(px)),drop=FALSE],function(x)x)
  else
    sumx <- lapply(px,function(v){
      lapply(v[,grep("n.risk",colnames(v))],function(x){x})})

  nlines <- length(sumx)
  if (missing(line)) line <- par()$mgp[2] + dist + (0:(2*nlines-1)) *interspace
  if (missing(cex)) cex <- 1
  if (missing(pos)) pos <- min(times)
  if (missing(adj)) adj <- 1.5

  if (missing(labels))
    if (length(names(sumx)==nlines))
      labels <- paste("[",names(sumx),"]",sep="")
    else
      labels <- c("No.   \nat-risk",rep("",nlines-1))
  
  # labeling the no. at-risk below plot
  # --------------------------------------------------------------------
  if (is.null(adjust.labels) || adjust.labels==TRUE){
    labels <- format(labels,justify="left")}
  
  lapply(1:nlines,function(y){
    mtext(text=as.character(sumx[[y]]),side=1,at=times,line=rep(line[y],length(times)),cex=cex,outer=FALSE,xpd=NA)
    mtext(text=labels[y],side=1,at=pos,line=line[y],adj=adj,cex=cex,outer=FALSE,xpd=NA)
  })
}

