plot.Hist <- function(x,
                      layout,
                      xbox.rule=.3,
                      ybox.rule=1.5,
                      state.lab,
                      state.cex=2,
                      rect.args,
                      args.state.lab,
                      arrow.lab,
                      arrow.lab.offset,
                      arrow.lab.side=NULL,
                      arrow.lab.cex=2,
                      arrow.lab.style,
                      arrows.args,
                      arrow.lab.args=NULL,
                      arrow.head.offset=3,
                      arrow.double.dist=1,
                      arrow.fix.code=NULL,
                      enumerate.boxes=FALSE,
                      box.numbers,
                      cex.boxlabs=1.28,
                      margin,
                      verbose=FALSE,
                      reverse,
                      ...){
  # symetric frame
  # --------------------------------------------------------------------
  
  oldmar <- par()$mar
  oldxpd <- par()$xpd
  if (!missing(margin))
    par(mar=margin,xpd=TRUE)
  else
    par(mar=c(2,2,2,2),xpd=TRUE)
  
  # find states and labs that appear inside the boxes
  # --------------------------------------------------------------------
  
  model.type <- attr(x,"model")
  states <- attr(x,"states")
  NS <- length(states)
  if (all(as.character(as.numeric(as.factor(states)))==states))
    states <- paste(switch(model.type,"survival"="Event","competing.risks"="Cause","State"),states)
  ##   print(states)
  
  x <- x[x[,"status"]!=attr(x,"cens.code"),]
  if (NROW(x)==0) stop("No uncensored transition.")
  
  if (missing(reverse)) reverse <- FALSE
  ##     reverse <- ifelse(model.type=="competing.risks",TRUE,FALSE)
  
  if (!missing(state.lab)){
    if (model.type!="multi.states"){## maybe add an initial state
      if (length(state.lab)==NS){
        states <- c("Initial",if(reverse) rev(state.lab) else state.lab)
        NS <- length(states)
      }
      else
        if(length(state.lab)-1==NS){
          states <- state.lab
          NS <- length(states)
          if (reverse==TRUE) states <- c(states[1],rev(states[-1]))
        }
        else stop("Wrong number of state names.")
    }
    else{
      if (length(state.lab)==NS){
        states <- state.lab
        if (reverse==TRUE) states <- rev(states) 
      }
      else
        stop("Wrong number of state names.")
    }
  }
  else{
    if (model.type!="multi.states"){## add an initial state
      states <- c("Initial",if(reverse) rev(states) else states)
      NS <- length(states)
    }
  }
  
  # find transitions between the states
  # --------------------------------------------------------------------
  
  if (model.type=="multi.states"){  
    
    ## `from' and `to' originate from
    ## factors with levels cens.code,states
    ## and are now integer valued 1 for censored
    ## and 2,..., NS otherwise
    ## thus we have to start counting from 2

    if (is.null(transitions <- summary(x)$trans.frame)){
      transitions <- data.frame(x[,c("from","to")]) 
      ##     print(unique(transitions))
      transitions$from <- factor(transitions$from)
      levels(transitions$from) <- states[as.numeric(levels(transitions$from))-1]
      transitions$to <- factor(transitions$to)
      levels(transitions$to) <- states[as.numeric(levels(transitions$to))-1]
      ordered.transitions <- unique(transitions)
      ##     print(ordered.transitions)
      ##     stop()
    }
    else{
      ordered.transitions <- transitions
    }
  }
  else{
    if (model.type=="survival")
      transitions <- table(x[,"status"])
    else
      transitions <- table(x[,"event"][x[,"status"]!=0])
    ordered.transitions <- data.frame(cbind(rep(states[1],length(states[-1])),states[-1]))
  }
  names(ordered.transitions) <- c("from","to")
  
  
  # arranging the boxes
  # --------------------------------------------------------------------
  
  if (missing(layout)){
    if (model.type=="multi.states"){
      sumx <- summary(x,verbose=FALSE)
      state.types <- unlist(sumx$states)
      state.types <- state.types[state.types>0]
      if (NS==3 && all(unlist(state.types)==1)){ # illness-death-model
        auto.col <- rep(FALSE,3)
        box.row <- c(3,1,3)
        box.col <- c(1,2,3)
        nrow <- 3
        ncol <- 3
      }
      else{
        ncol <- length(state.types)
        nrow <- max(state.types)
        auto.col <- rep(TRUE,ncol)
        box.col <- rep(1:ncol,state.types)
        box.row <- unlist(sapply(state.types,function(x)1:x))
        ## print(box.row)
      }
    }
    else{
      state.types <- list(initial=1,absorbing=NS-1)
      ncol <- 2
      nrow <- max(unlist(state.types))
      auto.col <- c(TRUE,TRUE)
      box.row <- c(1,1:state.types$absorbing)
      box.col <- c(1,rep(2,state.types$absorbing))
    }
  }
  else{
    nrow <- layout$nrow
    ncol <- layout$ncol
    box.row <- sapply(layout$box.pos,function(x)x[1])
    box.col <- sapply(layout$box.pos,function(x)x[2])
    #    box.col <- layout$box.col
    auto.col <- layout$auto.col
    if (is.null(auto.col))
      auto.col <- rep(FALSE,ncol)
  }
  
  ##   print(ordered.transitions)
  
  # plot an empty frame
  # --------------------------------------------------------------------
  
  Xlim <- 100
  Ylim <- 100
  plot(0,0,type="n",xlim=c(0,Xlim),ylim=c(0,Ylim),xlab="",ylab="",axes=FALSE)
  
  # fixed size boxes
  # --------------------------------------------------------------------

  state.width <- sapply(states,strwidth,cex=state.cex)
  max.width <- max(state.width)
  state.height <- sapply(states,strheight,cex=state.cex)
  max.height <- max(state.height)
  box.width <- max.width + xbox.rule * max.width
  
  if ((ncol * box.width) > Xlim) stop("The horizontal dimension of the boxes is too big -- change parameters `state.cex' and/or `xbox.rule'.")
  
  box.height <- max.height + ybox.rule * max.height
  
  if ((nrow * box.height) > Ylim)
    stop("The vertical dimension of the boxes is too big -- change parameters `state.cex' and/or `ybox.rule'.")
  
  # position of boxes
  # --------------------------------------------------------------------

  ## distribute the boxes uniformly 
  position.finder <- function(border,len,n){
    if (n==1)
      (border - len)/2
    else{
      seq(0,border-.5*len,len + (border-(n * len))/(n-1))
    }
  }
  
  ##   ybox.position <- unlist(lapply(state.types,position.finder,border=Ylim,len=box.height))
  ##   ybox.position <- unlist(lapply(1:ncol,position.finder,border=Ylim,len=box.height))
  ##   ybox.position <- rep(position.finder(Ylim,box.height,nrow),NS)

  
  Ypossible <- lapply(1:ncol,function(gc){
    if(auto.col[gc]>0)
      rev(position.finder(Ylim,box.height,sum(box.col==gc)))
    else rev(position.finder(Ylim,box.height,nrow))
  })
  
  ybox.position <- sapply(1:NS,function(s){
    Ypossible[[box.col[s]]][box.row[s]]
  })
  
  
  ##   if (NS==3 && length(grep("transient",names(state.types)))>0){ # illness-death-model
  ##     ybox.position <- c(0,Ylim-box.height,0)
  ##   }
  Xpossible <- position.finder(border=Xlim,len=box.width,n=ncol)
  xbox.position <- sapply(1:NS,function(s){
    Xpossible[box.col[s]]
  })
  ##   print(xbox.position)  
  ##   xbox.position <- rep(position.finder(border=Xlim,len=box.width,n=ncol),state.types)
  ##   print(ybox.position)
  ##   print(xbox.position)
  
  
  # draw the boxes
  # --------------------------------------------------------------------

  rect.args.defaults <- list(xpd=TRUE)
  if (missing(rect.args)) rect.args <- c(list(xleft=xbox.position,ybottom=ybox.position,xright=xbox.position+box.width,ytop=ybox.position+box.height),rect.args.defaults)
  else   rect.args <- c(list(xleft=xbox.position,ybottom=ybox.position,xright=xbox.position+box.width,ytop=ybox.position+box.height),rect.args,rect.args.defaults)
  rect.args <- rect.args[!duplicated(names(rect.args))]
  do.call("rect",rect.args)
  
  # center the  text in boxes
  # --------------------------------------------------------------------

  xtext.position <- xbox.position + (box.width - state.width)/2
  ytext.position <- ybox.position + (box.height - state.height)/2
  
  if (missing(args.state.lab)) {
    args.state.lab <- list(x=xtext.position,y=ytext.position,labels=states,adj=c(0,0),cex=state.cex)
  }
  else {
    args.state.lab <- c(list(x=xtext.position,y=ytext.position,labels=states,adj=c(0,0),cex=state.cex),args.state.lab)
    args.state.lab <- args.state.lab[!duplicated(args.state.lab)]
  }
  do.call("text",args.state.lab)
  

  # maybe put numbers in the upper left corner of the boxes
  # --------------------------------------------------------------------
  
  if (enumerate.boxes==TRUE){
    if (missing(box.numbers)) {
      if (reverse)
        box.numbers <- c(0,rev(1:(length(xbox.position)-1)))
      else
        box.numbers <- 0:(length(xbox.position)-1)
    }
    nix <- lapply(1:length(xbox.position),function(x) {
      lab <- box.numbers[x]
      text(x=xbox.position[x],
           y=ybox.position[x]+box.height,
           labels=lab,
           cex=cex.boxlabs,
           adj=c(-.5,1.43))})
  }
  
  if (verbose==TRUE) {
    print(data.frame(cbind(xbox.position,ybox.position,states)))
  }
  
  # find default arg values and labs for the arrows
  # --------------------------------------------------------------------
  ## print(arrows.args.defaults)

  N <- NROW(ordered.transitions)
  doubleArrow <- match(paste(ordered.transitions[,"to"],ordered.transitions[,"from"]),
                     paste(ordered.transitions[,"from"],ordered.transitions[,"to"]),nomatch=0)
  ordered.transitions <- cbind(ordered.transitions,doubleArrow)
  ##   print(ordered.transitions)
  if (missing(arrows.args)) arrows.args <- list()
  arrows.args.defaults <- list(lwd=2)
  
  if (!missing(arrows.args) && length(arrows.args)>0){
    arrows.args <- c(arrows.args,arrows.args.defaults)
    arrows.args <- arrows.args[!duplicated(names(arrows.args))]
  }
  else arrows.args <- arrows.args.defaults

  ## Added Fri Jun 20 20:01:56 CEST 2008
  if (!missing(arrow.lab) && arrow.lab==FALSE){
    arrow.lab.style <- "character"
    arrow.lab <- rep("",N)
  }
  ## Added Fri Jun 20 20:01:56 CEST 2008
  
    if (missing(arrow.lab.style))
      arrow.lab.style <- if (missing(arrow.lab)) "symbolic" else "character"
  arrow.lab.style <- match.arg(arrow.lab.style,c("symbolic","character","count","none"))
  if (arrow.lab.style=="character")
    if (missing(arrow.lab)) arrow.lab <- rep("",N)
    else
      if (length(arrow.lab)!=N) stop(paste("Number of arrow labels and the number of arrows dont match. Provide",N,"labels"))
  if (arrow.lab.style=="count"){
    ##     arrow.lab <- paste("n=",as.numeric(transitions),sep="")
    arrow.lab <- as.numeric(transitions)
  }
  if (arrow.lab.style!="none"){
    arrow.lab <- lapply(1:N,function(trans){
      from.state <- match(ordered.transitions[trans,1],states)
      to.state <- match(ordered.transitions[trans,2],states)
      if (arrow.lab.style=="symbolic")
        lab <- bquote(alpha[.(paste(from.state-1,to.state-1,sep=""))](t))
      else
        lab <- arrow.lab[trans]
      lab
    })
    if (reverse) arrow.lab <- rev(arrow.lab)
  }

  ## ------------------------------------------------------------- ##
  ## the default offset is computed as the string width and height ##
  ## of the labels                                                 ##
  ## ------------------------------------------------------------- ##
  if (!missing(arrow.lab.offset)){
    if (length(arrow.lab.offset)==1)
      arrow.lab.offset <- rep(arrow.lab.offset,N)
    else
      if(length(arrow.lab.offset)!=N)
        stop(paste("Length of arrow.lab.offset unequal to the no. of arrows which is",N))
  }
  else{
    arrow.lab.offset <- lapply(1:N,function(trans){
      c(strwidth(arrow.lab[[trans]],cex=arrow.lab.cex),
        strheight(arrow.lab[[trans]],cex=arrow.lab.cex))
    })
  }
  
  # draw the arrows
  # --------------------------------------------------------------------

  arrow.pos <- lapply(1:N,function(trans){
    from.state <- match(ordered.transitions[trans,1],states)
    to.state <- match(ordered.transitions[trans,2],states)
    
    ##     print(c(from.state,to.state))
    ## print(list(Box1=c(round(xbox.position[from.state],4),round(ybox.position[from.state],4)),Box2=c(round(xbox.position[to.state],4),round(ybox.position[to.state],4)),BoxDim=c(box.width,box.height),offset=arrow.head.offset,verbose=verbose))
    Apos <- findArrow(Box1=c(round(xbox.position[from.state],4),round(ybox.position[from.state],4)),Box2=c(round(xbox.position[to.state],4),round(ybox.position[to.state],4)),BoxDim=c(box.width,box.height),offset=arrow.head.offset,verbose=verbose)
    dir <- attr(Apos,"direction")
    names(Apos) <- c("x0","y0","x1","y1")
    if (Apos[1]==Apos[3]){ # vertical arrow
      x1 <- Apos[1]
      x2 <- Apos[3]
      if (Apos[2]<Apos[4]) {
        y1 <- Apos[2] + arrow.head.offset
        y2 <- Apos[4] - arrow.head.offset
      }
      else{
        y1 <- Apos[2] - arrow.head.offset
        y2 <- Apos[4] + arrow.head.offset
      }
    }
    else{
      thisform <- function(x,y){
        (y[4]-y[2])*x/(y[3]-y[1])-(y[4]-y[2])*y[1]/(y[3]-y[1])+y[2]
      }
      x1 <- Apos[1] + arrow.head.offset
      x2 <- Apos[3] - arrow.head.offset
      y1 <- thisform(x=x1,y=Apos)
      y2 <- thisform(x=x2,y=Apos)
    }
    Apos <- c(x1,y1,x2,y2)
    names(Apos) <- c("x0","y0","x1","y1")
    list(Apos=Apos,ADirection=dir)
  })
  arrow.dir <- sapply(arrow.pos,function(x)x$ADirection)
  arrow.pos <- data.frame(do.call("rbind",lapply(arrow.pos,function(x)x$Apos)))
  if (verbose==TRUE) {
    cat("\n")
    print(cbind(ordered.transitions,"Arrow.pos"=arrow.pos))
    cat("\n")
  }

  ordered.transitions <- cbind(ordered.transitions,arrow.pos)
  
  if (missing(arrow.double.dist)) arrow.double.dist <- NULL

  # draw the arrows and the labels
  # --------------------------------------------------------------------
  drawarrows <- lapply(1:N,function(trans){
    acode <- if (arrow.dir[trans]==2 && ordered.transitions[trans,"doubleArrow"]==0) 1 else 2
    dd <- ordered.transitions[trans,"doubleArrow"]

    # modify the arrow coordinates in case of  <==>
    # --------------------------------------------------------------------
    if (dd!=0)
      if (is.null(arrow.double.dist))
        dist <- strwidth("a",cex=arrow.lab.cex)
      else dist <- arrow.double.dist
    else
      dist <- 0
    if (dd==2)
      if (dist==0)
        acode <- 3
      else
        acode <- 1
    apos <- ordered.transitions[trans,c("x0","y0","x1","y1"),drop=FALSE]
    if (dist>0)
      if (abs(apos$x1-apos$x0) < abs(apos$y1-apos$y0)){
        if (dd==2) {
          if ((apos$x0+apos$x1)<=Xlim){ ## left side of the figure
            apos <- apos + c(dist,0,dist,0)
          }
          else{
            apos <- apos - c(dist,0,dist,0)
          }
        }
        else{
          if (dd==1){
            if ((apos$x0+apos$x1)<=Xlim){ ## left side of the figure
              apos <- apos - c(dist,0,dist,0)
            }
            else{
              apos <- apos + c(dist,0,dist,0)
            }
          }
        }
      }
      else{
        if (dd==2){
          apos <- apos + c(0,dist,0,dist)
        }
        else{ if (dd==1){
          apos <- apos - c(0,dist,0,dist)
        }}}
    
    ##     print(c(apos,c(arrows.args,code=acode)))
    if (!is.null(arrow.fix.code))
      acode <- arrow.fix.code[trans]
    ##       if (length(arrow.lab)!=N) stop(paste("Number of arrow labels and the number of arrows dont match. Provide",N,"labels"))
    do.call("arrows",c(apos,c(arrows.args,code=acode)))
    off <- arrow.lab.offset[[trans]]
    slope <- function(x,y,a,b){(y-b)/(x-a)}
    S <- slope(apos$x0,apos$y0,apos$x1,apos$y1)
    if (is.na(S)) S <- Inf
    amid <- c(unlist(apos[c(1,2)],use.names=FALSE)+unlist(apos[c(3,4)],use.names=FALSE))/2
    
    # place the lab
    # --------------------------------------------------------------------
    ## 
    ##       q1     |     q2
    ##              |
    ##     ---------|-----------
    ##              |
    ##       q4     |     q3
    
    q1 <- c(-1,1)
    q2 <- c(1,1)
    q3 <- c(1,-1)
    q4 <- c(-1,-1)
    
    ##     print(acode)
    ##     print(S)
    Q <- if (amid[1]<=Xlim/2 && amid[2]<Ylim/2) 4
    else if (amid[1]<=Xlim/2 && amid[2]>=Ylim/2) 1
    else if (amid[1]>Xlim/2 && amid[2]<Ylim/2) 3
    else 2 ## (amid[1]>Xlim/2 && amid[2]>=Ylim/2)
    
    ## print(paste("Q=",Q))
    ##     corners <- cbind(xbox.position,ybox.position)
    ##     D <- apply(corners,1,function(box){slope(amid[1],amid[2],box[1],box[2])})
    if (acode==2){
      if (is.infinite(S)) sign <- switch(Q,q1,q2,q2,q1)
      else if (S==0) sign <- switch(Q,q1,q2,q3,q4)
      else if (NS==2) if (S<0) sign <- q2 else sign <- q1
      else if (Q==1) if (S<0) sign <- q2 else sign <- q3
      else if (Q==2) if (S<0) sign <- q4 else sign <- q1
      else if (Q==3) if (S<0) sign <- q4 else sign <- q1
      else if (Q==4) if (S<0) sign <- q4 else sign <- q1
    }
    else{
      if (is.infinite(S)) sign <- switch(Q,q2,q1,q1,q2)
      else if (S==0) sign <- switch(Q,q3,q4,q1,q2)
      else if (NS==2) if (S<0) sign <- q1 else sign <- q3
      else if (Q==1) if (S<0) sign <- q2 else sign <- q1
      else if (Q==2) if (S<0) sign <- q2 else sign <- q3
      else if (Q==3) if (S<0) sign <- q2 else sign <- q3
      else if (Q==4) if (S<0) sign <- q2 else sign <- q3 
    }
    if (is.null(arrow.lab.side) || is.na(arrow.lab.side[trans]))
      as <- 1 else as <- arrow.lab.side[trans]
    ## if (acode==1) sign <- -1*sign
    sign <- as*sign
    ##          abline(h=Xlim/2)
    ##          abline(v=Ylim/2)
    ##     print(S)
    ##     points(x=amid[1],y=amid[2])
    ##     off <- off
    if (S==0){
      ## sign <- c(ifelse(amid[1]<Xlim/2,-1,1),ifelse(amid[2]<Ylim/2,-1,1))
      if (verbose) print("horizontal")
      xy <- amid + sign*c(0,1)*(off-dist)
    }
    else{
      ## sign <- c(ifelse(amid[1]<Xlim/2,-1,1),ifelse(amid[2]<Ylim/2,-1,1))
      if (is.infinite(S)){
        if (verbose) print("vertical")
        xy <- amid + sign*c(1,0)*(off-dist)
      }
      else{
        if (S<0){
          S <- -S
          if (S>1){
            if (verbose) print("down steep")
            ##             if (acode==2)
            ##               sign <- c(1,1)
            ##             else
            ##               sign <- c(-1,-1)
            xy <- amid + sign*c(1/S, 1)*(off-dist)
          }
          else{
            if (verbose) print("down not steep")
            ##             if(match(ordered.transitions[trans,"from"],states,nomatch=0)==1)
            ##               sign <- c(-1,-1)
            ##             else
            ##               sign <- c(1,1)
            xy <- amid + sign*c(S,1)*(off-dist)
          }
        }
        else{
          if (S>0)
            if (S>1){
              if (verbose) print("up steep")
              ##               if (acode==2)
              ##                 sign <- c(-1,-1)
              ##               else
              ##                 sign <- c(1,1)
              xy <- amid + sign*c(1/S,1)*(off-dist)
            }
            else{
              if (verbose) print("up not steep")
              ##               if (acode==2) sign <- c(1,-1)
              ##               else sign <- c(-1,1)
              xy <- amid + sign*c(S,1)*(off-dist)
            }
        }
      }
    }
    ## print(xy)
    if (is.null(arrow.lab.args)) arrow.lab.args <- list(x=xy[1],
                                                        y=xy[2],
                                                        labels=bquote(arrow.lab[[trans]]),
                                                        cex=arrow.lab.cex)
    else{
      arrow.lab.args <- c(arrow.lab.args,list(x=xy[1],
                                              y=xy[2],
                                              labels=bquote(arrow.lab[[trans]])))
      arrow.lab.args <- arrow.lab.args[!duplicated(arrow.lab.args)]
    }
    do.call("text",arrow.lab.args)
  })
  par(mar=oldmar,xpd=oldxpd)
  invisible(x)
}



