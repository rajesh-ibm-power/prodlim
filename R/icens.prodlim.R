icens.prodlim <- function(L,
                          R,
                          grid,
                          status=NULL,
                          ml=FALSE,
                          maxiter=1000,
                          tol=7){
  ntol <- 10^{-tol}
  N <- length(L)
  sorted <- order(L,R)
  R <- R[sorted]
  L <- L[sorted]
  Rna <- (is.infinite(R) | is.na(R))
  status <- status[sorted]
  
  ##   dyn.load("/home/bs/tag/R/packages/prodlim/src/icens_prodlim.so")
  
  if (ml==FALSE) {
    if (length(status==N)){
      status <- as.numeric(status>0)
    }
    else
      status <- as.numeric(!(Rna))
    R[Rna] <- L[Rna]
    if (missing(grid)) grid <- sort(unique(c(L,R)))
    indexR <- sindex(grid,R)
    indexL <- sindex(grid,L)
    NS <- length(grid)
    ##  maxtime <- max(c(L[!is.infinite(L)],R[!is.infinite(R)]))
    Ind <- iindex(L,R,grid)
    fit <- .C("icens_prodlim",as.double(L),as.double(R),as.double(grid),as.integer(indexL),as.integer(indexR),as.integer(Ind$iindex),as.integer(c(Ind$imax,0)),as.integer(status),as.double(N),as.double(NS),nrisk=double(NS),nevent=double(NS),ncens=double(NS),hazard=double(NS),varhazard=double(NS),surv=double(NS),oldsurv=double(NS),as.double(ntol),as.integer(maxiter),n.iter=integer(1),package="prodim")
    ##   fit <- .C("old_icens_prodlim",as.double(L),as.double(R),as.double(grid),as.integer(status),as.double(N),as.double(NS),nrisk=double(NS),nevent=double(NS),ncens=double(NS),hazard=double(NS),varhazard=double(NS),surv=double(NS))
    out <- list("call"=NULL,"formula"=NULL,"time"=c(0,grid),"n.risk"=c(N,round(pmax(0,fit$nrisk),tol)),"n.event"=c(0,round(pmax(0,fit$nevent),tol)),"n.lost"=c(0,round(fit$ncens,tol)),"hazard"=c(0,round(fit$hazard,tol)),"surv"=c(1,round(pmax(0,fit$surv),tol)),"model.response"=NULL,"model"="survival","cens.type"="intervalCensored","covariate.type"=1,"maxtime"=max(grid),"reverse"=FALSE,"n.iter"=fit$n.iter)
  }
  else{
    R[Rna] <- max(c(L,R)) + 1 ## artificial closure of right censored intervals 
    oL <- order(L)
    L <- L[oL]
    R <- R[oL]
    status <- status[oL]
    peto.intervals  <-  PetoInt(L,R,status)
    indices <- IntIndex(x=peto.intervals,L=L,R=R)
    ## print(unclass(indices))
    Mindex <- indices$Mindex
    Mstrata <- indices$Mstrata
    Iindex <- indices$Iindex
    Istrata <- indices$Istrata
    M <- length(Mstrata)
    N <- length(Istrata)
    Z <- rep(1/M,M)
    fit  <- .C('GMLE',as.integer(c(0,Mstrata)),as.integer(c(0,Istrata)),as.integer(Mindex),as.integer(Iindex),as.integer(N),as.integer(M),Z=as.double(Z),double(length(Z)),as.double(ntol),as.integer(maxiter),steps=integer(1),package="prodlim")
    n.event <- c(0,fit$Z*M)
    surv <- c(1,1-cumsum(fit$Z))
    hazard <- c(0,fit$Z)/surv
    out <- list("call"=NULL,"formula"=NULL,"time"=c(0,peto.intervals[2,]),"peto.intervals"=peto.intervals,"n.risk"=N-n.event,"n.event"=n.event,"n.lost"= rep(NA,M),"hazard"=round(hazard,tol),"surv"=round(surv,tol),"model.response"=NULL,"model"="survival","cens.type"="intervalCensored","covariate.type"=1,"maxtime"=max(c(peto.intervals)),"reverse"=FALSE,"n.iter"=fit$n.iter)
  }
  class(out) <- "prodlim"
  out
}


##     names(L) <- rep("L",N)
##     names(R) <- rep("R",N)
##     ddd <- min(diff(sort(unique(c(L,R)))))*ntol
##     tie <- L==R
##     L <- L+ddd
##     R <- R-ddd
##     L[tie] <- L[tie] - ddd - ddd/2
##     R[tie] <- R[tie] + ddd
##     ##     print(cbind(R,L),digits=10)
##     grid.L <- L[!duplicated(L)]
##     grid.R <- R[!duplicated(R)]
##     grid <- c(grid.L,grid.R)
##     grid <- grid[order(grid)]
##     lg <- length(grid)
##     ng <- names(grid)
##     test <- as.numeric(as.factor(ng))
##     ##     print(rbind(ng,grid,test,c(diff(test),"E")))
##     who.L <- grep("^1$",diff(test))
##     ## print(who.L)
##     mat <- grid[sort(c(who.L,who.L+1))]
##     ## print(mat)
##     ##     maxtime <- max(L,R)
##     petoL <- as.numeric(mat[names(mat)=="L"])
##     petoR <- as.numeric(mat[names(mat)=="R"])
##     grid <- c(t(mat),max(L,R))
##     mindiff <- min(diff(mat))/2
##     indexL <- sindex(petoL,L)
##     indexR <- sindex(petoR,R)
##     ## Ind <- iindex(L,R,c(mat))
##     NS <- length(petoR)
##     ##     print(cbind(L,R,indexL,indexR))
##     ##     stop()
##     fit <- .C("icens_prodlim_ml",
##               as.double(L),
##               as.double(R),
##               as.double(petoL),
##               as.double(petoR),
##               as.integer(indexL),
##               as.integer(indexR),
##               as.integer(status),
##               as.double(N),
##               as.double(NS),
##               nrisk=double(NS),
##               nevent=double(NS),
##               ncens=double(NS),
##               hazard=double(NS),
##               varhazard=double(NS),
##               surv=double(NS),
##               oldsurv=double(NS),
##               as.double(ntol),
##               as.integer(maxiter),
##               maxtime=as.double(max(c(petoR,petoL))),
##               n.iter=integer(1))
