ConfInt <- function(x,times,newdata,type,cause,col,lty,lwd,...){
  sumx <- summary(x,times=times,newdata=newdata,cause=cause,verbose=FALSE,surv=ifelse(type=="cuminc",FALSE,TRUE))
  if (x$model=="survival" && x$covariate.type==1) sumx <- list(sumx)
  if (x$model=="competing.risks" && x$covariate.type>1) sumx <- sumx[[1]]
  nlines <- length(sumx)
  ci <- lapply(sumx,function(u){
    uu <- data.frame(u[,c("time","lower","upper")])
    uu[(uu$upper-uu$lower)<1,]
  })
  nix <- lapply(1:nlines,function(i){
    if (type=="bars"){
      segments(x0=ci[[i]]$time,x1=ci[[i]]$time,y0=ci[[i]]$lower,y1=ci[[i]]$upper,lwd=lwd[i],col=col[i],lty=lty[i],...)
    }
    else{
      lines(x=ci[[i]]$time,ci[[i]]$lower,type="s",lwd=lwd[i],col=col[i],lty=lty[i],...)
      lines(x=ci[[i]]$time,ci[[i]]$upper,type="s",lwd=lwd[i],col=col[i],lty=lty[i],...)
    }
  })
}
