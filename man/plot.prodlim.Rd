\name{plot.prodlim}
\alias{plot.prodlim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting event probabilities over time}
\description{
  Function to plot survival and cumulative incidence curves against
  time. 
}
\usage{
\method{plot}{prodlim}(x, type, cause = 1, newdata, add = FALSE, col, lty, lwd,
             ylim, xlim, xlab = "Time", ylab, legend = TRUE,
             marktime = FALSE, confint = TRUE, automar, atrisk=ifelse(add,FALSE,TRUE), timeOrigin,
axes=TRUE, percent=FALSE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class `prodlim' as returned by the
    \code{prodlim} function. }
  \item{type}{controls what part of the object is plotted.
    Defaults to \code{"survival"} for the Kaplan-Meier estimate of
    the survival function in two state models and to
    \code{"incidence"} for the Aalen-Johansen estimate of the
    cumulative incidence function in competing risk models}
  \item{cause}{determines the cause of the cumulative incidence
    function. Currently one cause is allowed at a time, but you may call
    the function again with add=TRUE to add the lines of the
    other causes.}
  \item{newdata}{a data frame containing strata for
    which plotted curves are desired.}
  \item{add}{if 'TRUE' curves are added to an existing plot. }
  \item{col}{color for curves defaults to 1:number(curves)}
  \item{lty}{line type for curves defaults to 1}
  \item{lwd}{line width for all curves}
  \item{ylim}{limits of the y-axis}
  \item{xlim}{limits of the x-axis}
  \item{ylab}{label for the y-axis}
  \item{xlab}{label for the x-axis}
  \item{legend}{if TRUE a legend is plotted by calling the function legend. 
    Optional arguments of the function \code{legend} can be given in the
    form \code{legend.x=val} where x is the name of the argument and
    val the desired value. See also Details.}
  \item{marktime}{if TRUE the curves are tick-marked at right censoring
    times by invoking the function \code{MarkTime}. Optional arguments of the
    function \code{MarkTime} can be given in the form
    \code{confint.x=val} as with legend. See also Details.}
  \item{confint}{if TRUE pointwise confidence intervals are plotted by
    invoking the function \code{ConfInt}. Optional arguments of the
    function \code{ConfInt} can be given in the form
    \code{confint.x=val} as with legend. See also Details.}
  \item{automar}{If TRUE the function trys to get good values for
    figure margins around the main plotting region.}
  \item{atrisk}{if TRUE display numbers of subjects at risk by invoking the function \code{AtRisk}. Optional arguments of the
    function \code{AtRisk} can be given in the form
    \code{atrisk.x=val} as with legend. See also Details.}
  \item{timeOrigin}{Start of the time axis}
  \item{axes}{If true axes are drawn.}
  \item{percent}{If true the y-axis is labeled in percent.}
  \item{\dots}{graphical parameters that are passed to function
    \code{plot}. }
}
\details{
  From version 1.1.3 on the arguments legend.args, atrisk.args,
  confint.args are obsolete and only available for backward
  compatibility. Instead arguments for the invoked functions
  \code{AtRisk}, \code{legend},  \code{ConfInt}, \code{MarkTime},
  \code{axis}
  are simply specified as \code{atrisk.cex=2}. The specification is not
  case sensitive, thus
  \code{AtRisk.cex=2} or \code{atRISK.cex=2} will have the same effect.
  The function \code{axis} is called twice, and arguments of the form
  \code{axis1.labels}, \code{axis1.at} are used for the time axis
  whereas \code{axis2.pos}, \code{axis1.labels}, etc. are used for the
  y-axis.

  These arguments are processed via  \code{\dots} of \code{plot.prodlim}
  and inside by using the function \code{resolveSmartArgs}.
  Documentation of these arguments can be found 
  in the help pages of the corresponding
  functions. 
}
\value{
  The (invisible) object.
}
\author{Thomas Alexander Gerds <tag@biostat.ku.dk>}
\note{
Similar functionality is provided by the function
\code{\link{plot.survfit}} of the survival library 
}
\seealso{\code{\link{prodlim}},\code{\link{plot.Hist}},\code{\link{summary.prodlim}},
\code{\link{neighborhood}}, \code{\link{AtRisk}}, \code{\link{ConfInt}},
\code{\link{MarkTime}}}

\examples{

## simulate right censored data from a two state model 
dat <- data.frame(time=rexp(100),status=rbinom(100,1,.3),X=rbinom(100,1,.5),Z=rnorm(100,10,3),patnr=sample(1:10,size=100,replace=TRUE))
with(dat,plot(Hist(time,status)))

### marginal Kaplan-Meier estimator
kmfit <- prodlim(Hist(time, status) ~ 1, data = dat)
plot(kmfit)
plot(kmfit,percent=TRUE)
plot(kmfit,percent=TRUE,axis1.at=c(0,2,5),axis1.pos=0,axis2.pos=0)

### Kaplan-Meier in discrete strata
kmfitX <- prodlim(Hist(time, status) ~ X, data = dat)
plot(kmfitX)
plot(kmfitX,legend.x="bottomleft",AtRisk.cex=1.3)

### Kaplan-Meier in continuous strata
kmfitZ <- prodlim(Hist(time, status) ~ Z, data = dat)
plot(kmfitZ,newdata=data.frame(Z=c(5,7,12)))

### Cluster-correlated data
kmfitC <- prodlim(Hist(time, status) ~ cluster(patnr), data = dat)
plot(kmfitC,atrisk.labels=c("Units","Patients"))

## simulate right censored data from a competing risk model 
datCR <- data.frame(time=rexp(100),status=rbinom(100,2,.3),X=rbinom(100,1,.5),Z=rnorm(100,10,3))
with(datCR,plot(Hist(time,status)))

### marginal Aalen-Johansen estimator
ajfit <- prodlim(Hist(time, status) ~ 1, data = datCR)
plot(ajfit)

### conditional Aalen-Johansen estimator
ajfitXZ <- prodlim(Hist(time, status) ~ X+Z, data = datCR)
plot(ajfitXZ,newdata=data.frame(X=c(1,1,0),Z=c(4,10,10)))
plot(ajfitXZ,newdata=data.frame(X=c(1,1,0),Z=c(4,10,10)),cause=2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
