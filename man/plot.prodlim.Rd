\name{plot.prodlim}
\alias{plot.prodlim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting event probabilities over time}
\description{
  Function to plot survival and cumulative incidence curves against
  time. 
}
\usage{
\method{plot}{prodlim}(x, what, cause = 1, newdata, add = FALSE, col, lty, lwd,
             ylim, xlim, xlab = "Time", ylab, legend = TRUE,
             legend.args= NULL, mark.time = FALSE, cex.mark.time = 1.5,
             conf.int = TRUE, conf.int.args = list(col = "gray"),
             atrisk, atrisk.args, timeOrigin,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class `prodlim' as returned by the
          \code{prodlim} function. }
	\item{what}{controls what part of the object is plotted.
	  Defaults to \code{"survival"} for the Kaplan-Meier estimate of
	  the survival function in two state models and to
	  \code{"incidence"} for the Aalen-Johansen estimate of the
	  cumulative incidence function in competing risk models}
  \item{cause}{Determines the cause of the cumulative incidence
    function. When there are strata defined by newdata only one cause is
  possible.}
  \item{newdata}{a data frame containing strata  for
          which plotted curves are desired.}
  \item{add}{if 'TRUE' curves are added to an existing plot. }
  \item{col}{color for curves defaults to 1:number(curves)}
  \item{lty}{line type for curves defaults to 1}
  \item{lwd}{line width for all curves}
  \item{ylim}{limits of the y-axis}
  \item{xlim}{limits of the x-axis}
  \item{ylab}{label for the y-axis}
  \item{xlab}{label for the x-axis}
  \item{legend}{if TRUE a legend is be plotted automatically}
  \item{legend.args}{list of arguments that is  passed to the function \code{legend}}
  \item{mark.time}{if TRUE the curves are tick-marked at right censoring times}
  \item{cex.mark.time}{cex for the tick-marks if mark.time is set TRUE}
  \item{conf.int}{if TRUE pointwise confidence intervals are plotted}
  \item{conf.int.args}{list with extra arguments like col, lty, etc. used for the
    pointwise confidence intervals. Extra argument \code{time}
    determines time points at which confidence intervals are plotted.}
  \item{atrisk}{if TRUE display numbers of subjects at risk. if data are
  cluster correlated also the number of clusters with at least one
  subject at risk is given.}
  \item{atrisk.args}{list of extra arguments for plotting of the numbers
at risk. Recognized extra arguments are: \code{times}: at which time
points numbers are plotted. \code{labels} vector of labels,
e.g. "No. Patients".  }
  \item{timeOrigin}{Start of the time axis}
  \item{\dots}{graphical parameters that are passed to function \code{plot}}
}
\details{
  See examples.
}
\value{
  \code{NULL}
}
\author{Thomas Alexander Gerds <tag@biostat.ku.dk>}
\note{
Similar functionality is provided by the function
\code{\link{plot.survfit}} of the survival library 
}
\seealso{\code{\link{prodlim}},\code{\link{plot.Hist}},\code{\link{summary.prodlim}},
\code{\link{neighborhood}}}

\examples{

## simulate right censored data from a two state model 
dat <- data.frame(time=rexp(100),status=rbinom(100,1,.3),X=rbinom(100,1,.5),Z=rnorm(100,10,3),patnr=sample(1:10,size=100,replace=TRUE))
with(dat,plot(Hist(time,status)))

### marginal Kaplan-Meier estimator
kmfit <- prodlim(Hist(time, status) ~ 1, data = dat)
plot(kmfit)

### Kaplan-Meier in discrete strata
kmfitX <- prodlim(Hist(time, status) ~ X, data = dat)
plot(kmfitX)

### Kaplan-Meier in continuous strata
kmfitZ <- prodlim(Hist(time, status) ~ Z, data = dat)
plot(kmfitZ,newdata=data.frame(Z=c(5,7,12)))

### Cluster-correlated data
kmfitC <- prodlim(Hist(time, status) ~ cluster(patnr), data = dat)
plot(kmfitC,atrisk.args=list(labels=c("Units","Patients")))

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
