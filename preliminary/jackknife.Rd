\name{jackknife}
\Rdversion{1.1}
\alias{jackknife}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Jackknife pseudo-values based on the product-limit method
}
\description{
  Jackknife pseudo-values based on the Kaplan-Meier and Aalen Johannsen estimates.
}
\usage{
jackknife(object, times, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    %%     ~~Describe \code{object} here~~
}
  \item{times}{
%%     ~~Describe \code{times} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

## Survival
dd <- prodlim:::SimSurv(100)
ff <- prodlim(Hist(time,status)~1,data=dd)
jj <- jackknife.survival(ff,times=c(10,30,100))

\dontrun{
dd <- prodlim:::SimSurv(1000)
library(pseudo)
system.time(jack <- with(dd,pseudosurv(time,status,c(10,30,100))))
system.time({ff <- prodlim(Hist(time,status)~1,data=dd);jack2 <- jackknife(ff,times=c(10,30,100))})
all(round(jack2[,1],9)==round(jack[[3]],9))
}

## CompRisk
ddd <- data.frame(time=1:100,event=rbinom(100,2,.5),X=rbinom(100,1,.5))
fff <- prodlim(Hist(time,event)~1,data=ddd)
fk <- leaveOneOut.competing.risks(fff,times=c(10,30,100,500))
jk <- jackknife.competing.risks(fff,times=c(10,30,100))

\dontrun{

library(prodlim)
library(pseudo)
ddd <- data(sampleCompriskFrame)
fff <- prodlim(Hist(time,event)~1,data=ddd)
fff7 <- prodlim(Hist(time,event)~1,data=ddd[-7,])
fff1 <- prodlim(Hist(time,event)~1,data=ddd[-1,])
ttt <- ddd$time
p1 <- predict(fff1,times=ttt)
p7 <- predict(fff7,times=ttt)
summary(fff1,times=1:10,cause=1)
sf <- leaveOneOut.competing.risks(fff,times=ttt)
cbind(loo=sf[7,],true=p7)
cbind(loo=sf[1,],true=p1)
jack <- with(ddd,pseudoci(time,event,ttt))[,-c(1:2,seq(4,22,2))]
jack2 <- jackknife.competing.risks(fff,times=ttt)
all(sapply(1:10,function(x){
a <- round(jack[,x],4)
b <- round(jack2[,x],4)
# all(a[!is.na(a)]==b[!is.na(b)])
all(a[!is.na(a)]==b[!is.na(a)])
}))


N <- 2000
ddd <- data.frame(time=1:N,event=rbinom(N,2,.5),X=rbinom(N,1,.5))
ttt <- c(3,5,10)
# ttt <- ddd$time
fff <- prodlim(Hist(time,event)~1,data=ddd)
system.time(jack <- with(ddd,pseudoci(time,event,ttt))[,-c(1:2,seq(4,N+2,2))])
system.time({jack2 <- jackknife.competing.risks(fff,times=ttt)})
all(round(jack2[,1],9)==round(jack[,1],9))

all(sapply(1:N,function(x){
a <- round(jack[x,],4)
b <- round(jack2[x,],4)
# all(a[!is.na(a)]==b[!is.na(b)])
all(a[!is.na(a)]==b[!is.na(a)])

}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
