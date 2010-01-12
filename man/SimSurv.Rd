\name{SimSurv}
\alias{SimSurv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulating survival data }
\description{
  Censored event times are drawn from user specified conditional distributions
  given simulated covariates.
}
\usage{
SimSurv(N,
	surv,
	cens,
	cova,
	verbose=1,
	...)
}
\arguments{
  \item{N}{Sample size}
  \item{surv}{Dummy argument. The survival distribution is determined by
    arguments of the form \code{surv.arg}. See Details}
  \item{cens}{If FALSE data are left uncensored. Otherwise the censoring
    distribution is specified by using arguments of the form
    \code{cens.arg}. See Details.}
  \item{cova}{either a matrix with \code{N} rows,
    or a named list to generate the covariates. each entry
    is a list where the first element is the names of the
    function used to draw the covariate values and the remaining
    elements are arguments passed to that function. For example
    \code{cova=list(X = list(dist = "rlnorm", meanlog = 2, sdlog = 0.4))}
    generates a single log-normally distributed covariate called
    \code{X}.
  }
  \item{verbose}{Set to \code{FALSE} to shut up in simulations}
  \item{...}{used for convenient argument specification, e.g.
    \code{surv.transform.X2=function(x)x^2}
    will overwrite a corresponding entry
    of \code{surv} for the \code{transform} of the
    covariate \code{X2}
  }
}
\details{
  \code{surv.dist} name of the function to draw the survival distribution
  \code{surv.baseline} the baseline risk
  \code{surv.link} names of the link to the covariates
  \code{surv.coef} numeric vector of regression coefficients. the ith is
    used for the ith entry of \code{cova}
  \code{surv.transform} list of function names one for each covariate. the ith is applied
    to values of the ith covariate before it is linked to generate the
    survival times. 
  \code{cens.dist} name of the function to draw the censoring distribution
  \code{cens.args} list of extra arguments to \code{dist}
  \code{cens.baseline} the baseline risk
  \code{cens.link} names of the link to the covariates
  \code{cens.max} maximal value where all event times are right censored
  \code{cens.type} "right" for right censored data, "interval" for interval
    censored data
  \code{cens.coef} numeric vector of regression coefficients. the ith is
    used for the ith entry of \code{cova}
  \code{cens.transform}list of function names one for each covariate,
    where the ith element is applied
  to values of the ith covariate before it is linked to generate the
  survival times.
  
  Possible distributions to generate covariates:

  X.lognorm = list(dist = "rlnorm", meanlog = 2, sdlog = 0.4)
   
   X.unif = list(dist = "runif", min = 0, max = 10),
   
   X.exp = list(dist = "rexp", rate = 0.4)
   
   X.bernoulli = list(dist = "rbinom", size = 1, prob = 0.3)
   
   X.binom = list(dist = "rbinom", size = 4, prob = 0.8)
   
   X.nbinom = list(dist = "rnbinom", size = 3, mu = 2)
   
   X.poisson = list(dist = "rpois", lambda = 1.3) 
}
\value{
  A data.frame:
  \item{time}{the right censored event times}
  \item{status}{the survival status}
  \item{X}{the values of the covariate \code{X}}
  \item{f.X}{the transformed values of the covariate \code{X}}
  \item{time}{the uncensored event times}

  Attributes of the data.frame
  \code{formula} a formula to evaluate the generated event history object
}
\references{

  Ralf Bender, Thomas Augustin, and Maria Blettner. Generating 
survival times to simulate Cox proportional hazards models by Ralf Bender, 
Thomas Augustin and Maria Blettner, Statistics in Medicine 2005; 
24:1713-1723. Stat Med, 25(11):1978-9, 2006.

}
\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}
% \note{
% }
\seealso{\code{\link{prodlim}}}
\examples{
SimSurv(10)

SurvData=SimSurv(100,cens.baseline=1/10,surv.baseline=2)

Hist(SurvData$time,SurvData$status)

prodlim(Hist(time,status)~1,data=SurvData)

plot(prodlim(Hist(time,status)~1,data=SurvData))

plot(SurvData,atrisk=FALSE,legend=FALSE)

SurvData=SimSurv(100,cens.baseline=1/10,surv.baseline=2,surv.coef=c(-1,-2),
         cova=list( X.exp = list(dist = "rexp", rate = 0.4),
         X.bernoulli = list(dist = "rbinom", size = 1, prob = 0.3)))


SurvData=SimSurv(100,cens.baseline=1/10,surv.baseline=2,surv.coef=c(-1,-2),cens.coef=c(0,1),
          cova=list( X.exp = list(dist = "rexp", rate = 0.4),
          X.bernoulli = list(dist = "rbinom", size = 1, prob = 0.3)))

\dontrun{
coxph(Surv(time,status==0)~X.exp+X.bernoulli,data=SurvData)

coxph(Surv(time,status)~X.exp+X.bernoulli,data=SurvData)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}% at least one, from doc/KEYWORDS
