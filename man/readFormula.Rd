\name{readFormula}
\alias{readFormula}
\title{
Read a special formula 
}
\description{
Distangle formula with special terms. 
}
\usage{
readFormula(formula, specials, specialArgumentNames, alias, unspecified = "unSpec")
}
\arguments{
  \item{formula}{
    A formula of the form \code{Hist(time, event) ~ g(a) + f(b) +
      g(c)} where \code{f} and \code{g} are examples of special
    functions.
  }
  \item{specials}{
    Names of the special functions. 
  }
  \item{specialArgumentNames}{
    Further arguments of special functions.
  }
  \item{alias}{
 List of alias names for special functions 
  }
  \item{unspecified}{
    Special function for variables that are not specialized.
  }
}
\details{
  This function restricts the formula syntax to some extend. For example,
  it is not possible to construct variables inside the formula:
  
  This is an error:
  
  \code{readFormula(a~b+X(I(c>4)),specials=c("X"))}

  Instead, one need to construct the variable \code{data$c4  <-
    data$c>4} and then use

  \code{readFormula(a~b+X(c),specials=c("X"))}
}
\value{
A list 
  \item{Response}{A lhs-formula, i.e. a character-class formula where the right hand set to \code{~1}}
  \item{Intercept}{The Intercept of the formula}
  \item{special}{For each special name and the unspecificied name, a
    list with a character-class formula and specialArguments.}
 \item{allVars}{A character-class rhs-formula containing all variables
   of argument formula}
}
\author{Thomas Gerds <tag@biostat.ku.dk>}
\examples{

readFormula(y~x+f(z),specials=c("f"),unspecified="g")

\dontrun{
## suitable construction to serve cmprsk
qFun <- function(x){x^2}
## multiply X1 by time^2 and X2 by time:
readFormula(Hist(time,event)~cov2(X1,tf=qFun)+cov2(X2),
	    specials=c("cov2","cov1"),
	    specialArgumentNames=list("cov2"="tf"),
	    unspecified="cov1")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
