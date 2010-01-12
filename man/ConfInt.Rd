\name{ConfInt}
\Rdversion{1.1}
\alias{ConfInt}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Add point-wise confidence limits to the graphs of Kaplan-Meier and
  Aalen-Johannsen estimates of survival and cumulative incidence.
}
\description{
  This function is invoked and controlled by \code{plot.prodlim}. 
}
\usage{
ConfInt(x, times, newdata, type, cause, col, lty, lwd, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{an object of class `prodlim' as returned by the
    \code{prodlim} function. }
  \item{times}{
    where to compute point-wise confidence limits
  }
  \item{newdata}{
    see \code{plot.prodim}
  }
  \item{type}{
    If \code{"bars"} then confidence limits are drawn as vertical bars
    at specified \code{times}. Else (the default) pseudo confidence
    bands are drawn by connecting all the upper confidence
 limits and all the lower confidence
 limits.
  }
  \item{cause}{
    see \code{plot.prodim}
  }
  \item{col}{
    the colour of the lines.
  }
  \item{lty}{
    the line type of the lines.
  }
  \item{lwd}{
    the line thickness of the lines.
  }
  \item{\dots}{
    Further arguments that are passed to the function \code{segments} if
    \code{type=="bars"} and to \code{lines} else.
  }
}
\details{
  This function should not be called directly. The arguments can be
  specified as \code{Confint.arg} in the call to \code{plot.prodim}.
}
\value{
  Nil
}
\author{
Thomas Alexander Gerds <tag@biostat.ku.dk>
}
\seealso{
  \code{\link{plot.prodlim}}, \code{\link{AtRisk}}, \code{\link{MarkTime}}
}
\keyword{survival}
