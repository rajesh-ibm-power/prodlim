\name{plot.Hist}
\alias{plot.Hist}
\title{
Box-arrow diagrams for multi-state models
}
\description{
  Automated plotting of the states and transitions that characterize a
  multi states model.
}

\usage{

\method{plot}{Hist}(x,
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
	  ...)
}

\arguments{
  \item{x}{
    An object of class \code{Hist}.
  }
  
  \item{layout}{
   A list saying how to arrange the boxes.
   The list should contain
   the number of graphic columns (\code{ncol})
   the number of graphic rows (\code{nrow})
   the box positions (\code{box.pos}).
   See the examples below.
 }
 
 \item{xbox.rule}{
   Numeric value specifying the distance between the text and the box
   in x-direction.
 }
 
 \item{ybox.rule}{
   Numeric value specifying the distance between the text and the box
   in y-direction.  }
 
 \item{state.lab}{
   Vector of text that appears in the boxes.
   Defaults to attr(x,"state.names").
 }

 \item{state.cex}{
   Cex for the \code{state.lab}
 }

 \item{rect.args}{
   List of arguments that are
   passed to the call of the function \code{\link{rect}}.
 } 
 
 \item{args.state.lab}{
   List of arguments that are
   passed to the call of the function \code{\link{text}}
   for labeling the inside of the boxes.
 }

 \item{arrow.lab}{
   Vector of labels for the arrows.
 }
 
 \item{arrow.lab.offset}{
   List of vectors of length two: the horizontal and the vertical offset
   of the arrow label from the arrow mid. Length should be equal to the number of arrows.
   Can be a single numeric value which is then used repeatedly.
   If missing the strwidth and strheight of the current label is used.
 } 

 \item{arrow.lab.side}{
   Vector of integers one for each arrow label. Can 
   either be one or 
   minus one. Minus one indicates that the label
   should appear on the other side of the arrow. 
}


\item{arrow.lab.cex}{
  Cex for arrow labels.
}

\item{arrow.lab.style}{
  Choice between "count", "symbolic", "character","no".
  
  "count": the arrows are labeled using the
  number of transitions

  "symbolic":
  the expression \code{alpha[hj](t)} is used to label the
  arrow between state h and state j
  
  "character": expects labels from argument \code{arrow.lab}
  
  "no": omit arrow labels
}

\item{arrows.args}{
  A list with further arguments that are
  passed to the function \code{arrows}.
}

\item{arrow.lab.args}{
  A list with further arguments that are
  passed to the function \code{\link{text}} for labeling
  the arrows.
}

\item{arrow.head.offset}{
  Distance from the arrow heads to the boxes.
}


\item{arrow.double.dist}{
  Distance between arrows that go in both directions.
}

\item{arrow.fix.code}{
  Vector of integers. 
  If for some reason the direction of one or several arrows
  is wrong the nth value of \code{arrow.fix.code} will be
  used for the nth arrow. See argument \code{code} of
  \code{\link{arrows}}.
}

\item{enumerate.boxes}{
  Logical. If TRUE the boxes are numbered in the upper
  left corner.
}


\item{box.numbers}{
  Vector of strings that appear in the
  upper left corner of the boxes,
  if \code{enumerate.boxes==TRUE}. Defaults to
  \code{0,...,number.of.boxes}.
}

\item{cex.boxlabs}{
  Cex for the numbers in the upper left corners.
}


\item{margin}{
  Set the figure margin via \code{par(mar=margin)}
}

\item{verbose}{
  If TRUE echo various things. 
}

\item{reverse}{
  If true the order of the boxes is reversed.
  For competing risk models only the order of the
  absorbing states.
}

\item{\dots}{
  for compatibility with other plot functions
}
}
\details{
  The default layout works only for two-state, competing
  risk and illness-death models.
  All the cex parameters and xy.rules depend on the device. 
}
\value{
  The first argument but in the invisible cloak.
}

\author{
  Thomas Gerds \email{tag@biostat.ku.dk}
}

\note{
  A cool thing would be to use the functionality of
  the unix program `dot' http://www.graphviz.org/About.php to
  obtain more complex graphs.
}
\seealso{
\code{\link{Hist}}
}
\examples{

## A simple survival model
\dontrun{
library(survival)
data(pbc)
SurvFrame <- pbc
## or simply create some data
SurvFrame <- data.frame(time=1:10,status=sample(0:1,10,TRUE))
SurvHist <- with(SurvFrame,Hist(time,status))
plot(SurvHist,
     state.lab=c("Randomization","Death"),
     arrow.lab=expression(lambda(t)),
     layout=list(nrow=2,ncol=2,box.pos=list(c(1,1),c(2,2))))

## two competing risks
comprisk.model <- data.frame(time=1:3,status=1:3)
CRHist <- with(comprisk.model,Hist(time,status,cens.code=2))
plot(CRHist,
     cex=2,
     state.lab=c("Alive","Dead\n cause 1","Dead\n cause 2"),
     style="character",ybox.rule=1,
     arrow.lab=c(expression(gamma[1](t),gamma[2](t))),
     enumerate.boxes=TRUE,
     cex.boxlabs = 1.28)

## change the layout
plot(CRHist,
     cex=2,
     state.lab=c("Alive","Dead\n cause 1","Dead\n cause 2"),
     style="character",ybox.rule=1,
     arrow.lab=c(expression(gamma[1](t),gamma[2](t))),
     arrow.lab.side=c(1,-1),
     enumerate.boxes=TRUE,
     cex.boxlabs = 1.28,
     layout=list(nrow=2,ncol=3,box.pos=list(c(1,2),c(2,1),c(2,3))))

## more competing risks
comprisk.model2 <- data.frame(time=1:4,status=1:4)
CRHist2 <- with(comprisk.model2,Hist(time,status,cens.code=2))
plot(CRHist2,
     cex=2,
     state.lab=c("Alive","Dead\n cause 1","Dead\n cause 2","Dead\n cause 3"),
     arrow.lab=c(expression(alpha[1](t),alpha[2](t),alpha[3](t))),
     enumerate.boxes=TRUE,
     ybox.rule=.8,
     verbose=FALSE,
     cex.boxlabs = 1.28)

## change the layout
plot(CRHist2,
     cex=2,
     state.lab=c("Alive","Dead\n cause 1","Dead\n cause 2","Dead\n cause 3"),
     arrow.lab=c(expression(alpha[1](t),alpha[2](t),alpha[3](t))),
     enumerate.boxes=TRUE,
     ybox.rule=.8,
     verbose=FALSE,
     cex.boxlabs = 1.28,
     layout=list(nrow=3,ncol=3,box.pos=list(c(2,1),c(1,2),c(3,2),c(2,3))),
     reverse=FALSE,
     arrow.lab.side=c(-1,1,-1),
     arrow.lab.offset=c(7,9,8))

## illness-death models
illness.death.frame <-
         data.frame(time=1:4,
		    from=c("Disease-free","Disease-free",
		      "Diseased","Disease-free"),
		    to=c("0","Diseased","Dead","Dead"))
IDHist <- with(illness.death.frame,
	       Hist(time,event=list(from,to)))
plot(IDHist,
     ybox.rule=4,
     xbox.rule=.3,
     state.cex=1.3,
     enum=TRUE,
     arrow.lab.side=c(-1,-1,1))

## illness-death with recovery
illness.death.frame2 <- data.frame(time=1:5,
				   from=c("Disease\nfree","Disease\nfree",
                                     "Diseased","Diseased","Disease\nfree"),
				   to=c("0","Diseased","Disease\nfree",
                                                           "Dead","Dead"))
IDHist2 <- with(illness.death.frame2,Hist(time,event=list(from,to)))
plot(IDHist2,
     ybox.rule=1.3,
     xbox.rule=.3,
     state.cex=2,
     arrow.lab.offset=c(13,13,8,10),
     enum=TRUE,
     verbose=FALSE)
## change the layout

plot(IDHist2,
     ybox.rule=1.3,
     xbox.rule=.3,
     state.cex=2,
     enum=TRUE,
     verbose=FALSE,
     layout=list(ncol=3,nrow=2,box.pos=list(c(1,1),c(2,2),c(1,3))),
     arrow.lab.side=c(-1,1,-1,1),
     arrow.lab.offset=c(15,15,10,10))
 }
}

\keyword{survival}
