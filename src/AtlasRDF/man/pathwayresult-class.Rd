\name{pathwayresult-class}
\alias{pathwayresult}
\Rdversion{1.1}
\docType{class}

\title{Class \code{"pathwayresult"}}
\description{
Class to store information about a pathway and associated genes.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{pathwayresult()}.
}
\section{Slots}{
  \describe{
    \item{\code{pathwayuri}:}{Object of class \code{"character"} the URI of the pathway from Reactome. }
    \item{\code{label}:}{Object of class \code{"character"} The human readable label of the pathway. }
    \item{\code{numgenes}:}{Object of class \code{"numeric"} A count of the number of genes for this pathway. }
    \item{\code{genes}:}{Object of class \code{"vector"} List of the genes associated with this pathway. }
  }
}
\section{Methods}{
No methods defined with class "pathwayresult" in the signature.
}
\author{
James Malone
}

\examples{
showClass("pathwayresult")
}
\keyword{classes}
\usage{
pathwayresult()
}
