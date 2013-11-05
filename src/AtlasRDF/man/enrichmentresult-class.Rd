\name{enrichmentresult-class}
\alias{enrichmentresult}

\Rdversion{1.1}
\docType{class}


\title{Class \code{"enrichmentresult"}}
\description{
Class to store the enrichment for a single experimental factor for a given gene set.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{enrichmentresult()}.
}
\section{Slots}{
  \describe{
    \item{\code{factoruri}:}{Object of class \code{"character"} the URI of the experimental factor}
    \item{\code{label}:}{Object of class \code{"character"} human readable label of the experimental factor }
    \item{\code{p.value}:}{Object of class \code{"numeric"} p.value for enrichment for this factor }
    \item{\code{estimate}:}{Object of class \code{"numeric"} see Fishers exact test }
    \item{\code{alternative}:}{Object of class \code{"character"} see Fishers exact test }
    \item{\code{null.value}:}{Object of class \code{"numeric"} see Fishers exact test }
    \item{\code{method}:}{Object of class \code{"character"} Fishers exact test }
    \item{\code{enrichedgenes}:}{Object of class \code{"vector"} Vector of the genes enriched for this factor from the gene list }
  }
}
\section{Methods}{
No methods defined with class "enrichmentresult" in the signature.
}
\author{
James Malone
}

\examples{
showClass("enrichmentresult")
}

\keyword{classes}

\usage{
enrichmentresult()
}

