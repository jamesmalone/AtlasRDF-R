\name{generef-class}
\alias{generef}
\Rdversion{1.1}
\docType{class}

\title{Class \code{"generef"}}
\description{
Class to store gene references, including the URI, name and experimental factors associated with the gene.
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.

}
\author{
James Malone, Simon Jupp
}

\examples{
showClass("generef")
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{geneuri}:}{Object of class \code{character} the URI of the gene. }
    \item{\code{genelabel}:}{Object of class \code{character}Gene name. }
    \item{\code{geneensemblid}:}{Object of class \code{character} ENSEMBL ID of the gene, if available. }
    \item{\code{species}:}{Object of class \code{character} Species ID of the gene. }
    \item{\code{exfactoruris}:}{Object of class \code{vector} Exerimental factors associated to this gene. }
  }
}
\section{Methods}{
  \describe{
    \item{\code{mergeexfactoruris(value)}:}{ merge an experimental factor(s) into the list }
    \item{\code{getexfactoruris()}:}{ get the list of experimental factors for this gene }
    \item{\code{getspecies()}:}{ Get the species ID for this gene }
    \item{\code{setspecies(value)}:}{ Set the species ID for this gene }
    \item{\code{getensemblid()}:}{ Get the ENSEMBL ID for this gene. }
    \item{\code{setensemblid(value)}:}{ Set the ENSEMBL ID for this gene. }
    \item{\code{getgenelabel()}:}{ Get gene name. }
    \item{\code{setgenelabel(value)}:}{Set gene name. }
    \item{\code{getgeneuri()}:}{ Get the URI for the gene. }
    \item{\code{setgeneuri(value)}:}{ Set the URI for the gene. }
  }
}
