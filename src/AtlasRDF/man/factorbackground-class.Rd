\name{factorbackground-class}
\alias{factorbackground}
\Rdversion{1.1}
\docType{class}

\title{Class \code{"factorbackground"}}
\description{
Used in enrichment to represent an experimental factor and any associated genes for that factor as well as additional ontological information about the factor.
}
\section{Extends}{

All reference classes extend and inherit methods from \code{"\linkS4class{envRefClass}"}.
}

\author{
James Malone, Simon Jupp
}

\examples{
showClass("factorbackground")
}
\keyword{classes}
\section{Fields}{
  \describe{
    \item{\code{uri}:}{Object of class \code{character} URI of the experimental factor. }
    \item{\code{label}:}{Object of class \code{character} name of the experimental factor. }
    \item{\code{species}:}{Object of class \code{character} Species ID of the experimental factor. }
    \item{\code{geneuris}:}{Object of class \code{vector} URIs of the genes for this experimental factor. }
    \item{\code{numgenesexpressed}:}{Object of class \code{integer} Count of the genes expressed for this factor. }
    \item{\code{numgenesnotexpressed}:}{Object of class \code{integer} Count of the genes not expressed for this factor from the total genes for this species.}
    \item{\code{subclasses}:}{Object of class \code{vector} Subclass of this factor from ontology}
    \item{\code{superclasses}:}{Object of class \code{vector} Superclasses of this factor from ontology }
  }
}
\section{Methods}{
  \describe{
    \item{\code{mergegeneuris(value)}:}{ Merge gene URIs into the current vector }
    \item{\code{getsuperclasses()}:}{ Get superclasses for this factor from the ontology }
    \item{\code{setsuperclasses(value)}:}{ Set superclasses for this factor from the ontology }
    \item{\code{getsubclasses()}:}{ Get subclasses for this factor from the ontology }
    \item{\code{setsubclasses(value)}:}{ Set ontology subclasses for this factor }
    \item{\code{getnumgenesnotexpressed()}:}{ Get the number of genes not expressed for this factor as a total from total number of genes for this species }
    \item{\code{setnumgenesnotexpressed(value)}:}{ Set the number of genes not expressed for this factor as a total from total number of genes for this species }
    \item{\code{getnumgenesexpressed()}:}{ Get the number of genes expressed for this factor. }
    \item{\code{setnumgenesexpressed(value)}:}{ Set the number of genes expressed for this factor. }
    \item{\code{getgeneuris()}:}{ Get URIs for genes differentially expressed for this factor. }
    \item{\code{setgeneuris(value)}:}{ Set URIs for genes differentially expressed for this factor. }
    \item{\code{getspecies()}:}{ Get species ID for this factor. }
    \item{\code{setspecies(value)}:}{ Set species ID for this factor. }
    \item{\code{getlabel()}:}{ Get name for this factor. }
    \item{\code{setlabel(value)}:}{ Set name for this factor. }
    \item{\code{geturi()}:}{ Get the URI for this factor. }
    \item{\code{seturi(value)}:}{ Set the URI for this factor. }
  }
}
