\name{getLabel}
\alias{getLabel}

\title{
Get label of an entity based on the URI.
}
\description{
Get label of an entity base on the URI
}
\usage{
getLabel(uri, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{uri}{
URI of the class for which the label is required.
}
  \item{endpoint}{
        The location of the SPARQL endpoint - default value is live Atlas endpoint
  }

}

\value{
List of labels for the class URI
}

\author{
James Malone, Simon Jupp
}

\examples{
###get the label for the class efo:EFO_0001185 (HeLa)
#getLabel("efo:EFO_0001185")

}

