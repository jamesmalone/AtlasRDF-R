\name{getGeneUriFromEnsemblId}
\alias{getGeneUriFromEnsemblId}

\title{
Get gene URI from an ENSEMBL ID.
}
\description{
Get the URI for an ENSEMBL gene based on the ENSEMBL ID.
}
\usage{
getGeneUriFromEnsemblId(id, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{id}{
     ENSEMBL gene ID
}
  \item{endpoint}{
      The location of the SPARQL endpoint - default value is live Atlas endpoint
}

}

\value{
List of URIs for a given ENSEMBL ID - this is almost always one but it is feasible there are two so be aware.
}

\author{
James Malone, Simon Jupp
}


\examples{
###get the URI identifier for a gene based on the ENSEMBL ID
#getGeneUriFromEnsemblId("ENSMUSG00000031314")

}

