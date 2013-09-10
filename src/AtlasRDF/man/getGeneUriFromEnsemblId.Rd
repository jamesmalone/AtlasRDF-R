\name{getGeneUriFromEnsemblId}
\alias{getGeneUriFromEnsemblId}

\title{
Get gene URI from an ENSEMBL ID.
}
\description{
Get the URI for an ENSEMBL gene based on the ENSEMBL ID.
}
\usage{
getGeneUriFromEnsemblId(id, taxon, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{id}{
     ENSEMBL gene ID
}
  \item{taxon}{
        The URI of the species of interest
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.


}

