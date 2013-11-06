\name{getExperimentIdsForGeneURI}
\alias{getExperimentIdsForGeneURI}

\title{
Get ArrayExpress experiment ID which contain a gene specified by URI.
}
\description{
Get all of the ArrayExpress experiment IDs which contain a gene, specified by a URI.
}
\usage{
    getExperimentIdsForGeneURI(geneuri, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{geneuri}{
     Gene URI
}
  \item{endpoint}{
      The location of the SPARQL endpoint - default value is live Atlas endpoint
}

}

\value{
List of ArrayExpress experiment IDs for the specified gene.
}

\author{
James Malone, Simon Jupp
}


\examples{
###get the experiment IDs for a gene ENSG00000142556 (ENSEMBL ID)
    getExperimentIdsForGeneURI("http://identifiers.org/ensembl/ENSG00000142556")

}

