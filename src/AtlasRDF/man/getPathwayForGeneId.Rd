\name{getPathwayForGeneId}
\alias{getPathwayForGeneId}

\title{
Get pathways associated to a gene.
}
\description{
Get signalling pathways from Reactome that are associated to a gene, via associated protein(s) which are in turn connected to said signalling patwhay(s) in Reactome.
}
\usage{
getPathwayForGeneId(geneid, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{geneid}{
    The ENSEMBL ID for a gene.
    }
    \item{endpoint}{
        SPARQL endpoint to query data from. Default value is the live Atlas SPARQL endpoint but this can be overidden.     
    }
}

\value{
List of pathway URIs.
}

\author{
James Malone, Simon Jupp
}


\examples{
###get pathways for the gene SNHG5
getPathwayForGeneId("http://identifiers.org/ensembl/ENSG00000203875")
}

