\name{getGeneUriFromName}
\alias{getGeneUriFromName}

\title{
    Get URI of a gene based on the common gene name.
}
\description{
   Get the URI of an ENSEMBL gene based on a search by the common gene name.
}
\usage{
    getGeneUriFromName(genename, taxon, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
    \item{genename}{
        The common gene name by which to search
    }
    \item{taxon}{
        The URI of the taxon to search genes for, use function getTaxonURI to get the URI
    }
    \item{endpoint}{
        The location of the SPARQL endpoint - default value is live Atlas endpoint
    }
}

\value{
list of gene URIs matching the query: note this can be more than one
}

\author{
   James Malone, Simon Jupp
}

\examples{
###get the gene URI identifier based on the common gene name for "SRGN" in human
#getGeneUriFromName("SRGN", "obo:NCBITaxon_9606")
}

