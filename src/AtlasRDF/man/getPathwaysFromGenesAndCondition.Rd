\name{getPathwaysFromGenesAndCondition}
\alias{getPathwaysFromGenesAndCondition}

\title{
    Get pathways connected to genes which are differentially expressed for a user specified condition.
}
\description{
    Get list of pathway names that are connected to genes (via proteins which are encoded for by those genes) and which are reported as
    differentially expressed under a specified condition. The condition is specified as an ontology URI. 
}
\usage{
    getPathwaysFromGenesAndCondition(condition, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}
\arguments{
    \item{condition}{
        URI of a condition based on the EFO ontology class. 
    }
    \item{endpoint}{
        The location of the SPARQL endpoint - default value is live Atlas endpoint
    }
    
}
\value{
    Returns a data frame with the following:
    \item{pathwayname}{The name of the pathway}
    \item{pathway }{The URI identifier of the pathway}
    \item{expressionvalue}{The gene and expression for the pathway}
    \item{pvalue}{p-value of the differential expression of the gene for this pathway in this condition}

}

\author{
    James Malone, Simon Jupp
}

\examples{
    ###get pathways connected to genes which are differentially expressed for HeLa (efo:EFO_00011185)
    #getPathwaysFromGenesAndCondition("efo:EFO_0001185")
    
}

