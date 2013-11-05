\name{getGenesForPathwayURI}
\alias{getGenesForPathwayURI}

\title{
Get genes associated with a signalling pathway from Reactome.
}
\description{
Get genes that are associated with proteins which are in turn associated to a specific signalling pathway in Reactome.
}
\usage{
getGenesForPathwayURI(pathwayuri, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{pathwayuri}{
	The URI of a reactome pathway.
    }
    \item{endpoint}{
        SPARQL endpoint to query data from. Default value is the live Atlas SPARQL endpoint but this can be overidden.     
    }
}

\value{
List of gene URIs for the pathway
}

\author{
James Malone, Simon Jupp
}


\examples{
###get gene list for the pathway AKT Signaling in Cancer
getGenesForPathwayURI("<http://identifiers.org/reactome/REACT_147727.2>")
}

