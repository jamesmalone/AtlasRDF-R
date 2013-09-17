\name{getExFactorURIFromLabel}
\alias{getExFactorURIFromLabel}

\title{
    Get the EFO URI from a given label
}
\description{
    Find the EFO ontology class URI for a given search label. 
}
\usage{
    getExFactorURIFromLabel(label, endpoint="http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
    \item{label}{
        Search string to find the EFO class URI by
    }
    \item{endpoint}{
        SPARQL endpoint to query data from. Default value is the live Atlas SPARQL endpoint but this can be overidden.     
    }
}
\value{
    EFO class URI for the given search term
}

\author{
    James Malone, Simon Jupp
}

\examples{
    ####find the efo class uri for "asthma"
    #getExFactorURIFromLabel("asthma")   #<http://www.ebi.ac.uk/efo/EFO_0000270>
}
