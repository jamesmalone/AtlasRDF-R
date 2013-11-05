\name{searchForEFOTerms}
\alias{searchForEFOTerms}

\title{
    Search for EFO URIs from a given label text
}
\description{
    Finds any EFO ontology class whose label matches the given search label or part of it. 
}
\usage{
    searchForEFOTerms(label, endpoint="http://www.ebi.ac.uk/rdf/services/atlas/sparql")
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
    Data frame of EFO class URIs and their corresponding labels
}

\author{
    James Malone, Simon Jupp
}

\examples{
    ###find the efo classes which contain the word 'cancer' in their label
    searchForEFOTerms("cancer")   
}
