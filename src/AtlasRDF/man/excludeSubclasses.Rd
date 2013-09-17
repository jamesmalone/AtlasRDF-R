\name{excludeSubclasses}
\alias{excludeSubclasses}
\title{
Exclude subclasses of factors from enrichment results
}
\description{
Removes all factors which are subclasses of a specified type from Atlas gene set enrichment results. 
For instance, if all chemical compounds were to be excluded form the enrichment results or all subclasses of cancer use this function.
}
\usage{
excludeSubclasses(filterparentclass, resultset, endpoint="http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{filterparentclass}{
        The class from which all subclasses will be excluded - including this class.
    }
  \item{resultset}{
        Result set from an Atlas gene set enrichment.
}
    \item{endpoint}{
        SPARQL endpoint to query data from. Default value is the live Atlas SPARQL endpoint but this can be overidden.     
    }

}

\value{
new enrichment result set with the specified class and subclasses excluded. 
}

\author{
James Malone, Simon Jupp
}


\examples{
####filter out all subclasses of chemical compounds from the result set
#####enrichment results has been produced using the enrichment functions in AtlasRDF package
#filteredresults <- excludeSubclasses("obo:CHEBI_37577", enrichmentresults)
#vizPvalues(filteredresults, 0.01)
}

