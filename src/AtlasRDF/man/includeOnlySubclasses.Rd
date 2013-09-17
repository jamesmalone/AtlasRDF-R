\name{includeOnlySubclasses}
\alias{includeOnlySubclasses}
Filter enrichment results to only include experimental factors of a given parent class (e.g. all cancers)
\title{
Function which allows filtering of enrichment results to remove experimental factors apart from those factors which are a subclasses of a given parent class.
For instance, include only disease factors or only cell types.
}
\usage{
includeOnlySubclasses(includeparentclass, resultset, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{includeparentclass}{
The URI of the parent class for which all experimental factor subclasses will be included (removing all others).
}
  \item{resultset}{
The full enrichment result set
}
    \item{endpoint}{
        The location of the SPARQL endpoint - default value is live Atlas endpoint
    }
}

\value{
New Atlas enrichment result set, filtered to include only experimental factors specified in the function call. 
}

\author{
James Malone, Simon Jupp
}

\examples{
###include only experimental factors which are diseases in the enrichment results
#filteredgenes <- includeOnlySubclasses("efo:EFO_0000408", mouseresults, endpoint="http://wwwdev.ebi.ac.uk/rdf/services/atlas/sparql")
###view results
#vizPvalues(filteredgenes)      
}
