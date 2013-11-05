\name{orderEnrichmentResults}
\alias{orderEnrichmentResults}

\title{
Order the results of a gene set enrichment by p-value.
}
\description{
Following a gene set enrichment, order the results so that the most significantly enriched factors are first, in decreasing order of significane (i.e. where lowest p-value is most signfiicant)
}

\usage{
orderEnrichmentResults(resultset)
}

\arguments{
  \item{resultset}{
        A set of results following an AtlasRDF gene set enrichment. 
}
}

\value{
    Ordered set of results with most signficant factor first, least signficant last.
}
\author{
James Malone
}


\examples{
###order the results of an enrichment
#orderedresults <- orderEnrichmentResults(enrichmentresults)
###view the top 10 most significant
#vizPValues(orderedresults[1:20], cutoff = "1")

}

