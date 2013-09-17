\name{vizPvalues}
\alias{vizPvalues}

\title{
Visualize the results of an enrichment filtering for a specificed p-value.
}
\description{
Produce a plot of the most enriched experimental factors following an Atlas gene set enrichment. The cutoff allows a filter to be placed by p-value.
}
\usage{
vizPvalues(resultset, cutoff = "0.05")
}

\arguments{
  \item{resultset}{
        A set of results following an AtlasRDF gene set enrichment. 
}
  \item{cutoff}{ 
        The p-value threshold at which to include a gene for visualization
}
}

\value{
    Plot of results
}
\author{
James Malone
}


\examples{
###visualize the results from a gene set enrichment performed using AtlasRDF package
#vizPValues(enrichmentresults, cutoff = "0.02")

}

