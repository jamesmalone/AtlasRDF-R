\name{getTaxonURI}
\alias{getTaxonURI}

\title{
Get URI of a taxon class based on the common or taxonomic name.
}
\description{
Convenience method to the ontology URI, used in other queries, based on an input common or taxonomic name such as "human" or "Homo sapiens"
}
\usage{
getTaxonURI(taxonName)
}

\arguments{
  \item{taxonName}{
Common or taxonomic name
}
}

\value{
URI of the taxon
}

\author{
James Malone, Simon Jupp
}

\examples{
###get URI for human
getTaxonURI("human")
    
###or for arabidopsis
getTaxonURI("arabidopsis")

}

