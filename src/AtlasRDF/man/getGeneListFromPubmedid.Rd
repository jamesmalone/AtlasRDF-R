\name{getGeneListFromPubmedid}
\alias{getGeneListFromPubmedid}

\title{
Get gene lists in Atlas for a given pubmed id (assuming the data is in Atlas)
}
\description{
Get differentially expressed genes reported in an experiment in Atlas for which a given pubmed id is linked (assuming the data is in Atlas)
}
\usage{
getGeneListFromPubmedid(searchid, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{searchid}{
   PubMed ID to find data for 
    }
    \item{endpoint}{
        SPARQL endpoint to query data from. Default value is the live Atlas SPARQL endpoint but this can be overidden.     
    }
}

\value{
List of genes for experiment (via PubMed ID)
}

\author{
James Malone, Simon Jupp
}


\examples{
###get gene list for the publication with pubmed id 19668377
getGeneListFromPubmedid("19668377")
}

