\name{getExperimentURIsForGeneId}
\alias{getExperimentURIsForGeneId}

\title{
Get experiment URIs which contain a gene specified by ENSEMBL ID.
}
\description{
Get all of the experiment URIs which contain a gene, specified by an ENSEMBL ID, which has been reported as differentially expressed with said experiment. 
}
\usage{
    getExperimentURIsForGeneId(geneid, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{geneid}{
     ENSEMBL gene ID 
}
  \item{endpoint}{
      The location of the SPARQL endpoint - default value is live Atlas endpoint
}

}

\value{
List of experiment URIs for the specified gene. Note, the last part of the URI is the accession number used in Arrayexpress to identify an experiment. 
}

\author{
James Malone, Simon Jupp
}


\examples{
###get the experiment IDs for a gene ENSMUSG00000031314 (ENSEMBL ID)
    getExperimentURIsForGeneId("ENSMUSG00000031314")

}

