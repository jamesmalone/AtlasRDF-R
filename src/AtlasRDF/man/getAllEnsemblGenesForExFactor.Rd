\name{getAllEnsemblGenesForExFactor}
\alias{getAllEnsemblGenesForExFactor}
\title{
Get all ENSEMBL genes reported as differentially expressed for a given experimental factor.
}
\description{
Get all ENSEMBL genes reported as differentially expressed for a given experimental factor. 
    Returns the gene URI, ENSEMBL Id and common gene name.
}
\usage{
getAllEnsemblGenesForExFactor(exfactor, limit, endpoint)
}
\arguments{
  \item{exfactor}{
    The URI of the experimental factor e.g. "efo:EFO_0001420"
}
  \item{limit}{
    Limit the number of genes returned by the query to this, no limit specified or value of 0 returns all.
}
  \item{endpoint}{
    SPARQL endpoint to query data from. Default value is the live Atlas SPARQL endpoint but this can be overidden.     
}

}

\value{
Data frame describing genes for the specified experimental factor with the following slots:
\item{dbXref}{the URI of the gene}
\item{geneName}{the common name of the gene}
\item{ensemblid}{the ENSEMBL id of the gene}
\item{propertyValue}{the property value for this gene corresponding to the experimental factor in the query}
}
\author{
James Malone, Simon Jupp
}

\examples{
###get the genes for the experimental factor 'sepsis' (efo:EFO_0001420)
results <- getAllEnsemblGenesForExFactor("efo:EFO_0001420")
results[1,] #first gene result
}