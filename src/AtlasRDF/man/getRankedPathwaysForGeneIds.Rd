\name{getRankedPathwaysForGeneIds}
\alias{getRankedPathwaysForGeneIds}

\title{
    Get pathway names and URIs for given gene list, sorted into list with most common pathways first.
}
\description{
    Given a list of ENSEMBL gene IDs, the function will find pathways which these genes are involved in. 
    These pathways will be grouped and ranked according to those with highest number of genes in common.
}
\usage{
getRankedPathwaysForGeneIds(genelist, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql") 
}
\arguments{
    \item{genelist}{
       The list of ENSEMBL gene IDs to find pathways for.
    }
    \item{endpoint}{
        The location of the SPARQL endpoint - default value is live Atlas endpoint
    }
}

\value{
A ranked list (pathway with greatest number of genes first) of pathwayresult classes with the following slots:
        \item{pathwayuri}{The URI of the corresponding pathway}
        \item{label}{The textual human readable name of the pathway}
    \item{numgenes}{A count of the number of genes from the gene list associated with this pathway}
    \item{genes}{Vector of genes from the gene list that are associated with this pathway}

}

\author{
    James Malone, Simon Jupp
}

\examples{
    ##get pathways which have the phrase 'ligand binding' included in their name
    genes <- c("ENSRNOG00000034254", "ENSRNOG00000009325", "ENSRNOG00000033065", "ENSRNOG00000015518", "ENSRNOG00000048955", "ENSRNOG00000031230")
    #pathwaylist <- getRankedPathwaysForGeneIds(genes)
}

