\name{doFishersEnrichment}
\alias{doFishersEnrichment}

\title{
Perform a gene set enrichment using gene list based on gene URIs. 
}
\description{
Perform an enrichment of an input gene list using the Atlas as a background set for a specified species. The genes are given as a set of gene URIs (e.g. <http://identifiers.org/ensembl/ENSG00000229807>). The results are a list of experimental factors for which these genes are enriched.
    Access the backgroun data sets for enrichment from https://github.com/jamesmalone/AtlasRDF-R
}
\usage{
doFishersEnrichment(genelist, genelist_bg, genecounts)
}
\arguments{
\item{genelist}{
     Variable \code{genelist} a vector of genes using gene URIs
}
\item{genelist_bg}{
     \code{genelist_bg} background list appropriate for the species of interest should be obtained from svn
}
\item{genecounts}{
      \code{genecounts}  counts list appropriate for the species of interest should be obtained from svn
}
}

\value{
Returns a list of enrichmentresult objects which are described as follows:

\item{factoruri}{the URI of the experimental factor}
\item{label }{the human readable label for the experimental factor}
\item{p.value}{the p.value for the Fisher's exact test for the gene list with this factor. This can be seen as avidence for this factor being highly enriched (low p value) or not (high p value) for the input gene list}
\item{estimate}{an estimate of the odds ratio. Note that the conditional Maximum Likelihood Estimate (MLE) rather than the unconditional MLE (the sample odds ratio) is used. Only present in the 2 by 2 case.}
\item{alternative}{a character string describing the alternative hypothesis.}
\item{null.value}{the odds ratio under the null, or. Only present in the 2 by 2 case.}
\item{method}{the character string "Fisher's Exact Test for Count Data".}
\item{enrichedgenes}{list of the genes from the input gene list that were enriched for this factor as a list of Gene URIs. To get the gene names use the function \code{getClassLabel}}
}
\references{
Agresti, A. (1990) Categorical data analysis. New York: Wiley. Pages 59-66. 
}
\author{
James Malone, European Bioinformatics Institute
Simon Jupp, European Bioinformatics Institute
Maryam Soleimani, European Bioinformatics Institute
}

\examples{

###perform gene set enrichment using Atlas as background to find factors genes are enriched for
#load the background files for the species of interest (requires downloading)
#load("human_gene_list.RData")   #human_genelist_bg
#load("human_factor_counts.RData")    #human_factor_counts

###creat or load your genelist of interest
genelist <- c("<http://identifiers.org/ensembl/ENSG00000184674>", 
"<http://identifiers.org/ensembl/ENSG00000138193>", 
"<http://identifiers.org/ensembl/ENSG00000229807>", 
"<http://identifiers.org/ensembl/ENSG00000106546>", 
"<http://identifiers.org/ensembl/ENSG00000163430>", "<http://identifiers.org/ensembl/ENSG00000145817>", "<http://identifiers.org/ensembl/ENSG00000003756>",
"<http://identifiers.org/ensembl/ENSG00000177731>")

###do enrichment
#results <- doFishersEnrichment(genelist, human_genelist_bg, human_factor_counts)

}


