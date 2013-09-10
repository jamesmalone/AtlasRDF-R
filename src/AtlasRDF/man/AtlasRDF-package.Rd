\name{AtlasRDF-package}
\alias{AtlasRDF-package}
\alias{AtlasRDF}
\docType{package}
\title{
    Gene Expression Atlas query and gene set enrichment package.
}
\description{
    Query the Gene Expression Atlas RDF data at the European Bioinformatics Institute using genes, experimental factors (such as disease, cell type, compound treatments), pathways and proteins. Also contains a function to perform an enrichment of your gene list across Experimental Factor Ontology (EFO) using the Atlas background set. }
\details{
    \tabular{ll}{
        Package: \tab AtlasRDF\cr
        Type: \tab Package\cr
        Version: \tab 1.0\cr
        Date: \tab 2013-09-09\cr
        License: \tab Apache 2.0\cr
    }
    
}
\author{
    James Malone, Simon Jupp, Maryam Soleimani
    
    Maintainer: James Malone <malone@ebi.ac.uk>
}
\references{
    Gene Expression Atlas RDF - http://wwwdev.ebi.ac.uk/rdf/services/atlas/
    Gene Expression Atlas - http://www.ebi.ac.uk/gxa
    Experimental Factor Ontology - http://www.ebi.ac.uk/efo
    Kapushesky M, Adamusiak T, Burdett T, Culhane A, Farne A, Filippov A, Holloway E, Klebanov A, Kryvych N, Kurbatova N, Kurnosov P, <b>Malone J</b>, Melnichuk O, Petryszak R, Pultsin N, Rustici G, Tikhonov A, Travillian RS, Williams E, Zorin A, Parkinson H, Brazma A. (2012) Gene Expression Atlas update--a value-added database of microarray and sequencing-based functional genomics experiments. Nucleic Acids Res. 2012 Jan;40(Database issue):D1077-81.
}
\keyword{ package }

\examples{
    ###find experiments where the description contains a given search term
    #experimentlist <- getExperimentsByDescription("breast carcinoma")
    
    ###get gene list for given experiment ID
    #genelist <- getGenesForExperiment("E-GEOD-1561")
    
    ###get list of pathways associated with a gene for which the genes are differentially expressed for a given condition (Type I diabetes in this case)
    #pathwayinfo <- getPathwaysFromGenesAndCondition("efo:EFO_0001359")
    #pathwayinfo$pathwayname
    
    ###perform gene set enrichment using Atlas as background to find factors genes are enriched for
    #load the background files for the species of interest (requires downloading)
    #load("human_bg/human_gene_list.RData")   #human_genelist_bg
    #load("human_bg/human_factor_counts.RData")    #human_factor_counts
    
    ###creat or load your genelist of interest
    #genelist <- c("GSTT1", "PLCE1", "XIST", "ENOSF-1", "AHR", "FSTL1", "YIPF5", "RBM5", "TNP01", "CREB1", "NOX4", "FLII")
    
    ###get id for species 'human'
    #taxon <-getTaxonURI("human")  #obo:NCBITaxon_9606
    
    ###do enrichment
    #results <- doFishersEnrichmentForGeneNames(genelist, taxon, human_genelist_bg, human_factor_counts)
    
    ###visualize results, filtering by p-value
    #vizPvalues(results, "0.005")
}
