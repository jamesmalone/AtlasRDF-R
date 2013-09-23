\name{getOntologyMappings}
\alias{getOntologyMappings}

\title{
    Get mappings for a given ontology class URI to EFO using the NCBO BioPortal
}
\description{
    Get mappings for a given class URI from an ontology other than EFO using the NCBO BioPortal. Function is useful for querying the Atlas RDF using your own ontology. For instance, ICD-10 or SnoMed.
    The function will find mappings between these classes in other ontologies do EFO where they exist. The returned EFO class can be used to query Atlas.
}
\usage{
    getOntologyMappings(searchuri, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
    \item{searchuri}{
        The ontology class URI for which mappings to EFO are required.
    }
    \item{endpoint}{
        SPARQL endpoint to query data from. Default value is the live Atlas SPARQL endpoint but this can be overidden.     
    }
}

\value{
List of EFO URI(s) to which the class maps
}
\references{
NCBO BioPortal - http://bioportal.bioontology.org/
}
\author{
James Malone, Simon Jupp
}

\examples{
###find mappings in EFO from  <http://purl.bioontology.org/ontology/SNOMEDCT/87163000>  (leukemia in snomed)
getOntologyMappings("<http://purl.bioontology.org/ontology/SNOMEDCT/87163000>")
  
}
