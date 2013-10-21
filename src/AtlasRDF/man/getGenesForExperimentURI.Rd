\name{getGenesForExperimentURI}
\alias{getGenesForExperimentURI}

\title{
Get all of the genes reported in an experiment, speicifed by the Atlas experiment URI.
}
\description{
Get all of the gene as reported within an experiment. The experiment is specified as an Atlas experiment URI e.g. <http://rdf.ebi.ac.uk/resource/atlas/E-GEOD-13396>
}
\usage{
getGenesForExperimentURI(experiment, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}
\arguments{
  \item{experiment}{
    Atlas experiment URI.
}
    \item{endpoint}{
        The location of the SPARQL endpoint - default value is live Atlas endpoint
    }

}

\value{
data frame containg following elements
  \item{expressionValue}{The summarised expression value for this gene}
 \item{dbXref}{The URI identifier for this gene}
 \item{pvalue}{The p-value representing differentially expression for this gene}
    \item{propertyValue}{The property value for which this gene was reported as differentially expressed}
}

\author{
James Malone, Simon Jupp
}

\examples{
###get genes for experiment <http://rdf.ebi.ac.uk/resource/atlas/E-GEOD-13396>

res <- getGenesForExperimentURI("<http://rdf.ebi.ac.uk/resource/atlas/E-MEXP-3739>")
###list first gene in the experiment
res[1,]
}
