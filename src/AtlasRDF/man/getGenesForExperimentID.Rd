\name{getGenesForExperimentID}
\alias{getGenesForExperimentID}

\title{
Get all of the genes reported in an experiment, speicifed by the Atlas experiment ID.
}
\description{
Get all of the gene as reported within an experiment. The experiment is specified as an Atlas experiment ID e.g. E-GEOD-3865
}
\usage{
getGenesForExperimentID(experiment, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}
\arguments{
  \item{experiment}{
    Atlas experiment ID.
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
###get genes for experiment ID E-MEXP-3739
#Can take several minutes to complete
#res <- getGenesForExperimentID("E-MEXP-3739")
###list first gene in the experiment
#res[1,]
}