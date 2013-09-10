\name{getGenesForExperiment}
\alias{getGenesForExperiment}

\title{
Get all of the genes reported in an experiment, speicifed by the Atlas experiment ID.
}
\description{
Get all of the gene as reported within an experiment. The experiment is specified as an Atlas experiment ID e.g. E-GEOD-3865
}
\usage{
getGenesForExperiment(experiment, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
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
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}

\author{
James Malone, Simon Jupp
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

}
