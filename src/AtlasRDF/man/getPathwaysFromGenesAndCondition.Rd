\name{getPathwaysFromGenesAndCondition}
\alias{getPathwaysFromGenesAndCondition}

\title{
Get pathways connected to genes which are differentially expressed for a user specified condition.
}
\description{
Get list of pathway names that are connected to genes (via proteins which are encoded for by those genes) and which are reported as
    differentially expressed under a specified condition. The condition is specified as an ontology URI. 
}
\usage{
getPathwaysFromGenesAndCondition(condition, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}
\arguments{
  \item{condition}{
URI of a condition based on the EFO ontology class. 
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
\references{
%% ~put references to the literature/web site here ~
}
\author{
James Malone, Simon Jupp
}

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

}

