\name{getSpeciesSpecificEnsemblGenesForExFactor}
\alias{getSpeciesSpecificEnsemblGenesForExFactor}

\title{
Get ENSEMBL genes which have been differentially expressed for a given condition and species. 
}
\description{
Get the ENSEMBL gene IDs and names which are reported as differentially expressed in Atlas for a specified conditions and species.
}
\usage{
getSpeciesSpecificEnsemblGenesForExFactor(exfactor, taxon, limit = 0, endpoint="http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
 \item{exfactor}{
%%     ~~Describe \code{exfactor} here~~
}
 \item{taxon}{
%%     ~~Describe \code{taxon} here~~
}
  \item{limit}{
%%     ~~Describe \code{limit} here~~
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
