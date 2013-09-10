\name{getGeneUriFromName}
\alias{getGeneUriFromName}

\title{
    Get URI of a gene based on the common gene name.
}
\description{
   Get the URI of an ENSEMBL gene based on a search by the common gene name.
}
\usage{
    getGeneUriFromName(genename, taxon, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
    \item{genename}{
        %%     ~~Describe \code{genename} here~~
    }
    \item{taxon}{
        %%     ~~Describe \code{taxon} here~~
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

