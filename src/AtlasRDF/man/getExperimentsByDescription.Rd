\name{getExperimentsByDescription}
\alias{getExperimentsByDescription}

\title{
Search for experiments in Atlas based on a string, e.g. liver, cancer etc.
}
\description{
Find experiments in Atlas based on a search string such as a "liver" or "cancer". Returns the experiment IDs plus the text in which the search string was found.
}
\usage{
getExperimentsByDescription(searchterm, limit = 0, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}

\arguments{
  \item{searchterm}{
%%     ~~Describe \code{searchterm} here~~
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

