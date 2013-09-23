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
  The search term.
}
  \item{limit}{
 Optional limit on number of returned results (default is everything)
}
  \item{endpoint}{
    The location of the SPARQL endpoint - default value is live Atlas endpoint
}
 
}

\value{

 \item{$experiment}{Atlas Experiment ID of matched experiment to search term'}
 \item{$description}{Description of experiment which matched term}

}

\author{
James Malone, Simon Jupp
}

\examples{
###get experiments which contain the word liver in the description
getExperimentsByDescription("liver")

}

