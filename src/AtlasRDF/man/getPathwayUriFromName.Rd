\name{getPathwayUriFromName}
\alias{getPathwayUriFromName}

\title{
    Get pathway URI given a pathway name
}
\description{
    Given the name of a pathway find URI(s) that match this name
}
\usage{
    getPathwayUriFromName(name, endpoint = "http://www.ebi.ac.uk/rdf/services/atlas/sparql")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{name}{
       The textual name of the pathway to search for  
    }
    \item{endpoint}{
        The location of the SPARQL endpoint - default value is live Atlas endpoint
    }
}

\value{
data frame containing slots
    \item{pathwayuri}{The URI of the corresponding pathway}
    \item{label}{The textual human readable name of the pathway}

}

\author{
    James Malone, Simon Jupp
}

\examples{
    ##get pathways which have the phrase 'ligand binding' included in their name
    getPathwayUriFromName("ligand binding")
}

