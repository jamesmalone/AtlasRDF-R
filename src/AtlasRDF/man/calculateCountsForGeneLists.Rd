\name{calculateCountsForGeneLists}
\alias{calculateCountsForGeneLists}
\title{
    Internal function for counting genes in a gene list against the background
}
\description{
    Internal function for counting genes in a gene list against the background
}
\usage{
    calculateCountsForGeneLists(genelist, genelist_bg, genecounts)
}

\arguments{
    \item{genelist}{
        Vector of gene uris which are being analysed for enrichment purposes and for which the factor-gene counts will be produced
    }
    \item{genelist_bg}{
        Background list for species of interest
    }
    \item{genecounts}{
        Factor counts for species of interest
    }
}
\value{
    returns a hash of factorbackground objects which contain information about which genes are diff expressed for a given factor
}
\author{
    James Malone, Simon Jupp
}

\examples{
    
}
\keyword{internal}

