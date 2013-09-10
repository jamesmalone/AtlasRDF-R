\name{vizPvalues}
\alias{vizPvalues}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize the results of an enrichment filtering for a specificed p-value.
}
\description{
Produce a plot of the most enriched experimental factors following an Atlas gene set enrichment. The cutoff allows a filter to be placed by p-value.
}
\usage{
vizPvalues(resultset, cutoff = "0.05")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{resultset}{
        A set of results following an AtlasRDF gene set enrichment. 
}
  \item{cutoff}{ 
        The p-value threshold at which to include a gene for visualization
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
    Plot of results
}
\author{
James Malone
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--    or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
