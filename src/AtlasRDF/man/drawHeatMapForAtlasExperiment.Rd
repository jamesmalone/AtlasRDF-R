\name{drawHeatMapForAtlasExperiment}
\alias{drawHeatMapForAtlasExperiment}

\title{
Produce a heatmap of differentially expressed genes vs conditions for an Atlas experiment.
}
\description{
Draws a heatmap of genes that are differentially expressed beyond a given t stat signficance against conditions they are diff expressed for for a given Atlas experiment. 
}
\usage{
drawHeatMapForAtlasExperiment(experimentid, tstatsignificance, endpoint)
}
\arguments{
  \item{experimentid}{
Atlas experiment ID which is to be visualized
}
  \item{tstatsignificance}{
The significance beyond which a gene should be included for visualisation. A significance value of '5' will only include genes with a t-stat of >= 5 or <= -5.
}
  \item{endpoint}{
SPARQL endpoint to query data from. Default value is the live Atlas SPARQL endpoint but this can be overidden. 
}

}

\value{
  matrix of gene names with t-statistics for level of differential expression for each factor (one factor per column)

}

\author{
Simon Jupp, James Malone
}

\examples{

#data <- drawHeatMapForAtlasExperiment("E-MTAB-1670", 12)

}
