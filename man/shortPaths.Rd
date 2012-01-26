\name{shortPaths}
\alias{shortPaths}
\title{Find the Shortest Paths between Nodes in a Graph}
\usage{
shortPaths(graph, ...)
}
\arguments{
  \item{graph}{Graph object created using }
  \item{...}{Other parameters.}
}
\description{
    \code{shortPaths} is a wrapper function for \code{\link{get.shortest.paths}}
    that calculated the shortest path length between all nodes in a graph
    object. Additionally, if the graph has node labels those are retained. The
    full graph (ie from the adjacency matrix) is recommended rather than
    removing unconnected nodes. This is an alternative method to
    \code{\link{shortest.paths}}.
}
\details{
    The output matrix will contain the shortest path between each note pair. If
    the pair is unconnected, a warning will be returned and -1 value will be
    present for that pair.
    }
\value{
    \code{shortPaths} returns a matrix with rows/columns corresponding to the
    edges and the values the shortest path length.
}
\author{Shannon M. Bell}
\references{
  Bell SM, Last RL. NetComp: Large Scale Data to Graph Analysis.
  Bioinformatics, in review
  }
\seealso{
  \code{\link{igraph}}, \code{\link{get.shortest.paths}},
  \code{\link{shortest.paths}}.
}
\examples{
 #using the state.x77 and USArrests datasets
 library(igraph)
 arrestCor<-cor(t(USArrests))
 #to remove low weighted edges
 arrestCor[abs(arrestCor) <0.5 | is.na(arrestCor)]<-0
 g<-graph.adjacency(as.matrix(arrestCor), weighted=TRUE, mode='upper', diag=FALSE)
 paths<-shortPaths(g)
 paths[1:15,1:5]
}