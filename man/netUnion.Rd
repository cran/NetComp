\name{netUnion}
\alias{netUnion}
\title{Network Union}
\usage{
netUnion(matrix1, matrix2, cutoff=NULL, ...)
}
\arguments{
  \item{matrix1}{Square matrix (e.g. correlation or adjacency) containing
    row/column labels} 
  \item{matrix2}{Square matrix (e.g. correlation or adjacency) containing
    row/column labels}
  \item{cutoff}{The cutoff value. Edges less than this value (absolute value
    considered)are converted to zero. }
  \item{...}{Other parameters.}
}
\description{
    Returns the union of the input matrices. Output object is an
    adjacency matrix where the edge weight is 1 if an edge is present in both
    graphs, 0.5 if present in 1 graph, and 0 if present in neither graph. Edges
    in both graphs failing to meet the \code{cutoff}, if provided, are set to
    zero before performing the graph union.
}
\details{
    Matrices must be square and have row and column labels. Output adjacency
    matrix can be used directly for creating a graph object.
    }
\value{
    \code{netUnion} returns an adjacency matrix containing edges present either
    both graphs.
}
\author{Shannon M. Bell}
\references{
  Bell SM, Last RL. SimMeasure: A non-imputing approach to analyzing missing data.
  Bioinformatics, in review
  }
\examples{
 #using the state.x77 and USArrest datasets
 #remove data from states for illustration
 ssArrest<-subset(t(USArrests), select=-c(Alabama,Colorado,Delaware))
 ssState<-subset(t(state.x77), select=-c(Alabama, Arizona, Iowa))
 arrestCor<-cor(ssArrest)
 stateCor<-cor(ssState)
 dataUnion<-netUnion(stateCor, arrestCor)
 dataUnion[1:15,1:5]
 #Setting a cutoff to remove any edges that are below 0.6 
 dataUnion.6<-netUnion(stateCor, arrestCor, cutoff=0.6)
 dataUnion.6[1:15,1:5]
}