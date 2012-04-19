\name{netIntersect}
\alias{netIntersect}
\title{Network Intersection}
\usage{
netIntersect(matrix1, matrix2, cutoff=NULL, absolute=TRUE, ...)
}
\arguments{
  \item{matrix1}{Square matrix (e.g. correlation or adjacency) containing
    row/column labels} 
  \item{matrix2}{Square matrix (e.g. correlation or adjacency) containing
    row/column labels}
  \item{cutoff}{The cutoff value. Edges less than this value (absolute value
    considered)are converted to zero. }
  \item{absolute}{Logical, should the absolute values be considered for the
    cutoff value.}
  \item{...}{Other parameters.}
}
\description{
    Returns the intersection of the input matrices. Output object is an
    adjacency matrix where the edge weight is the average value of the absolute
    value of edges in BOTH matrices. Edges not present in either, or below the
    \code{cutoff} are removed or set to zero.
}
\details{
    Matrices must be square and have row and column labels. Output adjacency
    matrix can be used directly for creating a graph object.
    }
\value{
    \code{netIntersect} returns an adjacency matrix containing edged present in
    both graphs.
}
\author{Shannon M. Bell}
\references{
  Bell SM, Last RL. SimMeasure: A non-imputing approach to analyzing missing data.
  Bioinformatics, in review
  }
\examples{
 #using the state.x77 and USArrests datasets
 arrestCor<-cor(t(USArrests))
 stateCor<-cor(t(state.x77))
 isect<-netIntersect(stateCor, arrestCor)
 isect[1:15,1:5]
 #Setting a cutoff to remove any edges that are below 0.6 
 isect.6<-netIntersect(stateCor, arrestCor, cutoff=0.6)
 isect.6[1:15,1:5]
}