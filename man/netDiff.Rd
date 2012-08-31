\name{netDiff}
\alias{netDiff}
\title{Network Difference}
\usage{
netDiff(matrix1, matrix2, cutoff=NULL, ...)
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
    Returns the difference of the input matrices. Output object is an
    adjacency matrix where edges present in \code{matrix1} but not \code{matrix2}
    are returned. The edge value is maintained provided \code{cutoff} = NULL.
    Edges in both graphs failing to meet the \code{cutoff}, if provided, are set
    to zero before taking the graph difference.
}
\details{
    Matrices must be square and have row and column labels. Output adjacency
    matrix can be used directly for creating a graph object.
    }
\value{
    \code{netDiff} returns an adjacency matrix containing edges present in
    \code{matrix1} that are not present in \code{matrix2}. Edges below
    \code{cutoff} are set to zero.
}
\author{Shannon M. Bell}
\examples{
 #using the state.x77 and USArrest datasets
 #remove data from states for illustration
 ssArrest<-subset(t(USArrests), select=-c(Alabama,Colorado,Delaware))
 ssState<-subset(t(state.x77), select=-c(Alabama, Arizona, Iowa))
 arrestCor<-cor(ssArrest)
 stateCor<-cor(ssState)
 dataDiff<-netDiff(stateCor, arrestCor)
 dataDiff[1:15,1:5]
 #Setting a cutoff to remove any edges that are below 0.6 
 dataDiff.6<-netDiff(stateCor, arrestCor, cutoff=0.6)
 dataDiff.6[1:15,1:5]
}