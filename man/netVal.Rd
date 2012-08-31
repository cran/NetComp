\name{netVal}
\alias{netVal}
\title{Network Validation}
\usage{
  netVal(matrix1, matrix2, method='ward', k=200, h=NULL,...)
}
\arguments{
  \item{matrix1}{Adjacency matrix containing the true topology} 
  \item{matrix2}{Adjacency matrix containing the test topology} 
  \item{method}{Hierarchical clustering method to be passed to \code{\link{hclust}}.}
  \item{k}{Number of groups. Passed to \code{\link{cutree}}}
  \item{h}{Heights where tree is to be cut. Passed to \code{\link{cutree}}}
  \item{...}{Other parameters.}
}
\description{
    Returns vector with the True Positives, True Negatives, False Postives,
    False Negatives, the Adjusted Rand Index, and the Balanced Accuracy of the
    test network (\code{matrix2}) compared to the true network (\code{matrix1}).
    Network communities are based on hierachical clustering, with communities
    defind by \code{k} or \code{h}.
}
\details{
    The input matrices must be square but do not need to contain the same edges.
}
\value{
    \code{netVal} returns an vector with values True Positives, True Negatives, False Postives,
    False Negatives, Adjusted Rand Index, and Balanced Accuracy.
}
\author{Shannon M. Bell}
  \seealso{
  \code{\link[stats]{hclust}{cutree}}
  }
\examples{
 #using the state.x77 and USArrest datasets
 arrestCor<-cor(t(USArrests))
 stateCor<-cor(t(state.x77))
 netVal(stateCor, arrestCor, method='ward', k=10)
 netVal(stateCor, stateCor, method='ward', k=10)
}