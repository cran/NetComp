\name{netClass}
\alias{netClass}
\title{Network from Class List}
\usage{
  netClass(x, labels=NULL)
}
\arguments{
  \item{x}{Numeric vector containing the class/cluster assignments} 
  \item{labels}{Optional. Vector of labels for the individual observations} 
}
\description{
    Returns an unweighted adjacency matrix with a '1' where observations are from the same class,
	'0' if from different classes based on a vector of class assingments.
}
\details{
	Used for converting class/community assingments into an adjacency matrix.
}
\value{
    A square matrix is returned with the number of rows, columns equal to the length of \code{x} 
	and in the same order. If \code{lables} is supplied, they will be used as the row/column names.
}
\author{Shannon M. Bell}

\seealso{
  \code{\link[stats]{hclust}{cutree}}
  }
\examples{
 #using the USArrest dataset
 arrestCor<-cor(t(USArrests))
 tree<-cutree(hclust(as.dist(1-arrestCor), method='ward'), k=10)
 netClass(tree, labels=colnames(arrestCor))
 
}