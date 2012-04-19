\name{matrix_threshold}
\alias{matrix_threshold}
\title{Remove values from Matrix based on threshold}
\usage{
matrix_threshold(matrix, threshold=NULL, minval=0, maxval=NULL, abs=TRUE, rmna=FALSE,...)
}
\arguments{
  \item{matrix}{matrix } 
  \item{threshold}{Integer}
  \item{minval}{Integer. Replaces values below threshold. }
  \item{maxval}{Integer. Replaces values above or equal to \code{threshold}. }
  \item{abs}{Logical. If TRUE the absolute value will be considered. }
  \item{rmna}{Logical. If TRUE NA's will be replaced with \code{minval}. }
  \item{...}{Other parameters.}
}
\description{
    Returns matrix with values above or  below the \code{threshold} are replaced.
    Can be used for converting to a binary matrix or for thresholding an adjacency
    matrix.
}
\details{
    If \code{rm.na} is TRUE and \code{minval} is NULL, then NA's will be replaced
    with 0.
    }
\value{
    \code{max_threshold} returns an matrix with values removed according to parameters.
}
\author{Shannon M. Bell}
\references{
  Bell SM, Last RL. SimMeasure: A non-imputing approach to analyzing missing data.
  Bioinformatics, in review
  }
\examples{
 #using USArrest dataset
 arrestCor<-cor(t(USArrests))
 arrestCor.7<-matrix_threshold(arrestCor, threshold=0.7, minval=0, abs=TRUE, rm.na=FALSE)
 arrestCor[10:15,1:5]
 arrestCor.7[10:15,1:5]
}