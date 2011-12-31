#################################
##Name: shortPaths
##Description: Function to calculate the shortes path lenght between 2 nodes in a graph object
##O/S: for R
##Date: 12/9/2011
##Author: Shannon M. Bell
##Company: Michigan State University
##note the path length might not be the only path but it is the shortest
###################################

#library(igraph)
#this function will find the shortest paths from the input graph and return a matrix with the vertex names
#as the row/column names and the value being path length.
#NOTE! it is important to have your full graph, ie the graph made from the adjacency matrix
#if you want to match it back to the adjacency matrix because some nodes may be removed
#prelim tests suggest that you still get a fine graph output though
#nodes that are unconnected will return a -1 and throw a warning
shortPaths<-function(graph, ...){
    l<-length(V(graph))
    edgeNames<-unlist(V(graph)$name)
    dataMat<-NULL
    for(i in 1:l){
        temp<-get.shortest.paths(graph, from=(i-1), to=c(0:(l-1)), mode='all', weights=NA)
        paths<-unlist(lapply(temp, function(x) length(x)-1))
        dataMat<-cbind(dataMat, paths)
    }
    dataMat<-as.data.frame(dataMat)
    colnames(dataMat)<-edgeNames
    row.names(dataMat)<-edgeNames
    dataMat
}
