#################################
##Name: netDiff
##Description: Function to find the difference b/w 2 square matricies (comming out of correlation)
##O/S: for R
##Date: 12/9/2010
##Author: Shannon M. Bell
##Company: Michigan State University
###################################
#This code will take as input (required) 2 square matrixies that must have row/column names
#MUST be a square matrix!!!
#if you are interested in doing a subset, take those nodes from the corr matrix (still needs to be square)
#the output is a graph with nodes/edges that are in matrix 1, not in matrix2
#the values of the edges b/w nodes found in matrix 1 will be returned as is
#provided (if cuttoff !=NULL) they are above the abs(cuttoff)
#output is a matrix of the same dimension (and order) as matrix 1 with 1 in spots where
#edge is unique to matrix 1 and above cutoff (if provided), zero otherwise


netDiff<-function(matrix1, matrix2, cutoff=NULL, ...){
    #ordering the graph
    matrix1<-matrix1[sort(row.names(matrix1)), sort(colnames(matrix1))]
    matrix2<-matrix2[sort(row.names(matrix2)), sort(colnames(matrix2))]
    #this part transforms mat1 and mat2 into a 1/0 graph
    #save the origional data for getting weights later on
    if(is.null(cutoff)){
        matrix1[abs(matrix1) <=0 | is.na(matrix1)] <-0
        ori.m1<-matrix1
        matrix1[abs(matrix1) >0] <-1
        matrix2[abs(matrix2) <=0 | is.na(matrix2)] <-0
        matrix2[abs(matrix2) >0] <-1
    }
    if(!is.null(cutoff)){
        matrix1[abs(matrix1) <= cutoff | is.na(matrix1)] <-0
        ori.m1<-matrix1
        matrix1[abs(matrix1) > cutoff] <-1
        matrix2[abs(matrix2) <= cutoff | is.na(matrix2)] <-0
        matrix2[abs(matrix2) > cutoff] <-1
    }
    
    #need to put in error for in matrix1 or 2 is not square or is null
    shared.names<-intersect(colnames(matrix1), colnames(matrix2))
    #this is the total names, and the length of the final product
    total.names<-sort(union(colnames(matrix1), colnames(matrix2)))
    #need to get what is in first matrix and NOT in second
    missing2.names<-setdiff(colnames(matrix1), colnames(matrix2))
    missing1.names<-setdiff(colnames(matrix2), colnames(matrix1))
    #get the number of smaples
    l1<-length(colnames(matrix1))
    l2<-length(colnames(matrix2))
    #this goes through and addes rows/col of 0 for missing in matrix2
    if(length(missing2.names) >0){
        temp<-NULL
        for(i in 1:length(missing2.names)){
            temp<-cbind(temp, rep(0,l2))
        }
        temp<-as.data.frame(temp)
        colnames(temp)<-missing2.names
        m2.temp<-cbind(matrix2, temp)
        temp2<-NULL
        for(i in 1:length(missing2.names)){
            temp2<-rbind(temp2, rep(0,length(total.names)))
        }
        temp2<-as.data.frame(temp2)
        colnames(temp2)<-colnames(m2.temp)
        rownames(temp2)<-missing2.names
        m2<-rbind(m2.temp, temp2)
    } else{
        m2<-matrix2
    }
    #this goes through and addes rows/col of 0 for missing in matrix1
    if(length(missing1.names) >0){
        temp<-NULL
        for(i in 1:length(missing1.names)){
            temp<-cbind(temp, rep(0,l1))
        }
        temp<-as.data.frame(temp)
        colnames(temp)<-missing1.names
        m1.temp<-cbind(matrix1, temp)
        temp2<-NULL
        for(i in 1:length(missing1.names)){
            temp2<-rbind(temp2, rep(0,length(total.names)))
        }
        temp2<-as.data.frame(temp2)
        colnames(temp2)<-colnames(m1.temp)
        rownames(temp2)<-missing1.names
        m1<-rbind(m1.temp, temp2)
    } else{
        m1<-matrix1
    }
    matrix1b<-m1[total.names, total.names]
    matrix2b<-m2[total.names, total.names]
    #now i want to subtract graph 2 from graph 1
    #result will be negative if in mat 2 but not mat 1, 0 if in both, 1 if only in mat 1
    temp.mat<-matrix1b - matrix2b
    mat3<-temp.mat[colnames(ori.m1), colnames(ori.m1)]
    matrix3 <-mat3*ori.m1
    matrix3
}
