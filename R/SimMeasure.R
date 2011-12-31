SimMeasure<-function(data, threshold=NULL, ...){
    .Call("SimMeasure",data, threshold, pkg="NetComp" )
}