library(ROntoTools)
library(graph)
test_edgeData_keggPathwaygraphs <- function()
{

    kpg<-keggPathwayGraphs("hsa")
    kpgN<-keggPathwayNames(organism="hsa" )
    ed<-edgeData_keggPathwaygraphs("hsa:7029","hsa:898")
    id<-names(ed)

    KeggGraph2Matrix <- function(g)
    {
    ind <- which(as.vector(as(g, "matrix")) != 0)
    source <- (ind-1) %% length(nodes(g)) + 1
    destination <- (ind-1) %/% length(nodes(g)) + 1
    return(cbind(nodes(g)[source], nodes(g)[destination]))
    }

    mat <- KeggGraph2Matrix(kpg[[id]])
    mat <- cbind(mat , unlist(edgeData(kpg[["path:hsa04110"]], mat[,1], mat[,2], "subtype")))
    allEdges <- cbind(rownames(mat), mat)
    checkTrue("hsa:7029|hsa:898" %in% allEdges[,1])


}


test_nodeData_keggPathwaygraphs <- function()
{
    kpg<-keggPathwayGraphs("hsa")
    kpgN<-keggPathwayNames(organism="hsa" )
    nd<-nodeData_keggPathwaygraphs("hsa:7029")
    id<-names(nd)[1]
    checkTrue("hsa:7029" %in% nodes(kpg[["path:hsa04110"]]))

}
