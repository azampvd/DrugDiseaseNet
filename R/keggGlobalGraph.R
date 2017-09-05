
#' keggGlobalGraph: Construct global graph of KEGG pathways.
#'
#'
#' @param updateKEGG re-download KEGG data
#' @param organism organism code as defined by KEGG
#'
#' @return
#'
#' An object of class \code{\link{graphNEL}} encoding the KEGG
#' global graph information.
#'
#' @author
#'
#' Azam Peyvandipour and Sorin Draghici
#'
#' @references
#' #' Draghici S., et.al.: "A systems biology approach for
#' pathway level analysis". Genome Research, 17, 2007.
#'
#' @seealso \code{\link{keggPathwayGraphs}},
#' \code{\link{keggPathwayNames}}
#'
#' @examples
#' library(graph)
#' gg<-keggGlobalGraph()
#' class(gg)
#' head(nodes(gg))
#' class(edges(gg))
#' allEdges<-edges(gg)
#' allEdges$`hsa:2065`
#' @import ROntoTools
#' @import graph
#' @export
keggGlobalGraph<-function(organism="hsa",updateKEGG=FALSE){

    if(updateKEGG)
        kpg<-keggPathwayGraphs("hsa",updateCache=updateKEGG,verbose=FALSE)
    else
        {
        load(paste(system.file("extdata",package="DrugDiseaseNet"),
        "/globalGraph.RData", sep= ""))
        return(globalGraph)
        }

    kpg <- setEdgeWeights(kpg, edgeTypeAttr= "subtype",
    edgeWeightByType=list(activation= 1, inhibition= -1,
    expression=1, repression=-1),
    defaultWeight=0)

    KeggGraph2Matrix <- function(g)
    {
        ind <- which(as.vector(as(g, "matrix")) != 0)
        source <- (ind-1) %% length(nodes(g)) + 1
        destination <- (ind-1) %/% length(nodes(g)) + 1
        return(cbind(nodes(g)[source], nodes(g)[destination]))
    }
    allNodes <- unique(unlist(sapply(kpg, nodes)))
    allEdges <- sapply(kpg, function(g)
    {
    mat <- KeggGraph2Matrix(g)
    mat <- cbind(mat , unlist(edgeData(g, mat[,1], mat[,2], "subtype")))
    cbind(rownames(mat), mat)
    })
    allEdges <- do.call(rbind, allEdges)

    eT <- tapply(allEdges[,4], as.factor(allEdges[,1]), function(x) x )
    eT <- sapply(eT, function(x) paste(unique(unlist(strsplit(unlist(x)
    , ','))), collapse=','))
    sdMat <- do.call(rbind, strsplit(names(eT), '\\|'))
    gg <- ftM2graphNEL(sdMat, V=allNodes)
    edgeDataDefaults(gg, "subtype") <- NA
    edgeData(gg, sdMat[,1], sdMat[,2], "subtype") <- eT
    globalGraphL <- setEdgeWeights(list(gg))
    globalGraph<-globalGraphL[[1]]
    globalGraph<-updateGlobalGraph(globalGraph)
    return(globalGraph)
    }


#' @import graph
calculate.B <- function(g, non.zero=TRUE){

    if (non.zero)
    # Nds: number of downstream genes
    nds <- sapply(edgeWeights(g), function(x) sum(x != 0 ))
    else
    nds <- sapply(edges(g), length)

    # add 1 for all genes with no downstream genes
    nds[nds == 0] <- 1

    # compute B = (I - beta/nds)
    B <- diag(length(nodes(g))) - t(as(g, "matrix")) /
    matrix(nds, byrow=TRUE, nrow=length(nds), ncol=length(nds))
    return(B)}


#' keggGlobalGraph: Construct global graph of KEGG signaling pathways.
#'
#'
#' @param globalGraph KEGG global graph as object of type
#' \code{\link{graphNEL}}
#'
#' @param edgesToremove An object of class \code{\link{list}}
#' encoding edges need to be removed from the graph
#'
#' @details
#'
#' See cited document for more details.
#'
#' @return
#'
#' An object of class \code{\link{graphNEL}} encoding the
#' KEGG global graph information.
#'
#' @author
#'
#' Azam Peyvandipour and Sorin Draghici
#'
#' @references
#' #' Draghici S., et.al.: "A systems biology approach for
#' pathway level analysis". Genome Research, 17, 2007.
#'
#' @seealso \code{\link{keggPathwayGraphs}}, \code{\link{keggGlobalGraph}}
#' @examples
#' library(graph)
#' gg<-keggGlobalGraph()
#' class(edges(gg))
#' length(unlist(edges(gg)))
#' allEdges<-edges(gg)
#' allEdges$`hsa:2065`
#' ggUpdated<-updateGlobalGraph(gg,edgesToremove=
#' list(c("hsa:2065","hsa:2549"),c("hsa:2065","hsa:25759")))
#' length(unlist(edges(ggUpdated)))
#' @import graph
#' @export
    updateGlobalGraph <- function(globalGraph,edgesToremove=
    list(c("hsa:6657","hsa:79923"))){

    #message("Global graph must not have any cycle.")
    adjmatrix <- as(globalGraph, "matrix")
    #M <-calculate.B(globalGraph)
    #abs(det(M)) #if this is zero we need to remove
    #the nodes causing self cycles
    index<-which(diag(adjmatrix) != 0)

    #remove the nodes with self cycle
    if (length(index)>0){
        globalGraph <- removeNode(colnames(adjmatrix)[index],globalGraph)
        index <- NULL
    }
    #cycles need to be identified and edges causing the cycles need
    #to be removed from the graph
    #in this global graph, edge (hsa:6657','hsa:79923') was removed

    #check for edges causing loops
    #edge ("hsa:6657","hsa:79923") causes a loop- removed

    if(!is.null(edgesToremove))
        for (i in 1:length(edgesToremove))
            globalGraph <- removeEdge(edgesToremove[[i]][1],
            edgesToremove[[i]][2],globalGraph)

    #M <- calculate.B(globalGraph)
    #(det(M)) #determinant of M is not zero
    # if(det(M) == 0)
    # stop("Cycles need to be removed!")
    return (globalGraph)}

