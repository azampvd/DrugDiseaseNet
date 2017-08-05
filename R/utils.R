#' utils: utility functions
#' names in which there is an edge between the given nodes
#'
#' @param updateKEGG re-download KEGG data
#' @param edge_from A character of node name
#' @param edge_to A character of node name
#'
#' @details
#'
#' #' Draghici S., et.al.: "A systems biology approach for
#' pathway level analysis". Genome Research, 17, 2007.
#'
#' @return
#'
#' A character vector of KEGG pathway names in which
#' edge (from_node, to_node) exists
#'
#' @author
#'
#' Azam Peyvandipour and Sorin Draghici
#'
#'
#' @seealso \code{\link{keggPathwayGraphs},\link{keggPathwayNames}}
#'
#' @examples
#'pathways<-edgeData_keggPathwaygraphs(edge_from="hsa:208",
#'edge_to="hsa:2309",updateKEGG=FALSE)
#' @import graph
#' @import ROntoTools
#' @export
edgeData_keggPathwaygraphs<-function(edge_from,edge_to,updateKEGG=FALSE)
    {


    kpg<-keggPathwayGraphs("hsa",updateCache = updateKEGG)
    kpgN<-keggPathwayNames(organism ="hsa" ,updateCache =updateKEGG )
    kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype",
                        edgeWeightByType = list(activation = 1, inhibition = -1,
                        expression = 1, repression = -1),
                        defaultWeight = 0)

    kpg_matrix<-lapply(kpg,function(x){
        m<-as(x,"matrix")
        return(m)
        })

    edgeExist<-unlist(sapply(kpg_matrix,function(x)
        {
        if(all(c(edge_from,edge_to) %in% rownames(x)))
            {
            if(x[edge_from,edge_to]!=0)
                return(TRUE)
            }
        else
        return(FALSE)
        }))

    res<-kpgN[names(edgeExist)[which(edgeExist)]]
    if(!length(res))
        stop(paste("The edge (",edge_from," , " ,edge_to,")
            does not belong to KEGG pathways",sep=""))
        else
        return(res)
}


#' edgeData_keggPathwaygraphs: A character vector of
#' KEGG pathway names in which the given node exists
#'
#' @param updateKEGG re-download KEGG data
#' @param node A character of node name
#'
#' @details
#'
#' See cited document for more details.
#'
#' @return
#'
#' A character vector of KEGG pathway names in which the given node exists
#'
#' @author
#'
#' Azam Peyvandipour and Sorin Draghici
#'
#'
#' @seealso \code{\link{keggPathwayGraphs},\link{keggPathwayNames}}
#'
#' @examples
#'pathways<-nodeData_keggPathwaygraphs(node="hsa:208",updateKEGG=FALSE)
#' @import graph
#' @import ROntoTools
#' @export

nodeData_keggPathwaygraphs<-function(node,updateKEGG=FALSE)
    {

    kpg<-keggPathwayGraphs(organism ="hsa",updateCache = updateKEGG)
    kpgN<-keggPathwayNames(organism ="hsa" ,updateCache =updateKEGG )
    kpg <- setEdgeWeights(kpg, edgeTypeAttr = "subtype",
    edgeWeightByType = list(activation = 1, inhibition = -1,
    expression = 1, repression = -1),
    defaultWeight = 0)

    node_exist<-sapply(kpg,function(x)
    {
        if(node %in% nodes(x))
            return(TRUE)
        else
        return(FALSE)
    })

    res<-kpgN[names(node_exist)[node_exist]]

    if(!length(res))
        stop(paste("The Node '",node, "' does not belong to
            KEGG pathways",sep=""))
    else
    return(res)
    }
