#' utils: utility functions
#' names in which there is an edge between the given nodes
#'
#' @param updateKEGG re-download KEGG data
#' @param edgeFrom A character of node name
#' @param edgeTo A character of node name
#'
#' @details
#'
#' #' Draghici S., et.al.: "A systems biology approach for
#' pathway level analysis". Genome Research, 17, 2007.
#'
#' @return
#'
#' A character vector of KEGG pathway names in which
#' edge (edgeFrom, edgeTo) exists
#'
#' @author
#'
#' Azam Peyvandipour and Sorin Draghici
#'
#'
#' @seealso \code{\link{keggPathwayGraphs},\link{keggPathwayNames}}
#'
#' @examples
#'pathways<-edgeData_keggPathwaygraphs(edgeFrom="hsa:208",
#'edgeTo="hsa:2309",updateKEGG=FALSE)
#'class(pathways)
#'head(pathways)
#' @import graph
#' @import ROntoTools
#' @export
edgeData_keggPathwaygraphs<-function(edgeFrom,edgeTo,updateKEGG=FALSE)
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
        if(all(c(edgeFrom,edgeTo) %in% rownames(x)))
            {
            if(x[edgeFrom,edgeTo]!=0)
                return(TRUE)
            }
        else
        return(FALSE)
        }))

    res<-kpgN[names(edgeExist)[which(edgeExist)]]
    if(!length(res))
        stop(paste("The edge (",edgeFrom," , " ,edgeTo,")
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
#'class(pathways)
#'head(pathways)
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
