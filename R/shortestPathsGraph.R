#' shortestPathsGraph: Construct a subgraph of KEGG global
#' graph connecting drug targets and disease-related genes.
#'
#' @param ggraph KEGG global graph as object of type \code{\link{graphNEL}}
#' @param drug_targets the vector of drug target genes
#' @param disease_genes the vector of disease-related genes
#'
#' @references
#'
#' Draghici S., et.al.: "A systems biology approach for
#' pathway level analysis".
#' Genome Research, 17, 2007.
#'
#' @return
#'
#' An object of class \code{\link{graphNEL}} encoding
#' the a subgraph of KEGG global graph
#'
#' @author
#'
#' Azam Peyvandipour and Sorin Draghici
#'
#'
#' @seealso \code{\link{keggGlobalGraph}}
#'
#' @examples
#' #obtain KEGG global graph gg
#' gg<-keggGlobalGraph()
#' #drug target genes can be obtained from The Comparative
#' #Toxicogenomics Database (CTD)
#' drug_targets<-c("hsa:9368","hsa:2322", "hsa:3932", "hsa:4067", "hsa:6714")
#'
#' #disease-related genes can be obtained from The Comparative
#' #Toxicogenomics Database (CTD)
#' disease_genes<-c("hsa:1832" ,"hsa:5073" ,"hsa:5328")
#'
#' #obtain drug-disease network
#' mygraph<-shortestPathsGraph(drug_targets,disease_genes,gg)
#' @import graph
#' @import igraph
#' @import stringr
#' @export
shortestPathsGraph<-function(drug_targets,disease_genes,ggraph){
    subGlobalGraph<-NULL
    #negative weights are identified and converted to positive temporarily

    if(any(is.na(str_match(drug_targets,pattern="hsa:[0-9]+"))))
        stop("drug_targets must start with 'hsa' and end with numbers.")
    if(any(is.na(str_match(disease_genes,pattern="hsa:[0-9]+"))))
        stop("disease_genes must start with 'hsa' and end with numbers.")

    adjmatrix<-as(ggraph,"matrix")
    globalGraphNonNeg<-ggraph
    mat<-adjmatrixNeg<-which(adjmatrix == -1,arr.ind=TRUE)
    mat[,1]<-rownames(adjmatrix)[adjmatrixNeg[,1]]
    mat[,2]<-rownames(adjmatrix)[adjmatrixNeg[,2]]

    edgeData(globalGraphNonNeg,from=as.character(mat[,1]),
    to=as.character(mat[,2]),attr="weight")<-2
    globaliGraphNonNeg<-igraph.from.graphNEL(globalGraphNonNeg)

    drug_targets<-drug_targets[which(drug_targets %in% nodes(ggraph))]
    disease_genes<-disease_genes[which(disease_genes %in% nodes(ggraph))]
    options(warn=-1)
    if(length(drug_targets)>0 & length(disease_genes)>0){
        allSelectedNodes<-unique(c(drug_targets,drug_targets))
        allSelectedNodes<-allSelectedNodes[which(
        allSelectedNodes %in% nodes(globalGraphNonNeg))]
        globaliGraph<-igraph.from.graphNEL(globalGraphNonNeg)
        alshortestpathdis2drug<-list()
        alshortestpathdrug2dis<-list()
        for(j in 1:length(disease_genes)){
            ap<-get.all.shortest.paths(globaliGraphNonNeg,from=disease_genes[j],
            to=drug_targets,mode="out")
            ap$res<-lapply(ap$res,function(x)
            {return(nodes(globalGraphNonNeg)[x])})
            alshortestpathdis2drug<-c(ap$res,alshortestpathdis2drug)
        }

    for(j in 1:length(drug_targets)){
        ap<-get.all.shortest.paths(globaliGraph,
        from=drug_targets[j],to=disease_genes,mode="out")
        ap$res<-lapply(ap$res,function(x){return(nodes(globalGraphNonNeg)[x])})
        if(length(ap$res)>0)
        alshortestpathdrug2dis<-c(ap$res,alshortestpathdrug2dis)
    }
    alshortestpath<-unique(c(alshortestpathdis2drug,alshortestpathdrug2dis))
    if(length(alshortestpath)>0){
        maxl<-max(sapply(alshortestpath,length))
        if(maxl>0)
        {
            alshortestpath2<-lapply(alshortestpath,function(x){
            z<-c()
            z<-x
            z[(length(x)+1):maxl]<-"NULL"
            return(z)} )

    alshortestpathMat <- do.call(rbind, alshortestpath2)
    alshortestpathMat<-alshortestpathMat[,-dim(alshortestpathMat)[2]]
    options(warn=0)
    adjmatrixSubGlobalGraph<-c()
    for( k in 1:(dim(alshortestpathMat)[2]-1))
        {
        z<-alshortestpathMat[,k:(k+1)]
        adjmatrixSubGlobalGraph<-rbind(adjmatrixSubGlobalGraph,z)
        }

    adjmatrixSubGlobalGraph<-
    adjmatrixSubGlobalGraph[which(adjmatrixSubGlobalGraph[,1] != "NULL"),]
    adjmatrixSubGlobalGraph<-
    adjmatrixSubGlobalGraph[which(adjmatrixSubGlobalGraph[,2] != "NULL"),]
    adjmatrixSubGlobalGraph<-
    adjmatrixSubGlobalGraph[!duplicated(adjmatrixSubGlobalGraph),]
    allNodes<-unique(as.character(adjmatrixSubGlobalGraph))
    subGlobalGraph <- ftM2graphNEL(adjmatrixSubGlobalGraph,
    V=allNodes,edgemode="directed")

    edgeDataDefaults(subGlobalGraph, "subtype") <- NA
    edgeDataDefaults(subGlobalGraph, "weight") <- NA

    #We revert the edge weights
    edgeData(subGlobalGraph,from=as.character(adjmatrixSubGlobalGraph[,1]),
    to=as.character(adjmatrixSubGlobalGraph[,2]),attr="weight")<-
    edgeData(ggraph,from=as.character(adjmatrixSubGlobalGraph[,1]),
    to=as.character(adjmatrixSubGlobalGraph[,2]),attr="weight")
    edgeData(subGlobalGraph,from=as.character(adjmatrixSubGlobalGraph[,1]),
    to=as.character(adjmatrixSubGlobalGraph[,2]),attr="subtype")<-
    edgeData(ggraph,from=as.character(adjmatrixSubGlobalGraph[,1]),
    to=as.character(adjmatrixSubGlobalGraph[,2]),attr="subtype")
    }    else stop("Drug target genes and disease-related genes
        should be more than one.")

        }

    }
    else
    {
    if(length(drug_targets) == 0)
    stop("drug_targets must belong to KEGG.")

    if(length(disease_genes) == 0)
    stop("disease_genes must belong to KEGG.")
    }

    return(subGlobalGraph) }

