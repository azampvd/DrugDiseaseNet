### R code from vignette source 'DrugDiseaseNet.Rnw'

###################################################
### code chunk number 1: DrugDiseaseNet.Rnw:92-95
###################################################
require(ROntoTools)
require(DrugDiseaseNet)
gg<-keggGlobalGraph()


###################################################
### code chunk number 2: DrugDiseaseNet.Rnw:102-105
###################################################
gg<-keggGlobalGraph(updateKEGG=TRUE)
gg


###################################################
### code chunk number 3: DrugDiseaseNet.Rnw:113-117
###################################################
require(graph)
allEdges<-edges(gg)
allEdges$"hsa:2065"
head(nodes(gg))


###################################################
### code chunk number 4: DrugDiseaseNet.Rnw:125-126
###################################################
edgeData(gg,from="hsa:2065",to = "hsa:2549" )


###################################################
### code chunk number 5: DrugDiseaseNet.Rnw:132-135
###################################################
adjmatrix<-as(gg,"matrix")
adjmatrix [1:4,1:4]
adjmatrix ["hsa:2065","hsa:2549"]

###################################################
### code chunk number 6: DrugDiseaseNet.Rnw:147-154
###################################################
gg<-keggGlobalGraph()
#drug target genes can be obtained from The CTD
drug_targets<-c("hsa:9368","hsa:2322", "hsa:3932",
"hsa:4067", "hsa:6714")
#disease-related genes can be obtained from CTD
disease_genes<-c("hsa:1832" ,"hsa:5073" ,"hsa:5328")
drugDiseaseNetwork<-shortestPathsGraph(drug_targets,
disease_genes,gg)
drugDiseaseNetwork
