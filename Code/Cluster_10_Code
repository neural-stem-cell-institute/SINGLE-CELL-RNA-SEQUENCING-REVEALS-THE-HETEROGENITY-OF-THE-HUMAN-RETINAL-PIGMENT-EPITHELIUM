library(hypeR)
library(doParallel)
library(org.Hs.eg.db)
library(pCalibrate)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(riverplot)
library(foreach)
library(GO.db)
library(LSAfun)
library(rrvgo)
library(GOfuncR)
library(reshape2)
library(dplyr)
library(igraph)
library(mgcv)

###Analysis of cluster 10 sing the ICELL8 data
###first map integrated clusters onto ICELL8 object
obj<-ICELL8.nomac
x<-int.obj@active.ident[grep("cell_",names(int.obj@active.ident))]
obj@active.ident<-x

#########Identify genes associated with each cluster

I8.markers.int<-FindAllMarkers(obj,only.pos=T,
                           likeFC.threshold=0.5)

x<-I8.markers.int[I8.markers.int$p_val_adj<0.01,]
I8.clust.int<-foreach(i=0:12) %do% {x[x$cluster==i,]$gene}

#####run slingshot ananlysis
####slingshot with ICELL8 integrated clustering
obj<-RunUMAP(obj,dims=1:30)
z<-obj@reductions$umap@cell.embeddings
y<-(as.vector(obj@active.ident))
names(y)<-names(obj@active.ident)

fbi<-slingshot(z,y)

#####run enrichment analysis
cluster.enr.int.i8<-foreach(i=1:length(I8.clust)) %do% {
  Path_BIOCARTA<-hypeR(I8.clust[[i]], BIOCARTA, background=23459, fdr=0.01)
  Path_Reactome<-hypeR(I8.clust[[i]], REACTOME, background=23459, fdr=0.01)
  Path_KEGG<-hypeR(I8.clust[[i]], KEGG, background=23459, fdr=0.01)
  GOBioProc<-hypeR(I8.clust[[i]], GOBP, background=23459, fdr=0.01)
  GOMoleFunc<-hypeR(I8.clust[[i]], GOMF, background=23459, fdr=0.01)
  GOCellComp<-hypeR(I8.clust[[i]], GOCC, background=23459, fdr=0.01)
  list(Path_BIOCARTA,Path_Reactome,Path_KEGG,GOBioProc,GOMoleFunc,GOCellComp)
}
names(cluster.enr.int.i8)<-paste("Cluster",c(0:11),sep="_")
cluster.enr.i8.data<-foreach(i=1:12) %do% {
  Path1<-cluster.enr.int.i8[[i]][[1]]$data
  Path2<-cluster.enr.int.i8[[i]][[2]]$data
  Path3<-cluster.enr.int.i8[[i]][[3]]$data
  Pathways<-rbind(Path1,Path2,Path3)
  BP<-cluster.enr.int.i8[[i]][[4]]$data
  MF<-cluster.enr.int.i8[[i]][[5]]$data
  CC<-cluster.enr.int.i8[[i]][[6]]$data
  list(Pathways,BP,MF,CC)
}
names(cluster.enr.i8.data)<-paste("Cluster",c(0:11),sep="_")


#######select cluster 10 genes for network analysis
x<-I8.markers.int[I8.markers.int$p_val_adj<0.01,]
x<-x[x$cluster==10,]
x<-x[x$avg_logFC>1,]
write.csv(x,"I8_cl10.csv")
###get interactions from StringDb and load into R
cl10.net<-data.frame(read.delim("string_cl10.tsv",as.is=T))


#####prep node attributes
x<-unique(c(cl10.net$X.node1,cl10.net$node2))
y<-intersect(x,TFS)
y1<-foreach(i=1:length(x),.combine='c') %do%{
  if(length(intersect(x[i],y))>0) {"TF"} else {"Non"}
}
y<-data.frame(x,y1)
y<-y%>% distinct()
colnames(y)<-c("node","Type")

gph10<-graph_from_data_frame(cl10.net,directed=F,vertices=y)

######chose louvain clustering algorithm based on PMID:31029085
x<-cluster_louvain(gph10,weights=E(gph10)$combined_score)
foreach(i=1:length(x),.combine='c')%do% {length(x[[i]])}
###[1] 138 102  24  25   2 189   2
##drop communities smaller than 20 and save genes for later removal from edgelist
xr<-c(x[[5]],x[[7]])
x<-x[-7]
x<-x[-5]

#####enrichments per community in network
community_cl10<-foreach(i=1:length(x)) %do% {
  Path_BIOCARTA<-hypeR(x[[i]], BIOCARTA, background=23459, fdr=0.01)
  Path_Reactome<-hypeR(x[[i]], REACTOME, background=23459, fdr=0.01)
  Path_KEGG<-hypeR(x[[i]], KEGG, background=23459, fdr=0.01)
  GOBioProc<-hypeR(x[[i]], GOBP, background=23459, fdr=0.01)
  GOMoleFunc<-hypeR(x[[i]], GOMF, background=23459, fdr=0.01)
  GOCellComp<-hypeR(x[[i]], GOCC, background=23459, fdr=0.01)
  list(Path_BIOCARTA,Path_Reactome,Path_KEGG,GOBioProc,GOMoleFunc,GOCellComp)
}
names(community_cl10)<-paste("Community",c(1:length(x)),sep="_")
community_cl10.data<-foreach(i=1:length(x)) %do% {
  Path1<-community_cl10[[i]][[1]]$data
  Path2<-community_cl10[[i]][[2]]$data
  Path3<-community_cl10[[i]][[3]]$data
  Pathways<-rbind(Path1,Path2,Path3)
  BP<-community_cl10[[i]][[4]]$data
  MF<-community_cl10[[i]][[5]]$data
  CC<-community_cl10[[i]][[6]]$data
  list(Pathways,BP,MF,CC)
}
names(community_cl10.data)<-paste("Community",c(1:length(x)),sep="_")
#######make a list of transcription factors based on GO categroies
nuc<-get_anno_genes("GO:0005634")
TFS<-get_anno_genes("GO:0006366")
RIBO<-get_anno_genes("GO:0005840")
PROT<-get_anno_genes("GO:0000502")
TFS<-setdiff(TFS[,2],RIBO[,2])
TFS<-setdiff(TFS,PROT[,2])
TFS<-intersect(TFS,nuc[,2])

############make a new graph object for export to cytoscape
z<-foreach(i=1:length(x)) %do% {
  print(i)
  y<-rownames(community_cl10.data[[i]][[2]])
  if(length(y)==0) {} else{
  y<-unlist(substring(y,4))
  y<-breakdown(y)
  y<-goterms[y]
  scores<-community_cl10.data[[i]][[2]]$overlap
  scores<-scores[!is.na(y)]
  sizes<-community_cl10.data[[i]][[2]]$geneset
  sizes<-sizes[!is.na(y)]
  y<-y[!is.na(y)]
  names(sizes)<-y
  names(scores)<-y
  try(sm<-calculateSimMatrix(y,orgdb="org.Hs.eg.db",ont="BP",method="Rel"))
  sizes<-sizes[intersect(rownames(sm),names(sizes))]
  scores<-scores[intersect(rownames(sm),names(scores))]
  o <- rev(order(scores, sizes, na.last = FALSE))
  sm <- sm[o, o]
  cluster <- cutree(hclust(as.dist(1 - sm)), h = .7)
  clusterRep <- tapply(rownames(sm), cluster, function(x) x[which.max(scores[x])])
  red<-data.frame(go = rownames(sm), cluster = cluster, parent = clusterRep[cluster],
                  parentSimScore = unlist(Map(seq_len(nrow(sm)), 
                                              clusterRep[cluster], f = function(i, j) sm[i,j])), 
                  score = scores[match(rownames(sm), names(scores))], 
                  size = sizes[match(rownames(sm), names(sizes))], 
                  term = rrvgo:::getGoTerm(rownames(sm)), 
                  parentTerm = rrvgo:::getGoTerm(clusterRep[cluster]))
  }
  
}  

z1<-foreach(i=1:length(x),.combine='rbind') %do% {
  y<-names(sort(table(z[[i]]$parentTerm),decreasing=T))
  if(is.null(y)) {
    x2<-x[[i]]
    pTerm<-rep("No_Cat",length(x2))
    fin<-data.frame(x2,pTerm)
    fin$cm<-i
    fin
  }else {
  x1<-x[[i]]
  b<-1
  y1<-get_anno_genes(z[[i]][z[[i]]$parentTerm==y[b],]$go)$gene
  x2<-intersect(x1,y1)
  x1<-setdiff(x1,x2)
  pTerm<-rep(y[b],length(x2))
  fin<-data.frame(x2,pTerm)
  repeat{
    b<-b+1
    print(b)
    y1<-get_anno_genes(z[[i]][z[[i]]$parentTerm==y[b],]$go)$gene
    x2<-intersect(x1,y1)
    x1<-setdiff(x1,x2)
    pTerm<-rep(y[b],length(x2))
    m<-data.frame(x2,pTerm)
    fin<-rbind(fin,m)
    if(length(x1)==0||is.na(y[b+1])) {break}
  }
  x2<-x1
  pTerm<-rep("No_Cat",length(x2))
  m<-data.frame(x2,pTerm)
  fin<-rbind(fin,m)
  fin$cm<-i
  fin
  }
}

y<-intersect(z1$x2,TFS)
y1<-foreach(i=1:length(z1$x2),.combine='c') %do%{
  if(length(intersect(z1$x2[i],y))>0) {"TF"} else {"Non"}
}

z1$Type<-y1
colnames(z1)<-c("node","pTerm","cm","Type")
####remove rows from cl.net corresponding to the removed communities

x<-cl10.net
for(i in 1:length(xr)) {
  x<-x[-(which(x$X.node1==xr[i])),]
  x<-x[-(which(x$node2==xr[i])),]
}



gph10<-graph_from_data_frame(x,directed=F,vertices=z1)

write_graph(gph10,"GPH10.graphml",format="graphml")
####use cytoscape to generate figure for network

######make Vlnplots for figure
StackedVlnPlot(obj,features=c("GATA2","GATA1","MYC","BMI1","REST"))
StackedVlnPlot(obj,features=c("FOXO3","MEF2C","STAT3","CHD7","MYB"))








