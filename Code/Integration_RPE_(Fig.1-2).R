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
library(slingshot)
library(igraph)
library(gam)
library(dplyr)



load("10x.RData")
load("ICELL8-nomac.RData")
options(future.globals.maxSize=2000 * 1024 ^ 2)

#########functions
source("stackedViolinPlot.R")
#########ggplot color replicator
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#####make a riverplot figure

quickriver<-function(obj.int,obj1,obj2) {
  r1<-foreach(i=0:max(as.numeric(levels(obj1@active.ident))),.combine='rbind') %do% {
    x<-table(obj.int@active.ident[names(obj1@active.ident[obj1@active.ident==i])])
    x<-x/sum(x)
    cbind(paste("A",rep(i,length(x)),sep=""),paste("I",
                                                   as.numeric(levels(obj.int@active.ident)),sep=""),x)
  }
  
  
  r2<-foreach(i=0:max(as.numeric(levels(obj2@active.ident))),.combine='rbind') %do% {
    x<-table(obj.int@active.ident[names(obj2@active.ident[obj2@active.ident==i])])
    x<-x/sum(x)
    cbind(paste("I",as.numeric(levels(obj.int@active.ident)),sep=""),
          paste("B",rep(i,length(x)),sep=""),x)
  }
  edges<-rbind(r1,r2)
  edges<-data.frame(edges[,1],edges[,2],as.numeric(edges[,3]))
  colnames(edges)<-c("N1","N2","Value")
  nodes<-data.frame(ID=unique(c(as.character(edges$N1),as.character(edges$N2))), stringsAsFactors = F)
  nodes$x<-c(rep(0,max(as.numeric(levels(obj1@active.ident)))+1),
             rep(5,max(as.numeric(levels(obj.int@active.ident)))+1),
             rep(10,max(as.numeric(levels(obj2@active.ident)))+1))
  
  
  pal1<-paste0(colorRampPalette(c("royalblue1","lawngreen"))(max(as.numeric(levels(obj1@active.ident)))+1),"80")
  pal2<-paste0(colorRampPalette(c("mediumorchid1","goldenrod3"))(max(as.numeric(levels(obj.int@active.ident)))+1),"80")
  pal3<-paste0(colorRampPalette(c("cyan","red4"))(max(as.numeric(levels(obj2@active.ident)))+1),"80")
  palette<-c(pal1, pal2,pal3)
  
  nodes$col<-palette
  styles<-foreach(i=1:length(palette)) %do% {
    list(col=palette[i],textcol="black",srt=0,nodestyle="regular",lty=0,nodewidth=1,
         yscale=0.5)}
  names(styles) = nodes$ID
  
  edges<-edges[edges$Value>0,]
  rp<-makeRiver(nodes,edges,node_style=styles)
  
  riverplot(rp,plot_area=1,add_mid_points=T, gravity="center")
}
####annotation and set up

registerDoParallel(makeCluster(3))
SYMBOLS<-unlist(as.list(org.Hs.egSYMBOL))

x <- Term(GOTERM)
goterms<-names(x)
names(goterms)<-x

BIOCARTA <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:BIOCARTA")
KEGG     <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:KEGG")
REACTOME <- msigdb_gsets(species="Homo sapiens", category="C2", subcategory="CP:REACTOME")
paths<-c(BIOCARTA, KEGG, REACTOME)
GOBP <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory="BP")
GOMF <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory="MF")
GOCC <- msigdb_gsets(species="Homo sapiens", category="C5", subcategory="CC")

#######integrating data from 10x and ICELL8 using Seurat Package 
obj.list<-list()

obj.list[[1]]<-SCTransform(ICELL8.nomac,verbose=F)
obj.list[[2]]<-SCTransform(RNA.10x,verbose=F)

m1<-SelectIntegrationFeatures(object.list = obj.list, nfeatures = 4000)
m<-PrepSCTIntegration(object.list = obj.list, anchor.features = m1 , 
                      verbose = FALSE)
m4<-FindIntegrationAnchors(object.list = m, normalization.method = "SCT", 
                           anchor.features = m1, verbose = FALSE,reduction="cca")
int.obj<-IntegrateData(anchorset = m4, normalization.method = "SCT", 
                       verbose = FALSE)
int.obj<-RunPCA(int.obj,verbose=F)
int.obj<-RunUMAP(int.obj,dims=1:30)

int.obj<- FindNeighbors(int.obj, dims = 1:10)
int.obj <- FindClusters(int.obj, resolution = 0.5)
DefaultAssay(int.obj)<-"RNA"
int.obj<-NormalizeData(int.obj,verbose=F)
int.obj<-ScaleData(int.obj,features=rownames(int.obj@assays$RNA@counts))
DimPlot(int.obj,group.by=c("tech","ident"),label=T)
DimPlot(int.obj,group.by="seurat_clusters",label=T)


quickriver(int.obj,RNA.10x,ICELL8.nomac)


#########cluster markers
ser.m.markers<-FindAllMarkers(int.obj,only.pos=T,likeFC.threshold=0.05)
x<-ser.m.markers[ser.m.markers$p_val_adj<0.1,]
ser.clust<-foreach(i=0:12) %do% {x[x$cluster==i,]$gene}
write.csv(x,"Ser.int.markers.csv")

########RPE signiture analysis

x<-read.csv("human_RPE_signature_genes_26517551.csv",as.is=T)
y<-intersect(rownames(int.obj@assays$SCT@scale.data),x[,1])
test<-probs(int.obj)
m<-test[[1]][which((apply(test[[1]][intersect(rownames(test[[1]]),x[,1]),],1,max)>0)),]

z<-heatmap(int.obj@assays$SCT@scale.data[y,])
DoHeatmap(int.obj,features=rev(y[z$rowInd]))+scale_fill_distiller(palette="RdYlBu")

foreach(i=0:max(as.numeric(levels(int.obj@active.ident))),.combine='c') %do% {
  length(which(int.obj@active.ident==i))}

m<-foreach(i=1:dim(m)[2]) %do% {names(which(m[,i]>0))}

names(m)<-paste("Cluster",0:12,sep=" ")
y<-fromList(m)

upset(y,nsets=13,nintersects = NA,
      matrix.dot.alpha=0.5,order.by="freq")

######enrichment analysis
cluster.enr.int.ser<-foreach(i=1:length(ser.clust)) %do% {
  Path_BIOCARTA<-hypeR(ser.clust[[i]], BIOCARTA, background=23459, fdr=0.01)
  Path_Reactome<-hypeR(ser.clust[[i]], REACTOME, background=23459, fdr=0.01)
  Path_KEGG<-hypeR(ser.clust[[i]], KEGG, background=23459, fdr=0.01)
  GOBioProc<-hypeR(ser.clust[[i]], GOBP, background=23459, fdr=0.01)
  GOMoleFunc<-hypeR(ser.clust[[i]], GOMF, background=23459, fdr=0.01)
  GOCellComp<-hypeR(ser.clust[[i]], GOCC, background=23459, fdr=0.01)
  list(Path_BIOCARTA,Path_Reactome,Path_KEGG,GOBioProc,GOMoleFunc,GOCellComp)
}
names(cluster.enr.int.ser)<-paste("Cluster",c(0:12),sep="_")
cluster.enr.ser.data<-foreach(i=1:13) %do% {
  Path1<-cluster.enr.int.ser[[i]][[1]]$data
  Path2<-cluster.enr.int.ser[[i]][[2]]$data
  Path3<-cluster.enr.int.ser[[i]][[3]]$data
  Pathways<-rbind(Path1,Path2,Path3)
  BP<-cluster.enr.int.ser[[i]][[4]]$data
  MF<-cluster.enr.int.ser[[i]][[5]]$data
  CC<-cluster.enr.int.ser[[i]][[6]]$data
  list(Pathways,BP,MF,CC)
}
names(cluster.enr.ser.data)<-paste("Cluster",c(0:12),sep="_")
enr<-c("Pathways","GO:BP","GO:MF","GO:CC")
x1<-foreach(i=1:length(cluster.enr.ser.data), .combine='rbind') %do% {
  y<-cluster.enr.ser.data[[i]]
  
  z<-foreach(b=1:4, .combine='rbind') %do% {
    m<-y[[b]]
    if(dim(m)[1]==0) {} else {m$Type<-enr[b]}
    m
  }
  
  if(dim(z)[1]==0) {} else {z$Cluster<-names(cluster.enr.int.ser)[i]}
  z
}
write.csv(x1,"Supplemental_Table_int_enr_ser.csv",row.names=F)

######GO visualization of parent terms after semantic simialirty analysis
Go.list.viz<-foreach(i=1:13) %do% {
  y<-rownames(cluster.enr.ser.data[[i]][[2]])
  if(length(y)==0) {"No Enrichments"} else{
    y<-unlist(substring(y,4))
    y<-breakdown(y)
    y<-goterms[y]
    scores<-cluster.enr.ser.data[[i]][[2]]$overlap
    scores<-scores[!is.na(y)]
    sizes<-cluster.enr.ser.data[[i]][[2]]$geneset
    sizes<-sizes[!is.na(y)]
    y<-y[!is.na(y)]
    names(sizes)<-y
    names(scores)<-y
    sm<-calculateSimMatrix(y,orgdb="org.Hs.eg.db",ont="BP",method="Rel")
    sizes<-sizes[intersect(rownames(sm),names(sizes))]
    scores<-scores[intersect(rownames(sm),names(scores))]
    o <- rev(order(scores, sizes, na.last = FALSE))
    sm <- sm[o, o]
    cluster <- cutree(hclust(as.dist(1 - sm)), h = .7)
    clusterRep <- tapply(rownames(sm), cluster, function(x) x[which.max(scores[x])])
    red<-data.frame(go = rownames(sm), cluster = cluster, parent = clusterRep[cluster],parentSimScore = unlist(Map(seq_len(nrow(sm)), 
                                                                                                                   clusterRep[cluster], f = function(i, j) sm[i,j])), score = scores[match(rownames(sm), 
                                                                                                                                                                                           names(scores))], size = sizes[match(rownames(sm), names(sizes))], term = rrvgo:::getGoTerm(rownames(sm)), 
                    parentTerm = rrvgo:::getGoTerm(clusterRep[cluster]))
  }
}  
x<-Go.list.viz[-4]  

z<-unique(x[[1]]$parentTerm)
for( i in 2:12) {
  z<-c(z,unique(x[[i]]$parentTerm))
}
z<-sort(table(z))

y<-c(names(sort(z)[143:149]),"cell cycle","detection of stimulus","ion transport","growth","sensory perception of light stimulus")

x<-foreach(i=1:length(z)) %do% {
  fin<-foreach(b=1:12,.combine='c') %do% {
    g1<-x[[b]]
    g1<-g1[which(g1$parentTerm==names(z)[i]),]$go
    
  }
  unique(fin)
}


cats<-GOBP$genesets
y<-substring(names(cats),4)
y<-tolower(y)
y<-gsub("_"," ",y)
y<-goterms[y]
names(cats)<-y

cats2<-foreach(i=1:length(x)) %do% {
  x1<-x[[i]]
  x2<-foreach(b=1:length(x1),.combine='c') %do% {
    cats[[x1[b]]]
  }
  unique(x2)
}



names(cats2)<-names(z)


G1<-foreach(i=1:13) %do% {
  g3<-hypeR(ser.clust[[i]], cats2, background=23459, fdr=1)
  Percentage<-g3$data$overlap/g3$data$geneset*100
  FDR<-g3$data$fdr
  data.frame(g3$data$label,Percentage,FDR)
}

G2<-foreach(i=1:13,.combine='rbind') %do% {
  GO_Names<-G1[[i]]$g3.data.label
  Values<-G1[[i]]$Percentage
  Clusters<-as.character(rep(i-1,length(Values)))
  FDR<-G1[[i]]$FDR
  data.frame(GO_Names,Values,Clusters,FDR)
}
G2$Clusters<-factor(G2$Clusters,
                    levels=c("0","1","2","3","4","5","6","7","8","9","10","11","12"))

y<-c(rev(names(sort(z)[143:149])),"ion transport","cell cycle","detection of stimulus","growth","sensory perception of light stimulus")
y<-c(y,"regulation of stem cell differentiation")

G3<-foreach(i=1:length(y),.combine='rbind') %do% {
  G2[G2$GO_Names==y[i],]
}
G3$GO_Names<-factor(G3$GO_Names,levels=y)


p<-ggplot(G3, aes(GO_Names,-log10(FDR)))
p<-p+theme_light()
p<-p+geom_jitter(height=0,width=0.3,aes(color=Clusters,size=Values))
p<-p +theme(axis.text.x = element_text(angle=90,size=10),axis.title.x = element_blank()) 
p<-p+ylab("-log10(FDR)")
p<-p+geom_hline(yintercept=-log10(0.1),color="red",linetype="dashed")

p

###Reactome pathway analysis and visualization
RP<-read.delim("ReactomePathways.txt",as.is=T,header=F)
RP<-RP[grep("Homo",RP[,3]),]
rownames(RP)<-RP$V1

###data prep of hypeR enrichment data
clx<-foreach(i=1:13,.combine='rbind') %do% {
  x<-cluster.enr.ser.data[[i]][[1]]
  x$Cluster<-i-1
  x[grep("REACTOME",x$label),]
}
x<-strsplit(clx$label,"REACTOME_")
x<-unlist(x)[seq(2,length(unlist(x)),2)]
x<-gsub("_", " ", x)
clx$label1<-x
RP$caps<-toupper(RP$V2)
RP$caps<-gsub("\\(","",RP$caps)
RP$caps<-gsub("\\)","",RP$caps)
RP$caps<-gsub("-"," ",RP$caps)
RP$caps<-gsub(",","",RP$caps)
RP$caps<-gsub("\\+","",RP$caps)
RP$caps<-gsub(":"," ",RP$caps)
RP$caps<-gsub("\\/"," ",RP$caps)
RP$caps<-trimws(RP$caps)
clx$PathID<-foreach(i=1:length(clx$label1),.combine="c") %dopar% {
  x<-RP[which(RP$caps==clx$label1[i]),1]
  if(length(x)==0) {x<-NA} else{x}
}

x<-c("R-HSA-983170","R-HSA-9006936","R-HSA-380994","R-HSA-450282",
     "R-HSA-163200","R-HSA-163200","R-HSA-163200","R-HSA-163200","R-HSA-163200",
     "R-HSA-187577","R-HSA-389958","R-HSA-72165","R-HSA-5632684","R-HSA-5632684",
     "R-HSA-983169","R-HSA-380994","R-HSA-75035","R-HSA-429958",
     "R-HSA-163200","R-HSA-163200","R-HSA-5576892","R-HSA-75035")

clx$PathID[which(is.na(clx$PathID))]<-x

##prep of Reactome Pathway relation data to either root or 2nd level pathways in the tree
RPH<-read.delim("ReactomePathwaysRelation.txt",as.is=T,header=F)
RPH<-RPH[grep("HSA",RPH[,1]),]
RPH.adj<-RPH[-(grep("R-HSA-8953897",RPH[,1])),]
RPH.adj<-RPH.adj[-(grep("R-HSA-168256",RPH.adj[,1])),]
RPH.adj<-RPH.adj[-(grep("R-HSA-74160",RPH.adj[,1])),]
RPH.adj<-RPH.adj[-(grep("R-HSA-1266738",RPH.adj[,1])),]
RPH.adj<-RPH.adj[-(grep("R-HSA-5357801",RPH.adj[,1])),]

###generating a data frame of enriched pathways with  high level terms 
x<-foreach(i=1:length(clx$PathID),.combine='c') %do% {
  y<-RPH.adj[which(RPH.adj[,2]==clx$PathID[i]),1]
  if(length(y)==0) {y<-clx$PathID[i]} else {y[1]}
}

repeat{
  z1<-length(unique(x))
  x<-foreach(i=1:length(x),.combine='c') %do% {
    y<-RPH.adj[which(RPH.adj[,2]==x[i]),1]
    if(length(y)==0) {y<-x[i]} else {y[1]}
  }
  if(length(unique(x))==z1) {break}
}

clx$rootPath<-x
clx$rootName<-RP[x,]$V2

#####prepping a data frame for visualization
x1<-table(x)
x<-foreach(i=1:length(names(x1))) %do% {
  y<-RPH[which(RPH[,1]==names(x1[i])),2]
  if(length(y)==0) {y<-names(x1)[i]} else {y}
}
x<-foreach(b=1:length(x)) %do% {
  x2<-x[[b]]
  x3<-x2
  repeat{
    
    x3<-foreach(i=1:length(x3),.combine='c') %do% {
      y<-RPH.adj[which(RPH.adj[,1]==x3[i]),2]
    }
    if(length(x3)==0) {break} else {x2<-c(x2,x3)}
    
  }
  x2
}

x<-foreach(i=1:length(x),.combine='c') %do% {length(x[[i]])}
names(x)<-RP[names(x1),]$V2
x1<-x




x<-unique(clx$rootName)
x<-c(x[4],x[9],x[15],x[20],
     x[6],x[17],x[8],x[13],x[26],
     x[7],x[1],x[10],
     x[22:24],x[25],x[29],
     x[18],x[2],x[30],x[11:12],x[16],x[3],
     x[19],x[14],x[21],x[28],x[27],x[5])
z<-foreach(i=0:12,.combine='cbind') %do% {
  y<-table(clx[clx$Cluster==i,]$rootName)
  z<-setdiff(x,names(y))
  z1<-rep(0,length(z))
  names(z1)<-z
  y<-c(y,z1)
  y<-y[x]
  
}
colnames(z)<-0:12
z1<-apply(z,2,sum)
z2<-apply(z,2,function(x) x/sum(x)*100)


y<-melt(z)
colnames(y)<-c("Path","Cluster","P.N")
y1<-melt(z2)
y$percent<-y1$value
y1<-z/x1*100
y1<-melt(y1)
y$Tree.per<-y1$value
y<-y[-(which(y$`P.N`==0)),]
y1<-y$Cluster
for( i in 0:12) {y1[which(y1==i)]<-rev(gg_color_hue(13))[i+1]}
y$colors<-y1  

###some clean up, removing top level terms where 2nd level terms used, and 
###fixing the "GPCR downstream signaling" to be included in parent term
y<-y[-(which(y$Path=="Developmental Biology")),]
y<-y[-(which(y$Path=="Cellular responses to external stimuli")),]
y<-y[-(which(y$Path=="Programmed Cell Death")),]
y["379",]$P.N<-8
y["379",]$percent<-8/71*100
y["379",]$Tree.per<- 8/399*100
y<-y[-115,]

p<-ggplot(y, aes(y=factor(Cluster), x=Path, size=percent, alpha=Tree.per,color=colors))
p<-p + geom_point() + geom_text(aes(label=P.N,vjust=-0.5), alpha=1.0, size=4) 
p<-p + scale_alpha_continuous(range=c(0.5, 1)) + scale_size_area(max_size = 6)
p<-p + theme_bw() + theme(axis.line = element_blank(),        
                          axis.title = element_blank(),         
                          panel.border = element_blank(),         
                          panel.grid.major.x = element_blank(),  
                          panel.grid.minor.x = element_blank(),
                          axis.text.x= element_text(angle=90,vjust=1,hjust=1),aspect.ratio=1/2 )  

p<-p + ggtitle("Pathways")

p<- p + theme(plot.title = element_text(size=60,margin = margin(10, 0, 10, 0)))

p

#######gene heatmap

z<-grep("MT-",ser.m.markers$gene)
z<-ser.m.markers[-z,]
z<-z[-(grep("hsa-mir",z$gene)),]
z<-z[-(grep("RP5-",z$gene)),]
x<-foreach(i=0:12,.combine='c') %do% {
  y<-z[z$cluster==i,]
  y<-y[y$p_val_adj<0.0001,]
  b<-y$pct.1/y$pct.2
  y<-y[which(b>0.5),]
  y<-y[which(y$pct.1>0.5),]
  y<-y[order(y$avg_logFC,decreasing=T),]
  y$gene[1:5]
}
DoHeatmap(int.obj,features=x)+scale_fill_distiller(palette="RdYlBu")
#######stem cell GO ananlysis pmid:18511597 
###Get GO genes for stem cell cats, "Nucleus", "RNA binding", and diff cats
nuc<-get_anno_genes("GO:0005634")
nuc<-unique(nuc$gene)
nuc<-unique(nuc$hgnc_symbol)
int_mem<-get_anno_genes("GO:0016021")
int_mem<-unique(int_mem$gene)


x<-foreach(i=1:13,.combine='rbind') %do% {
  x1<-ser.clust[[i]]
  nucs<-length(intersect(x1,nuc))/length(x1)
  mems<-length(intersect(x1,int_mem))/length(x1)
  c(mems,nucs)
}
colnames(x)<-c("Mems","Nuc")
rownames(x)<-0:12
plot(x,pch=19,cex=2,col=gg_color_hue(13))
text(Nuc~Mems,labels=rownames(x),data=x,pos=4,offset=0.5)


############sligshot analysis
z<-int.obj@reductions$umap@cell.embeddings
y<-(as.vector(int.obj@active.ident))
names(y)<-names(int.obj@active.ident)
fb<-slingshot(z,y)

fbp<-slingPseudotime(fb,na=F)
cw <- slingCurveWeights(fb)
assignment <- apply(cw, 1, which.max)
ptAs <- c() 
for(ii in 1:nrow(fbp)) ptAs[ii] <- fbp[ii,assignment[ii]]
t1 <- ptAs[assignment == 1]
cond1<-factor(y[which(assignment==1)])
m1s<-gam(cond1 ~ s(t1), family="quasibinomial")
summary(m1s)

t1 <- ptAs[assignment == 2]
cond1<-factor(y[which(assignment==2)])
m2s<-gam(cond1 ~ s(t1), family="quasibinomial")
summary(m2s)

t1 <- ptAs[assignment == 3]
cond1<-factor(y[which(assignment==3)])
m3s<-gam(cond1 ~ s(t1), family="quasibinomial")
summary(m3s)

fbgr<-graph_from_adjacency_matrix(fb@adjacency,mode="upper")
plot.igraph(fbgr)

########focusing on cluster 10 using ICELL8 data only
obj<-ICELL8.nomac
x<-int.obj@active.ident[grep("cell_",names(int.obj@active.ident))]
obj@active.ident<-x
I8.int.clust.markers<-FindAllMarkers(obj,only.pos=T,
                           likeFC.threshold=0.5)

x<-I8.int.clust.markers[I8.int.clust.markers$p_val_adj<0.01,]
I8.int.clust<-foreach(i=0:12) %do% {x[x$cluster==i,]$gene}
names(I8.int.clust)<-foreach(i=0:12,.combine='c') %do% {paste("Cluster_",i,sep="")}
##############Enrichment analysis
cluster.enr.int.i8<-foreach(i=1:length(I8.int.clust)) %do% {
  Path_BIOCARTA<-hypeR(I8.int.clust[[i]], BIOCARTA, background=23459, fdr=0.01)
  Path_Reactome<-hypeR(I8.int.clust[[i]], REACTOME, background=23459, fdr=0.01)
  Path_KEGG<-hypeR(I8.int.clust[[i]], KEGG, background=23459, fdr=0.01)
  GOBioProc<-hypeR(I8.int.clust[[i]], GOBP, background=23459, fdr=0.01)
  GOMoleFunc<-hypeR(I8.int.clust[[i]], GOMF, background=23459, fdr=0.01)
  GOCellComp<-hypeR(I8.int.clust[[i]], GOCC, background=23459, fdr=0.01)
  list(Path_BIOCARTA,Path_Reactome,Path_KEGG,GOBioProc,GOMoleFunc,GOCellComp)
}
names(cluster.enr.int.i8)<-paste("Cluster",c(0:12),sep="_")
cluster.enr.i8.data<-foreach(i=1:length(I8.int.clust)) %do% {
  Path1<-cluster.enr.int.i8[[i]][[1]]$data
  Path2<-cluster.enr.int.i8[[i]][[2]]$data
  Path3<-cluster.enr.int.i8[[i]][[3]]$data
  Pathways<-rbind(Path1,Path2,Path3)
  BP<-cluster.enr.int.i8[[i]][[4]]$data
  MF<-cluster.enr.int.i8[[i]][[5]]$data
  CC<-cluster.enr.int.i8[[i]][[6]]$data
  list(Pathways,BP,MF,CC)
}
names(cluster.enr.i8.data)<-paste("Cluster",c(0:12),sep="_")

#########stem cell GO ananlysis pmid:18511597 
x<-foreach(i=1:13,.combine='rbind') %do% {
  x1<-I8.int.clust[[i]]
  nucs<-length(intersect(x1,nuc))/length(x1)
  mems<-length(intersect(x1,int_mem))/length(x1)
  c(mems,nucs)
}
colnames(x)<-c("Mems","Nuc")
rownames(x)<-0:12
plot(x[,1:2])
text(Nuc~Mems,labels=rownames(x),data=x,pos=4)

####slingshot with ICELL8 integrated clustering
obj<-RunUMAP(obj,dims=1:30)
z<-obj@reductions$umap@cell.embeddings
y<-(as.vector(obj@active.ident))
names(y)<-names(obj@active.ident)

fb_int_i8<-slingshot(z,y)

############
x<-I8.int.clust.markers[I8.int.clust.markers$p_val_adj<0.01,]
x<-x[x$cluster==10,]
x<-x[x$avg_logFC>1,]
write.csv(x,"I8_cl10.csv")

####Go to StringDB find interactions from genes in"I8_cl10.csv"
cl10.net<-data.frame(read.delim("string_cl10.tsv",as.is=T))

##### generate a vector of transcription factors using GO categories
TFS<-get_anno_genes("GO:0006366")
RIBO<-get_anno_genes("GO:0005840")
PROT<-get_anno_genes("GO:0000502")
KINE<-get_anno_genes("GO:0016301")
TFS<-setdiff(TFS[,2],RIBO[,2])
TFS<-setdiff(TFS,PROT[,2])
TFS<-setdiff(TFS,KINE[,2])
TFS<-intersect(TFS,nuc)
###prep node attributes and make a graph object
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

#####generate graph object to export to cytoscape for figure
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
  if(length(which(x$X.node1==xr[i]))>0){x<-x[-(which(x$X.node1==xr[i])),]}
  if(length(which(x$node2==xr[i]))>0) {x<-x[-(which(x$node2==xr[i])),]}
}


gph10<-graph_from_data_frame(x,directed=F,vertices=z1)


####generate other figures
StackedVlnPlot(obj,features=c("GATA2","GATA1","MYC","BMI1","REST"))
StackedVlnPlot(obj,features=c("FOXO3","MEF2C","STAT3","CHD7","MYB"))
DotPlot(int.obj,features=c("RHO","SAG","PDE6A",
                           "MERTK","GAS6","MFGE8","ITGAV",
                           "GATA1","MYB","CHEK1","ID2","MKI67","DNMT1",
                           "MT1G","MT1E","MT3",
                           "B2M","VEGFA","FGF2","BMP2","BMP4","IGFBP5",
                           "GNB1","ABCA4","PCP4","MAPT",
                           "DDX5","SF1","SRSF1"),
        cols="RdYlBu")

save.image("Integrated RPESC.RData")
##########
sessionInfo()
#R version 4.0.2 (2020-06-22)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)

#Matrix products: default

#locale:
#  [1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
#[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
#[5] LC_TIME=English_United States.1252    

#attached base packages:
#  [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#  [1] Homo.sapiens_1.3.1                      TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
#[3] OrganismDbi_1.30.0                      GenomicFeatures_1.40.1                 
#[5] GenomicRanges_1.40.0                    GenomeInfoDb_1.24.2                    
#[7] patchwork_1.0.1                         slingshot_1.6.1                        
#[9] princurve_2.1.4                         mgcv_1.8-31                            
#[11] nlme_3.1-148                            igraph_1.2.5                           
#[13] dplyr_1.0.1                             reshape2_1.4.4                         
#[15] GOfuncR_1.8.0                           vioplot_0.3.5                          
#[17] zoo_1.8-8                               sm_2.2-5.6                             
#[19] rrvgo_1.0.1                             LSAfun_0.6.1                           
#[21] rgl_0.100.54                            lsa_0.73.2                             
#[23] SnowballC_0.7.0                         GO.db_3.11.4                           
#[25] riverplot_0.6                           RColorBrewer_1.1-2                     
#[27] Seurat_3.2.0                            ggplot2_3.3.2                          
#[29] pCalibrate_0.2-1                        MCMCpack_1.4-8                         
#[31] MASS_7.3-51.6                           coda_0.19-3                            
#[33] exact2x2_1.6.5                          exactci_1.3-3                          
#[35] ssanv_1.1                               org.Hs.eg.db_3.11.4                    
#[37] AnnotationDbi_1.50.3                    IRanges_2.22.2                         
#[39] S4Vectors_0.26.1                        Biobase_2.48.0                         
#[41] BiocGenerics_0.34.0                     doParallel_1.0.15                      
#[43] iterators_1.0.12                        foreach_1.5.0                          
#[45] hypeR_1.4.0                             UpSetR_1.4.0

#loaded via a namespace (and not attached):
#  [1] rappdirs_0.3.1              SparseM_1.78                rtracklayer_1.48.0          cellTree_1.18.0            
#[5] tidyr_1.1.1                 bit64_4.0.2                 knitr_1.29                  irlba_2.3.3                
#[9] DelayedArray_0.14.1         data.table_1.13.0           rpart_4.1-15                RCurl_1.98-1.2             
#[13] generics_0.0.2              cowplot_1.0.0               RSQLite_2.2.0               RANN_2.6.1                 
#[17] future_1.18.0               bit_4.0.3                   spatstat.data_1.4-3         webshot_0.5.2              
#[21] xml2_1.3.2                  httpuv_1.5.4                SummarizedExperiment_1.18.2 assertthat_0.2.1           
#[25] xfun_0.16                   hms_0.5.3                   evaluate_0.14               promises_1.1.1             
#[29] progress_1.2.2              caTools_1.18.0              dbplyr_1.4.4                DBI_1.1.0                  
#[33] htmlwidgets_1.5.1           mcmc_0.9-7                  purrr_0.3.4                 ellipsis_0.3.1             
#[37] RSpectra_0.16-0             crosstalk_1.1.0.1           gridBase_0.4-7              biomaRt_2.44.1             
#[41] deldir_0.1-28               vctrs_0.3.2                 SingleCellExperiment_1.10.1 quantreg_5.61              
#[45] ROCR_1.0-11                 abind_1.4-5                 withr_2.2.0                 ggforce_0.3.2              
#[49] treemap_2.4-2               sctransform_0.2.1           GenomicAlignments_1.24.0    prettyunits_1.1.1          
##[53] goftest_1.2-2               cluster_2.1.0               ape_5.4                     lazyeval_0.2.2             
#[57] crayon_1.3.4                pkgconfig_2.0.3             slam_0.1-47                 labeling_0.3               
#[61] tweenr_1.0.1                wordcloud_2.6               rlang_0.4.7                 globals_0.12.5             
#[65] lifecycle_0.2.0             miniUI_0.1.1.1              MatrixModels_0.4-1          BiocFileCache_1.12.1       
#[69] rsvd_1.0.3                  tcltk_4.0.2                 polyclip_1.10-0             matrixStats_0.56.0         
#[73] lmtest_0.9-37               graph_1.66.0                Matrix_1.2-18               ggridges_0.5.2             
#[77] pheatmap_1.0.12             png_0.1-7                   viridisLite_0.3.0           bitops_1.0-6               
#[81] KernSmooth_2.23-17          visNetwork_2.0.9            Biostrings_2.56.0           blob_1.2.1                 
#[85] stringr_1.4.0               manipulateWidget_0.10.1     readr_1.3.1                 mapplots_1.5.1             
#[89] scales_1.1.1                memoise_1.1.0               magrittr_1.5                plyr_1.8.6                 
#[93] ica_1.0-2                   gplots_3.0.4                gdata_2.18.0                zlibbioc_1.34.0            
#[97] compiler_4.0.2              kableExtra_1.1.0            fitdistrplus_1.1-1          Rsamtools_2.4.0            
#[101] XVector_0.28.0              listenv_0.8.0               pbapply_1.4-2               tidyselect_1.1.0           
#[105] stringi_1.4.6               yaml_2.2.1                  GOSemSim_2.14.1             askpass_1.1                
#[109] ggrepel_0.8.2               grid_4.0.2                  tools_4.0.2                 future.apply_1.6.0         
#[113] rstudioapi_0.11             gridExtra_2.3               farver_2.0.3                Rtsne_0.15                 
#[117] digest_0.6.25               BiocManager_1.30.10         shiny_1.5.0                 Rcpp_1.0.5                 
#[121] later_1.1.0.1               RcppAnnoy_0.0.16            httr_1.4.2                  colorspace_1.4-1           
#[125] rvest_0.3.6                 XML_3.99-0.5                tensor_1.5                  reticulate_1.16            
#[129] splines_4.0.2               uwot_0.1.8                  RBGL_1.64.0                 conquer_1.0.1              
#[133] spatstat.utils_1.17-0       plotly_4.9.2.1              xtable_1.8-4                jsonlite_1.7.0             
#[137] spatstat_1.64-1             NLP_0.2-0                   R6_2.4.1                    tm_0.7-7                   
#[141] pillar_1.4.6                htmltools_0.5.0             mime_0.9                    glue_1.4.1                 
#[145] fastmap_1.0.1               reactable_0.2.0             BiocParallel_1.22.0         codetools_0.2-16           
#[149] lattice_0.20-41             tibble_3.0.3                curl_4.3                    leiden_0.3.3               
#[153] gtools_3.8.2                zip_2.0.4                   openxlsx_4.1.5              openssl_1.4.2              
#[157] survival_3.1-12             limma_3.44.3                rmarkdown_2.3               munsell_0.5.0              
#[161] GenomeInfoDbData_1.2.3      gtable_0.3.0                msigdbr_7.1.1               


