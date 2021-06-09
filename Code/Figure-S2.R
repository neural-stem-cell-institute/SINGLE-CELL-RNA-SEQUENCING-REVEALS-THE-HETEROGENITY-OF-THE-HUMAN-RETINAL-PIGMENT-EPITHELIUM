library(hypeR)
library(doParallel)
library(org.Hs.eg.db)
library(pCalibrate)
library(ggplot2)
library(Seurat)
library(RColorBrewer)
library(UpSetR)
library(liger)
library(dplyr)
library(plyr)

load("/Volumes/Farhad_RPE/R/Paper_figures/ICELL8-nomac.RData")

c.names<-ICELL8.nomac$seurat_clusters
# c.names<-mapvalues(c.names, from = 0:11, to = LETTERS[9:(8+length(levels(c.names)))])
c.names<-mapvalues(c.names, from = 0:11, to = 8:19)
ICELL8.nomac$clusters<-c.names

####loaded code "single.cell.samp.diff.R"
registerDoParallel(makeCluster(3))
SYMBOLS<-unlist(as.list(org.Hs.egSYMBOL))
obj<-ICELL8.nomac
DimPlot(obj, group.by = "clusters", label = TRUE)
nomac.top <- markers.nomac %>% group_by(cluster) %>% top_n(n = 1, wt = avg_logFC)

# DotPlot(obj,features=c(nomac.top$gene,"RPE65","BEST1"),scale.min = 0, 
        # group.by = "seurat_clusters", split.by = "patient", cols = c("coral2","darkcyan","orange2"))
DotPlot(obj,features=c(nomac.top$gene,"RPE65","BEST1"),scale.min = 0, 
        group.by = "clusters")
foreach(i=0:(length(levels(obj$seurat_clusters))-1),.combine='c') %do% {
  x<-markers.nomac[markers.nomac$cluster==i,]
  length(which(x$p_val_adj<0.1))
}


x<-read.csv("/Users/Feri/Downloads/human_RPE_signature_genes_26517551.csv",as.is=T)
y<-intersect(rownames(obj@assays$SCT@scale.data),x[,1])
z<-setdiff(x[,1],y)

test<-probs(obj)
length(which((apply(test[[1]][intersect(rownames(test[[1]]),z),3:11],1,max)>0)))
length(which((apply(test[[1]][intersect(rownames(test[[1]]),x[,1]),3:11],1,max)>0)))

m<-test[[1]][which((apply(test[[1]][intersect(rownames(test[[1]]),x[,1]),3:11],1,max)>0)),]

z<-heatmap(obj@assays$SCT@scale.data[y,])
DoHeatmap(obj,features=rev(y[z$rowInd]), group.by = "clusters")+scale_fill_distiller(palette="RdYlBu")

foreach(i=0:max(as.numeric(levels(obj$seurat_clusters))),.combine='c') %do% {
  length(which(obj@active.ident==i))}

m<-foreach(i=1:dim(m)[2]) %do% {names(which(m[,i]>0))}

names(m)<-paste("Cluster",levels(ICELL8.nomac$clusters),sep=" ")
y<-fromList(m)

upset(y,nsets=12,nintersects = NA,
      matrix.dot.alpha=0.5,order.by="freq",show.numbers=F, set_size.scale_max=T)
###get databases
msigdb_info <- msigdb_download_all(species="Homo sapiens")
BIOCARTA <- msigdb_fetch(msigdb_info, "C2.CP.BIOCARTA")
KEGG     <- msigdb_fetch(msigdb_info, "C2.CP.KEGG")
REACTOME <- msigdb_fetch(msigdb_info, "C2.CP.REACTOME")

paths<-c(BIOCARTA, KEGG, REACTOME)
GOBP <- msigdb_fetch(msigdb_info, "C5.BP")
GOMF <- msigdb_fetch(msigdb_info, "C5.MF")
GOCC <- msigdb_fetch(msigdb_info, "C5.CC")

cluster.genes.10x<-foreach(i=0:max(as.numeric(levels(obj@meta.data$seurat_clusters)))) %do% {
  x<-markers.nomac[markers.nomac$cluster==i,]
  rownames(x[x$p_val_adj<0.1,])
}
names(cluster.genes.10x)<-paste("Cluster",0:11,sep="_")


foreach(i=1:length(cluster.genes.10x),.combine='c') %do% {
  length(cluster.genes.10x[[i]])}


cluster.enr.10x<-foreach(i=1:length(cluster.genes.10x)) %do% {
  
  PATHWAYS<-hypeR(cluster.genes.10x[[i]], paths, bg=23459, pval_cutoff =.05)$data #[,1:7]
  GOBioProc<-hypeR(cluster.genes.10x[[i]], GOBP, bg=23459, fdr_cutoff =05)$data #[,1:7]
  GOMoleFunc<-hypeR(cluster.genes.10x[[i]], GOMF, bg=23459, fdr_cutoff = 05)$data #[,1:7]
  GOCellComp<-hypeR(cluster.genes.10x[[i]], GOCC, bg=23459, fdr_cutoff = 05)$data #[,1:7]
  list(PATHWAYS,GOBioProc,GOMoleFunc,GOCellComp)
}

names(cluster.enr.10x)<-paste("Cluster",LETTERS[seq(9:9+length(m))],sep="_")

enr<-c("Pathways","GO:BP","GO:MF","GO:CC")


x<-foreach(i=1:length(cluster.enr.10x), .combine='rbind') %do% {
  y<-cluster.enr.10x[[i]]

  z<-foreach(b=1:4, .combine='rbind') %do% {
    m<-y[[b]]
    if(dim(m)[1]==0) {} else {m$Type<-enr[b]}
    m
  }
  
  if(dim(z)[1]==0) {} else {z$Cluster<-names(cluster.enr.10x)[i]}
  z
}
write.csv(x,"Supplemental_Table_2.csv",row.names=F)


#####

y<-c("Pathways","GOBP","GOMF","GOCC")


enr.df<-foreach(j=1:length(y),.combine='rbind') %do% {
  
  fin<-foreach(i=1:length(cluster.enr.10x),.combine='rbind') %do% {
    x<-cluster.enr.10x[[i]]
    if(dim(x[[j]])[1]==0){} else {
      x<-data.frame(x[[j]],rep(y[j],dim(x[[j]])[1]),rep(i-1,dim(x[[j]])[1]),rownames(x[[j]]))
    }
    
  }
  colnames(fin)<-c(colnames(cluster.enr.10x[[1]][[j]]),"Category","Cluster","Name")
  fin
}

words<-function(df1,freq,col_pal){
  x<-df1$label
  x<-unlist(strsplit(x,"_"))
  
  # remove numeric characters
  x<-x[is.na(as.numeric(x))]
  
  # remove words smaller than 2 letters
  x<-x[nchar(x)>2]
  
  #remove unwanted words
  
  stopwords<-c("REACTOME","OF","IN","VIA","AND","BIOCARTA","PATHWAY","OF","AND","KEGG","BY","GO","TO",
               "PROCESS", "CELL", "REGULATION", "POSITIVE", "NEGATIVE", "PROTEIN", "FACTOR", "B",
               "RESPONSE","INVOLVED","CELLULAR","ACTIVITY", "SYSTEM", "LIKE", "TYPE","NON","FROM")
  stopwords<-paste("^",stopwords,"$", sep="")
  x<-x[!grepl(paste(stopwords, collapse = "|"), x)]
  
  x<-data.frame(table(x))
  x<-x[order(x$Freq,decreasing=T),]
  
  colnames(x)<-c("word","freq")
  wordcloud(x$word,x$freq,min.freq=freq,random.order=F,random.color=F,max.words=500,
            colors=brewer.pal(8,col_pal), scale=c(.9,.1))
}


x<-enr.df[(enr.df$Category=="GOBP"),]
words(x,5,"PRGn")























