#######libraries
library(EnsDb.Hsapiens.v75)
library(Seurat)

####function to QC cells and choose ones for further analysis
choose.cells<-function(mat,library.size=NULL,gene.number=NULL) {
  
  x<-apply(mat,2,sum)
  z<-apply(mat,2,function(x) length(which(x>5)))
  pdf("Libray.vs.gene.pdf",onefile=T,width=11)
  plot(x,z,xlab="Library Size",ylab="Genes with greater than 5 counts")
  hist(x,xlab="Library size")
  hist(z,xlab="Genes with greater than 5 counts")
  dev.off()
  
  if(is.null(library.size)) {
    b<-intersect(which(x>quantile(x,seq(0,1,0.1))[2]),which(x<quantile(x,seq(0,1,0.1))[10]))
    RNA.sub<-mat[,b]
    crit.lib<-c(quantile(x,seq(0,1,0.1))[2],quantile(x,seq(0,1,0.1))[10])
    names(crit.lib)<-c("Min_lib_size","Max_lib_size")
  }else{
    print("Library size entered by user")
    b<-intersect(which(x>library.size[1]),which(x<library.size[2]))
    RNA.sub<-mat[,b]
    crit.lib<-library.size
    names(crit.lib)<-c("Min_lib_size","Max_lib_size")
  }
  if(is.null(gene.number)) {
    z<-apply(RNA.sub,2,function(x) length(which(x>5)))
    y<-log(z)
    RNA.sub<-RNA.sub[,which((median(y)-mad(y)*2)<=y)] 
    x<-apply(RNA.sub,2,sum)
    z<-apply(RNA.sub,2,function(x) length(which(x>5)))
    pdf("Sub.Libray.vs.gene.pdf",onefile=T,width=11)
    plot(x,z,xlab="Library Size",ylab="Genes with greater than 5 counts")
    hist(x,xlab="Library size")
    hist(z,xlab="Genes with greater than 5 counts")
    hist(y,breaks=16,xlab="Genes with greater than 5 counts(log)",main="Before gene size drops are made")
    hist(log(z),breaks=16,xlab="Genes with greater than 5 counts(log)",main="After gene size drops are made")
    dev.off()
    plot(x,z,xlab="Library Size",ylab="Genes with greater than 5 counts")
    hist(x,xlab="Library size")
    hist(z,xlab="Genes with greater than 5 counts")
    hist(y,breaks=16,xlab="Genes with greater than 5 counts(log)",main="Before gene size drops are made")
    hist(log(z),breaks=16,xlab="Genes with greater than 5 counts(log)",main="After gene size drops are made")
    
    crit.gene<-exp(1)^(median(y)-mad(y)*2)
    names(crit.gene)<-"Gene_number_cutoff"
  }else{
    print("Minimal gene size entered by user")
    RNA.sub<-RNA.sub[,which(gene.number<=z)]
    crit.gene<-gene.number
    names(crit.gene)<-"Gene_number_cutoff"
  }
  
  crit<-c(crit.lib,crit.gene)
  write.csv(crit,"Criteria_for_choosing_cells.csv")
  return(RNA.sub)
}

###########read in raw data from RData files generated with the mapping pipeline

load("/Volumes/Farhad_RPE/R/114380/rpe-114380.RData")

x<-data.frame(rpe.114380@assays$data$counts, row.names = rownames(rpe.114380))
RNA<-x[,grep("1_R2",colnames(x))];

colnames(RNA) <- sprintf("cell_80_%d",seq(1:ncol(RNA)))

#####QC of cells to choose subset for analysis
z<-apply(RNA,2,sum)
quantile(z)
RNA<-choose.cells(RNA,library.size=c(25000,max(z)))

######set up annotation files

x<-gsub("\\..*","",rownames(RNA))
x <- mapIds( EnsDb.Hsapiens.v75, keys=x, column="SYMBOL", keytype="GENEID", multiVals="first")
x <- ifelse(!duplicated(x) & !is.na(x), x, names(x))

rownames(RNA)<-x

RNA <- CreateSeuratObject(counts = RNA, project = "114380", min.cells = 3, min.features = 200)

# store mitochondrial percentage in object meta data

RNA <- PercentageFeatureSet(RNA, pattern = "^MT-", col.name = "percent.mt")

# run sctransform, clustering and dimensionality reduction with the Seural workflow
RNA <- SCTransform(RNA, vars.to.regress = "percent.mt", verbose = FALSE)
RNA <- FindVariableFeatures(RNA, selection.method = "vst", nfeatures = 2000)
RNA <- RunPCA(RNA, features = VariableFeatures(RNA), verbose = FALSE)
RNA <- RunUMAP(RNA, dims = 1:30, verbose = FALSE)
RNA <- FindNeighbors(RNA, dims = 1:30, verbose = FALSE)
RNA <- FindClusters(RNA, resolution = .8, verbose = FALSE)

DimPlot(RNA, reduction = "umap", label = TRUE) + NoLegend()
