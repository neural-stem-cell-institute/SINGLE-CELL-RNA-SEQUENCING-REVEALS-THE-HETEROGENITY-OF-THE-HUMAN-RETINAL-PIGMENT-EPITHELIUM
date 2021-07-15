# This code generates a Seurat object and performs clustering for
# ICELL8 single cell RNA sequencing data from raw read counts
# Code developed by Farhad Farjood

#load libraries
library(EnsDb.Hsapiens.v75)
library(Seurat)

#function to QC cells and choose ones for further analysis
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

#read in raw data from RData files generated with the mapping pipeline

load("/ICELL8-Donorx.RData")

x<-data.frame(ICELL8-Donorx.RData@assays$data$counts, row.names = rownames(ICELL8-Donorx.RData))
RNA<-x[,grep("_R2",colnames(x))];

colnames(RNA) <- sprintf("Donorx_%d",seq(1:ncol(RNA)))

#QC of cells to choose subset for analysis
z<-apply(RNA,2,sum)
quantile(z)
RNA<-choose.cells(RNA,library.size=c(25000,max(z)))

#set up annotation files

x<-gsub("\\..*","",rownames(RNA))
x <- mapIds( EnsDb.Hsapiens.v75, keys=x, column="SYMBOL", keytype="GENEID", multiVals="first")
x <- ifelse(!duplicated(x) & !is.na(x), x, names(x))

rownames(RNA)<-x

RNA <- CreateSeuratObject(counts = RNA, project = "Donorx", min.cells = 3, min.features = 200)

# store mitochondrial percentage in object meta data
RNA <- PercentageFeatureSet(RNA, pattern = "^MT-", col.name = "percent.mt")

# run sctransform
RNA <- SCTransform(RNA, vars.to.regress = "percent.mt", verbose = FALSE)
RNA <- FindVariableFeatures(RNA, selection.method = "mean.var.plot", nfeatures = 2000)
RNA <- RunPCA(RNA, features = VariableFeatures(RNA), verbose = FALSE)
RNA <- RunUMAP(RNA, dims = 1:30, verbose = FALSE)

RNA <- FindNeighbors(RNA, dims = 1:30, verbose = FALSE)
RNA <- FindClusters(RNA, resolution = 0.5, verbose = FALSE)

DimPlot(RNA, reduction = "umap", label = TRUE) + NoLegend()


#session info:

# R version 4.0.3 (2020-10-10)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Big Sur 10.16
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] SeuratObject_4.0.0        Seurat_4.0.1              EnsDb.Hsapiens.v75_2.99.0
# [4] ensembldb_2.14.0          AnnotationFilter_1.14.0   GenomicFeatures_1.42.2   
# [7] AnnotationDbi_1.52.0      Biobase_2.50.0            GenomicRanges_1.42.0     
# [10] GenomeInfoDb_1.26.2       IRanges_2.24.1            S4Vectors_0.28.1         
# [13] BiocGenerics_0.36.0      
# 
# loaded via a namespace (and not attached):
#   [1] BiocFileCache_1.14.0        plyr_1.8.6                  igraph_1.2.6               
# [4] lazyeval_0.2.2              splines_4.0.3               BiocParallel_1.24.1        
# [7] listenv_0.8.0               scattermore_0.7             ggplot2_3.3.3              
# [10] digest_0.6.27               htmltools_0.5.1.1           magrittr_2.0.1             
# [13] memoise_2.0.0               tensor_1.5                  cluster_2.1.0              
# [16] ROCR_1.0-11                 globals_0.14.0              Biostrings_2.58.0          
# [19] matrixStats_0.58.0          askpass_1.1                 spatstat.sparse_2.0-0      
# [22] prettyunits_1.1.1           colorspace_2.0-0            blob_1.2.1                 
# [25] rappdirs_0.3.2              ggrepel_0.9.1               xfun_0.20                  
# [28] dplyr_1.0.3                 crayon_1.3.4                RCurl_1.98-1.2             
# [31] jsonlite_1.7.2              spatstat.data_2.1-0         survival_3.2-7             
# [34] zoo_1.8-8                   glue_1.4.2                  polyclip_1.10-0            
# [37] gtable_0.3.0                zlibbioc_1.36.0             XVector_0.30.0             
# [40] leiden_0.3.7                DelayedArray_0.16.3         future.apply_1.7.0         
# [43] abind_1.4-5                 scales_1.1.1                DBI_1.1.1                  
# [46] miniUI_0.1.1.1              Rcpp_1.0.6                  viridisLite_0.3.0          
# [49] xtable_1.8-4                progress_1.2.2              reticulate_1.18            
# [52] spatstat.core_2.2-0         bit_4.0.4                   htmlwidgets_1.5.3          
# [55] httr_1.4.2                  RColorBrewer_1.1-2          ellipsis_0.3.1             
# [58] ica_1.0-2                   pkgconfig_2.0.3             XML_3.99-0.6               
# [61] uwot_0.1.10                 dbplyr_2.1.0                deldir_0.2-9               
# [64] tidyselect_1.1.0            rlang_0.4.10                reshape2_1.4.4             
# [67] later_1.1.0.1               munsell_0.5.0               tools_4.0.3                
# [70] cachem_1.0.1                generics_0.1.0              RSQLite_2.2.3              
# [73] ggridges_0.5.3              stringr_1.4.0               fastmap_1.1.0              
# [76] goftest_1.2-2               bit64_4.0.5                 fitdistrplus_1.1-3         
# [79] purrr_0.3.4                 RANN_2.6.1                  nlme_3.1-149               
# [82] pbapply_1.4-3               future_1.21.0               mime_0.9                   
# [85] xml2_1.3.2                  biomaRt_2.46.3              compiler_4.0.3             
# [88] rstudioapi_0.13             plotly_4.9.3                curl_4.3                   
# [91] png_0.1-7                   spatstat.utils_2.2-0        tibble_3.0.5               
# [94] stringi_1.5.3               lattice_0.20-41             ProtGenerics_1.22.0        
# [97] Matrix_1.2-18               vctrs_0.3.6                 pillar_1.4.7               
# [100] lifecycle_1.0.0             spatstat.geom_2.2-0         lmtest_0.9-38              
# [103] RcppAnnoy_0.0.18            data.table_1.13.6           cowplot_1.1.1              
# [106] bitops_1.0-6                irlba_2.3.3                 httpuv_1.5.5               
# [109] patchwork_1.1.1             rtracklayer_1.50.0          R6_2.5.0                   
# [112] promises_1.1.1              KernSmooth_2.23-17          gridExtra_2.3              
# [115] parallelly_1.23.0           codetools_0.2-16            MASS_7.3-53                
# [118] assertthat_0.2.1            SummarizedExperiment_1.20.0 openssl_1.4.3              
# [121] GenomicAlignments_1.26.0    sctransform_0.3.2           Rsamtools_2.6.0            
# [124] GenomeInfoDbData_1.2.4      mgcv_1.8-33                 hms_1.0.0                  
# [127] rpart_4.1-15                grid_4.0.3                  tidyr_1.1.2                
# [130] MatrixGenerics_1.2.1        Rtsne_0.15                  shiny_1.6.0                
# [133] tinytex_0.30