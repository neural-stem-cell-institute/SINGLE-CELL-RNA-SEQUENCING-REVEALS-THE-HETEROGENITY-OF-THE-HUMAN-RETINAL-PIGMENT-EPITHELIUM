# This code replicates the analyses and figures for scRNA-Seq analysis on 10x and ICELL
# scRNA-Seq datasets for RPE cells from 4 adult human donors.
# Code developed by Farhad Farjood

#------------------------------------ sessionInfo ------------------------------------
# R version 3.6.1 (2019-07-05)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Sierra 10.12.6
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] Seurat_3.1.2
# 
# loaded via a namespace (and not attached):
#   [1] tsne_0.1-3          nlme_3.1-143        bitops_1.0-6        RcppAnnoy_0.0.14    RColorBrewer_1.1-2 
# [6] httr_1.4.1          numDeriv_2016.8-1.1 sctransform_0.2.1   tools_3.6.1         R6_2.4.1           
# [11] irlba_2.3.3         KernSmooth_2.23-16  uwot_0.1.5          lazyeval_0.2.2      BiocGenerics_0.32.0
# [16] colorspace_1.4-1    sn_1.5-5            npsurv_0.4-0        gridExtra_2.3       tidyselect_1.0.0   
# [21] mnormt_1.5-6        compiler_3.6.1      Biobase_2.46.0      TFisher_0.2.0       plotly_4.9.1       
# [26] sandwich_2.5-1      labeling_0.3        caTools_1.18.0      scales_1.1.0        lmtest_0.9-37      
# [31] mvtnorm_1.0-12      ggridges_0.5.2      pbapply_1.4-2       rappdirs_0.3.1      stringr_1.4.0      
# [36] digest_0.6.23       R.utils_2.9.2       htmltools_0.4.0     pkgconfig_2.0.3     bibtex_0.4.2.2     
# [41] plotrix_3.7-7       htmlwidgets_1.5.1   rlang_0.4.4         rstudioapi_0.10     farver_2.0.3       
# [46] zoo_1.8-7           jsonlite_1.6.1      ica_1.0-2           gtools_3.8.1        dplyr_0.8.4        
# [51] R.oo_1.23.0         magrittr_1.5        Matrix_1.2-18       Rcpp_1.0.3          munsell_0.5.0      
# [56] ape_5.3             reticulate_1.14     lifecycle_0.1.0     R.methodsS3_1.7.1   stringi_1.4.5      
# [61] multcomp_1.4-12     yaml_2.2.1          gbRd_0.4-11         MASS_7.3-51.5       gplots_3.0.1.2     
# [66] Rtsne_0.15          plyr_1.8.5          grid_3.6.1          parallel_3.6.1      gdata_2.18.0       
# [71] listenv_0.8.0       ggrepel_0.8.1       crayon_1.3.4        lattice_0.20-38     cowplot_1.0.0      
# [76] splines_3.6.1       multtest_2.42.0     SDMTools_1.1-221.2  pillar_1.4.3        igraph_1.2.4.2     
# [81] reshape2_1.4.3      future.apply_1.4.0  codetools_0.2-16    stats4_3.6.1        leiden_0.3.3       
# [86] mutoss_0.1-12       glue_1.3.1          lsei_1.2-0          metap_1.3           RcppParallel_4.4.4 
# [91] data.table_1.12.8   png_0.1-7           vctrs_0.2.2         Rdpack_0.11-1       tidyr_1.0.2        
# [96] gtable_0.3.0        RANN_2.6.1          purrr_0.3.3         future_1.16.0       assertthat_0.2.1   
# [101] ggplot2_3.2.1       rsvd_1.0.2          viridisLite_0.3.0   survival_3.1-8      tibble_2.1.3       
# [106] cluster_2.1.0       globals_0.12.5      fitdistrplus_1.0-14 TH.data_1.0-10      ROCR_1.0-7         


options(future.globals.maxSize = 4000 * 1024^2)

# load seurat objects

load("/Volumes/Farhad_RPE/R/114380/RNA_80.RData")
RNA.80@meta.data$patient <- "318"
RNA.80@meta.data$stage <- "Fresh"
RNA.80@meta.data$tech <- "ICELL8"

load("/Volumes/Farhad_RPE/R/114315/RNA_15.RData")
RNA.15@meta.data$patient <- "319"
RNA.15@meta.data$stage <- "Fresh"
RNA.15@meta.data$tech <- "ICELL8"

load("/Volumes/Farhad_RPE/R/114305/RNA_05.RData")
RNA.05@meta.data$patient <- "322"
RNA.05@meta.data$stage <- "Fresh"
RNA.05@meta.data$tech <- "ICELL8"

# Integrate ICELL8 datasets
ICELL8.list<-c(RNA.80,RNA.15,RNA.05)

for (i in 1:length(ICELL8.list)) {
  ICELL8.list[[i]] <- SCTransform(ICELL8.list[[i]], verbose = FALSE)
}

ICELL8.features <- SelectIntegrationFeatures(object.list = ICELL8.list, nfeatures = 3000)
ICELL8.list <- PrepSCTIntegration(object.list = ICELL8.list, anchor.features = ICELL8.features, verbose = FALSE)


ICELL8.anchors <- FindIntegrationAnchors(object.list = ICELL8.list, normalization.method = "SCT", dims = 1:30,
                                         anchor.features = ICELL8.features, verbose = FALSE)
ICELL8.intg<-IntegrateData(anchorset = ICELL8.anchors, normalization.method = "SCT", dims = 1:30)

ICELL8.intg <- RunPCA(ICELL8.intg, features = VariableFeatures(ICELL8.intg), verbose = FALSE)
ICELL8.intg <- RunUMAP(ICELL8.intg, dims = 1:30, reduction='pca', verbose = FALSE)
ICELL8.intg <- FindNeighbors(ICELL8.intg, dims = 1:30, reduction='pca', verbose = FALSE)
ICELL8.intg <- FindClusters(ICELL8.intg, resolution = 0.8, verbose = FALSE)

# Remove macrophage cluster (cluster 16)

c16<-subset(ICELL8.intg, idents = 16)

RNA80<-subset(RNA.80, cells = setdiff(colnames(RNA.80), colnames(c16)))
RNA15<-subset(RNA.15, cells = setdiff(colnames(RNA.15), colnames(c16)))
RNA05<-subset(RNA.05, cells = setdiff(colnames(RNA.05), colnames(c16)))


ICELL8.list<-c(RNA80,RNA15,RNA05)

for (i in 1:length(ICELL8.list)) {
  ICELL8.list[[i]] <- SCTransform(ICELL8.list[[i]], verbose = FALSE, variable.features.n = 3000)
}

options(future.globals.maxSize = 4000 * 1024^2)

ICELL8.features <- SelectIntegrationFeatures(object.list = ICELL8.list, nfeatures = 3000)
ICELL8.list <- PrepSCTIntegration(object.list = ICELL8.list, anchor.features = ICELL8.features, verbose = FALSE)


ICELL8.anchors <- FindIntegrationAnchors(object.list = ICELL8.list, normalization.method = "SCT", dims = 1:30,
                                         anchor.features = ICELL8.features, verbose = FALSE)

ICELL8.nomac<-IntegrateData(anchorset = ICELL8.anchors, 
                            features.to.integrate = ICELL8.features,
                            normalization.method = "SCT", dims = 1:30)

#DefaultAssay(object = ICELL8.nomac)<-"integrated"

ICELL8.nomac <- RunPCA(ICELL8.nomac, features = VariableFeatures(ICELL8.nomac), verbose = FALSE)
ICELL8.nomac <- RunUMAP(ICELL8.nomac, dims = 1:30, reduction='pca', verbose = FALSE)
ICELL8.nomac <- FindNeighbors(ICELL8.nomac, dims = 1:30, reduction='pca', verbose = FALSE)
ICELL8.nomac <- FindClusters(ICELL8.nomac, resolution = 0.5, verbose = FALSE)

DimPlot(ICELL8.nomac, reduction = "umap", label=TRUE , group.by = c("seurat_clusters"), pt.size = .1)
