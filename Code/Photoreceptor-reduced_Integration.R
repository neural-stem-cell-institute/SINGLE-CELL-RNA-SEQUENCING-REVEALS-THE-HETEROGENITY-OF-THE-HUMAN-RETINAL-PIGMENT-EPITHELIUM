# This code generates a reduced integrated RPE Seurat object from ICELL8 and 10x datasets lacking
# photoreceptor genes and makes a riverplot for the redistribution of cells in the reduced dataset
# Code developed by Farhad Farjood

library(doParallel)
library(Seurat)
library(RColorBrewer)
library(riverplot)
library(foreach)
library(rrvgo)
library(dplyr)
library(corrplot)
library(ggplot2)

options(future.globals.maxSize=2000 * 1024 ^ 2)

# Load ICELL8 abd 10x datasets

load("/Users/farhad/Desktop/R/Paper_figures/ICELL8-nomac.RData")
load("/Users/farhad/Desktop/R/Paper_figures/10X.RData")

# Remove photoreceptor genes
list.retina<-c("RHO","PDE6A","PDE6B","PDE6G","SAG","RIMS2","MAP2","SCG2","RP1","PDC","DCT","CNGA1","NR2E3",
               "ARR3","PDE6H","OPN1SW","OPN1MW","OPN1MW2","OPN1LW","MAP1B","AIPL1","NRL","VSX1","VSX2",
               "RCVRN","GNGT1","GNGT2","GRK1","GUCY2D","GUCY2F","RGS9","PRPH2","NEUROD1","CNGB1",
               "RPGRIP1","EYS","PPEF2","GUCA1B","PPFIA2","PROM1","IMPG1","CPLX4","RDH8","GNAT1","GNB1",
               "ROM1","KCNV2","MPP4","RBP3","GUCA1A","SAMD7","YPEL2","AHI1","KCNB1","IMPG2","RDH12","PPP3R1",
               "RAX2")

ICELL8.nomac<-subset(ICELL8.nomac, features = setdiff(rownames(ICELL8.nomac), list.retina))
RNA.10x<-subset(RNA.10x, features = setdiff(rownames(RNA.10x), list.retina))

# Integrate datasets
obj.list<-list()

obj.list[[1]]<-SCTransform(ICELL8.nomac,verbose=F)
obj.list[[2]]<-SCTransform(RNA.10x,verbose=F)

rm(ICELL8.nomac,RNA.10x,markers.10x,markers.nomac)

m1<-SelectIntegrationFeatures(object.list = obj.list, nfeatures = 4000)
m<-PrepSCTIntegration(object.list = obj.list, anchor.features = m1 , 
                      verbose = FALSE)
m4<-FindIntegrationAnchors(object.list = m, normalization.method = "SCT",
                           anchor.features = m1, verbose = FALSE,reduction="cca")
Int.noretina<-IntegrateData(anchorset = m4, normalization.method = "SCT", 
                       verbose = FALSE)
rm(m,m4,obj.list)

Int.noretina<-RunPCA(Int.noretina,verbose=F)
Int.noretina<-RunUMAP(Int.noretina,dims=1:30)

Int.noretina<- FindNeighbors(Int.noretina, dims = 1:10)
Int.noretina <- FindClusters(Int.noretina, resolution = 0.51)

DefaultAssay(Int.noretina)<-"RNA"
Int.noretina<-NormalizeData(Int.noretina,verbose=F)
Int.noretina<-ScaleData(Int.noretina,features=rownames(Int.noretina@assays$RNA@counts))

DimPlot(Int.noretina,group.by="tech",label=F, pt.size = .1, shuffle=T)
DimPlot(Int.noretina,group.by="seurat_clusters",label=T, label.size = 6, pt.size = .1)

DefaultAssay(Int.noretina)<-"RNA"

# Function to make a riverplot
quickriver2<-function(obj1,obj2) {
  edges<-foreach(i=0:max(as.numeric(levels(obj1@active.ident))),.combine='rbind') %do% {
    x<-table(obj2@active.ident[names(obj1@active.ident[obj1@active.ident==i])])
    x<-x/sum(x)
    cbind(paste("o",rep(i,length(x)),sep=""),paste("n",
                                                   as.numeric(levels(obj2@active.ident)),sep=""),x)
  }
  
  colnames(edges)<-c("N1","N2","Value")
  edges<-as.data.frame(edges)
  edges$Value<-as.numeric(edges$Value)
  nodes<-data.frame(ID=unique(c(as.character(edges$N1),as.character(edges$N2))), stringsAsFactors = F)
  nodes$x<-c(rep(0,max(as.numeric(levels(obj1@active.ident)))+1),
             rep(5,max(as.numeric(levels(obj2@active.ident)))+1))
  
  
  pal1<-paste0(colorRampPalette(c("royalblue1","green3"))(max(as.numeric(levels(obj1@active.ident)))+1),"80")
  pal2<-paste0(colorRampPalette(c("mediumorchid1","goldenrod3"))(max(as.numeric(levels(obj2@active.ident)))+1),"80")
  palette<-c(pal1, pal2)
  
  nodes$col<-palette
  styles<-foreach(i=1:length(palette)) %do% {
    list(col=palette[i],textcol="black",srt=0,nodestyle="regular",lty=0,nodewidth=1,
         yscale=0.5)}
  names(styles) = nodes$ID
  
  edges<-edges[edges$Value>0,]
  rp<-makeRiver(nodes,edges,node_style=styles)
  
  riverplot(rp,plot_area=1,add_mid_points=T, gravity="center")
}

# Load the original integrated object and make riverplot
Integrated_Primary_RPE <- readRDS("Integrated_Primary_RPE.rds")
quickriver2(Integrated_Primary_RPE,Int.noretina)

# session info:
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
#   [1] ggplot2_3.3.3             corrplot_0.90             dplyr_1.0.3              
# [4] rrvgo_1.2.0               riverplot_0.10            RColorBrewer_1.1-2       
# [7] doParallel_1.0.16         iterators_1.0.13          foreach_1.5.1            
# [10] SeuratObject_4.0.0        Seurat_4.0.1              EnsDb.Hsapiens.v75_2.99.0
# [13] ensembldb_2.14.0          AnnotationFilter_1.14.0   GenomicFeatures_1.42.2   
# [16] AnnotationDbi_1.52.0      Biobase_2.50.0            GenomicRanges_1.42.0     
# [19] GenomeInfoDb_1.26.2       IRanges_2.24.1            S4Vectors_0.28.1         
# [22] BiocGenerics_0.36.0      
# 
# loaded via a namespace (and not attached):
#   [1] BiocFileCache_1.14.0        plyr_1.8.6                  igraph_1.2.6               
# [4] lazyeval_0.2.2              splines_4.0.3               BiocParallel_1.24.1        
# [7] listenv_0.8.0               scattermore_0.7             gridBase_0.4-7             
# [10] digest_0.6.27               GOSemSim_2.16.1             htmltools_0.5.1.1          
# [13] GO.db_3.12.1                magrittr_2.0.1              memoise_2.0.0              
# [16] tm_0.7-8                    tensor_1.5                  cluster_2.1.0              
# [19] ROCR_1.0-11                 globals_0.14.0              Biostrings_2.58.0          
# [22] wordcloud_2.6               matrixStats_0.58.0          askpass_1.1                
# [25] spatstat.sparse_2.0-0       prettyunits_1.1.1           colorspace_2.0-0           
# [28] treemap_2.4-2               blob_1.2.1                  rappdirs_0.3.2             
# [31] ggrepel_0.9.1               xfun_0.20                   crayon_1.3.4               
# [34] RCurl_1.98-1.2              jsonlite_1.7.2              spatstat.data_2.1-0        
# [37] survival_3.2-7              zoo_1.8-8                   glue_1.4.2                 
# [40] polyclip_1.10-0             gtable_0.3.0                zlibbioc_1.36.0            
# [43] XVector_0.30.0              leiden_0.3.7                DelayedArray_0.16.3        
# [46] future.apply_1.7.0          abind_1.4-5                 scales_1.1.1               
# [49] pheatmap_1.0.12             DBI_1.1.1                   miniUI_0.1.1.1             
# [52] Rcpp_1.0.6                  viridisLite_0.3.0           xtable_1.8-4               
# [55] progress_1.2.2              reticulate_1.18             spatstat.core_2.2-0        
# [58] bit_4.0.4                   htmlwidgets_1.5.3           httr_1.4.2                 
# [61] ellipsis_0.3.1              ica_1.0-2                   pkgconfig_2.0.3            
# [64] XML_3.99-0.6                uwot_0.1.10                 dbplyr_2.1.0               
# [67] deldir_0.2-9                tidyselect_1.1.0            rlang_0.4.10               
# [70] reshape2_1.4.4              later_1.1.0.1               munsell_0.5.0              
# [73] tools_4.0.3                 cachem_1.0.1                generics_0.1.0             
# [76] RSQLite_2.2.3               ggridges_0.5.3              stringr_1.4.0              
# [79] fastmap_1.1.0               goftest_1.2-2               bit64_4.0.5                
# [82] fitdistrplus_1.1-3          purrr_0.3.4                 RANN_2.6.1                 
# [85] nlme_3.1-149                pbapply_1.4-3               future_1.21.0              
# [88] mime_0.9                    slam_0.1-48                 xml2_1.3.2                 
# [91] biomaRt_2.46.3              compiler_4.0.3              rstudioapi_0.13            
# [94] plotly_4.9.3                curl_4.3                    png_0.1-7                  
# [97] spatstat.utils_2.2-0        tibble_3.0.5                stringi_1.5.3              
# [100] lattice_0.20-41             ProtGenerics_1.22.0         Matrix_1.2-18              
# [103] vctrs_0.3.6                 pillar_1.4.7                lifecycle_1.0.0            
# [106] spatstat.geom_2.2-0         lmtest_0.9-38               RcppAnnoy_0.0.18           
# [109] data.table_1.13.6           cowplot_1.1.1               bitops_1.0-6               
# [112] irlba_2.3.3                 httpuv_1.5.5                patchwork_1.1.1            
# [115] rtracklayer_1.50.0          R6_2.5.0                    promises_1.1.1             
# [118] KernSmooth_2.23-17          gridExtra_2.3               parallelly_1.23.0          
# [121] codetools_0.2-16            MASS_7.3-53                 assertthat_0.2.1           
# [124] SummarizedExperiment_1.20.0 openssl_1.4.3               withr_2.4.1                
# [127] GenomicAlignments_1.26.0    sctransform_0.3.2           Rsamtools_2.6.0            
# [130] GenomeInfoDbData_1.2.4      mgcv_1.8-33                 hms_1.0.0                  
# [133] rpart_4.1-15                grid_4.0.3                  tidyr_1.1.2                
# [136] MatrixGenerics_1.2.1        Rtsne_0.15                  NLP_0.2-1                  
# [139] shiny_1.6.0                 tinytex_0.30 
