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

#sessionInfo()
#R version 4.0.2 (2020-06-22)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)

#Matrix products: default

#locale:
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
#[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

#attached base packages:
#[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
# [1] mgcv_1.8-31          nlme_3.1-148         igraph_1.2.5         dplyr_1.0.1          reshape2_1.4.4      
# [6] GOfuncR_1.8.0        vioplot_0.3.5        zoo_1.8-8            sm_2.2-5.6           rrvgo_1.0.1         
#[11] LSAfun_0.6.1         rgl_0.100.54         lsa_0.73.2           SnowballC_0.7.0      GO.db_3.11.4        
#[16] riverplot_0.6        wordcloud_2.6        RColorBrewer_1.1-2   Seurat_3.2.0         ggplot2_3.3.2       
#[21] pCalibrate_0.2-1     MCMCpack_1.4-8       MASS_7.3-51.6        coda_0.19-3          exact2x2_1.6.5      
#[26] exactci_1.3-3        ssanv_1.1            org.Hs.eg.db_3.11.4  AnnotationDbi_1.50.3 IRanges_2.22.2      
#[31] S4Vectors_0.26.1     Biobase_2.48.0       BiocGenerics_0.34.0  doParallel_1.0.15    iterators_1.0.12    
#[36] foreach_1.5.0        hypeR_1.4.0         

#loaded via a namespace (and not attached):
 # [1] reticulate_1.16         tidyselect_1.1.0        RSQLite_2.2.0           htmlwidgets_1.5.1      
 # [5] grid_4.0.2              Rtsne_0.15              munsell_0.5.0           codetools_0.2-16       
 # [9] ica_1.0-2               future_1.18.0           miniUI_0.1.1.1          withr_2.2.0            
 #[13] colorspace_1.4-1        GOSemSim_2.14.1         NLP_0.2-0               knitr_1.29             
 #[17] rstudioapi_0.11         ROCR_1.0-11             tensor_1.5              listenv_0.8.0          
 #[21] slam_0.1-47             GenomeInfoDbData_1.2.3  polyclip_1.10-0         bit64_4.0.2            
 #[25] farver_2.0.3            pheatmap_1.0.12         vctrs_0.3.2             generics_0.0.2         
 #[29] xfun_0.16               R6_2.4.1                GenomeInfoDb_1.24.2     rsvd_1.0.3             
 #[33] msigdbr_7.1.1           manipulateWidget_0.10.1 bitops_1.0-6            spatstat.utils_1.17-0  
 #[37] promises_1.1.1          scales_1.1.1            gtable_0.3.0            globals_0.12.5         
 #[41] conquer_1.0.1           goftest_1.2-2           mcmc_0.9-7              rlang_0.4.7            
 #[45] MatrixModels_0.4-1      splines_4.0.2           lazyeval_0.2.2          yaml_2.2.1             
 #[49] abind_1.4-5             crosstalk_1.1.0.1       httpuv_1.5.4            tools_4.0.2            
 #[53] tcltk_4.0.2             gridBase_0.4-7          ellipsis_0.3.1          kableExtra_1.1.0       
 #[57] ggridges_0.5.2          Rcpp_1.0.5              plyr_1.8.6              zlibbioc_1.34.0        
 #[61] visNetwork_2.0.9        purrr_0.3.4             RCurl_1.98-1.2          rpart_4.1-15           
 #[65] deldir_0.1-28           pbapply_1.4-2           cowplot_1.0.0           mapplots_1.5.1         
 #[69] ggrepel_0.8.2           cluster_2.1.0           magrittr_1.5            data.table_1.13.0      
 #[73] openxlsx_4.1.5          SparseM_1.78            lmtest_0.9-37           RANN_2.6.1             
 #[77] fitdistrplus_1.1-1      matrixStats_0.56.0      hms_0.5.3               patchwork_1.0.1        
 #[81] mime_0.9                reactable_0.2.0         evaluate_0.14           xtable_1.8-4           
 #[85] gridExtra_2.3           compiler_4.0.2          tibble_3.0.3            KernSmooth_2.23-17     
 #[89] crayon_1.3.4            htmltools_0.5.0         later_1.1.0.1           tidyr_1.1.1            
 #[93] DBI_1.1.0               tweenr_1.0.1            rappdirs_0.3.1          Matrix_1.2-18          
 #[97] readr_1.3.1             GenomicRanges_1.40.0    pkgconfig_2.0.3         plotly_4.9.2.1         
#[101] xml2_1.3.2              XVector_0.28.0          webshot_0.5.2           rvest_0.3.6            
#[105] stringr_1.4.0           digest_0.6.25           sctransform_0.2.1       RcppAnnoy_0.0.16       
#[109] spatstat.data_1.4-3     tm_0.7-7                rmarkdown_2.3           leiden_0.3.3           
#[113] uwot_0.1.8              shiny_1.5.0             gtools_3.8.2            quantreg_5.61          
#[117] lifecycle_0.2.0         jsonlite_1.7.0          viridisLite_0.3.0       pillar_1.4.6           
#[121] lattice_0.20-41         fastmap_1.0.1           httr_1.4.2              survival_3.1-12        
#[125] treemap_2.4-2           glue_1.4.1              zip_2.0.4               spatstat_1.64-1        
#[129] png_0.1-7               bit_4.0.3               ggforce_0.3.2           stringi_1.4.6          
#[133] blob_1.2.1              memoise_1.1.0           irlba_2.3.3             future.apply_1.6.0     
#[137] ape_5.4  






