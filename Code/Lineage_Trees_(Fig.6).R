####make plots for lineage trees (figure 6) 

library(foreach)
library(ggplot2)

trees<-data.frame(read.csv("result_trees.csv",as.is=T))
trees$ID<-as.character(trees$X)
trees<-trees[,-1]
trees$Asymmetric<-as.character(trees$Asymmetric)

y<-expand.grid(y=1:10,x=1:16)
z<-round(table(as.factor(trees$Generations)) * ((10*16)/(length(as.factor(trees$Generations)))))

y$Generations<-factor(rep(names(z), z)) 
y$Cells<-trees[order(trees$Generations),]$Cells_produced



p<-ggplot(y, aes(x=x, y=y, fill=Generations))
p<-p + geom_tile(color = "black", size = 0.5)
p<-p + geom_text(aes(label=Cells))
p<-p + scale_x_continuous(expand = c(0, 0)) 
p<-p +scale_y_continuous(expand = c(0, 0), trans = 'reverse')
p<-p + theme(plot.title = element_text(size = rel(1.2)),
      axis.text = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right")
p<-p + ggtitle("Lineage Trees")


p

z<-trees[trees$Generations>0,]
z<-z[-4,]
z$Ave_CD<-z$Ave_CD/6
z$Slowest_CD<-z$Slowest_CD/6
z$Fastest_CD<-z$Fastest_CD/6


p<-ggplot(z, aes(x=Generations, y=Ave_CD))
p<-p + geom_point(aes(color=Asymmetric), size=3)
p<-p + scale_color_manual(values=c("gray58","green3"))
p<-p + scale_x_continuous(breaks=1:9)
p<-p + theme_bw() + theme(axis.line = element_blank(),        
                          panel.grid.minor.x = element_blank(),
                          panel.grid.minor.y = element_blank())


p

#sessionInfo()
#R version 4.0.2 (2020-06-22)
#Platform: x86_64-w64-mingw32/x64 (64-bit)
#Running under: Windows >= 8 x64 (build 9200)

#Matrix products: default

#locale:
#[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
#[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

#attached base packages:
#[1] stats     graphics  grDevices utils     datasets  methods   base     

#other attached packages:
#[1] ggplot2_3.3.2 foreach_1.5.0

#loaded via a namespace (and not attached):
# [1] codetools_0.2-16 withr_2.2.0      dplyr_1.0.1      crayon_1.3.4     grid_4.0.2       R6_2.4.1        
# [7] lifecycle_0.2.0  gtable_0.3.0     magrittr_1.5     scales_1.1.1     pillar_1.4.6     rlang_0.4.7     
#[13] rstudioapi_0.11  generics_0.0.2   vctrs_0.3.2      ellipsis_0.3.1   iterators_1.0.12 tools_4.0.2     
#[19] glue_1.4.1       purrr_0.3.4      munsell_0.5.0    yaml_2.2.1       compiler_4.0.2   pkgconfig_2.0.3 
#[25] colorspace_1.4-1 tidyselect_1.1.0 tibble_3.0.3  
