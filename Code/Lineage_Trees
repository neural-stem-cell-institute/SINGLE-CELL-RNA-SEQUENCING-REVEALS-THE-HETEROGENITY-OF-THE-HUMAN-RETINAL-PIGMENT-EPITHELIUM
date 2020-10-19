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


