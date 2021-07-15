# This code performs co-localization analysis on the output of "particle correlation" test
# from the GDSC plugin for ImageJ for confocal images of RPE tissues immunostained for
# combinations of HMGN5 and GATA2 or CMYC and GATA2

# load libraries
library(foreach)
library(doParallel)
library(ggplot2)
library(ggpubr)

registerDoParallel(makeCluster(5))

# Read in the co-localization results
raw.table<-read.csv("/HMGN5-GATA2_Results.csv", header=T)

# Convert the distance to mm
raw.table$x<-raw.table$x*0.0006919
raw.table$y<-raw.table$y*0.0006919

# Trim data to only include only the first 17mm and particle smaller than 200 pixels
all.table<-raw.table[raw.table$x<17 & raw.table$N<200,]

# Define 17 bins (1mm each)
breaks<-17
cuts <- cut(all.table$x,breaks = breaks, include.lowest=FALSE,)

# Find positively stained nuclei. A threshold of >50% for pixel overlap between
# DAPI staining and protrein signal was used
exp.results<-foreach(i = 1:length(levels(cuts)),.combine=rbind) %dopar% {
  nucs<-all.table[cuts==levels(cuts)[i],]
  nucs$bin<-paste("bin",as.character(i))
  
  all.nucs<-nrow(nucs)
  HMGN5.nucs<-length(which((nucs$Sum2/255>nucs$N/2)))
  GATA2.nucs<-length(which((nucs$Sum1/255>nucs$N/2)))
  coexp<-length(which(nucs$Sum1/255>nucs$N/2 & nucs$Sum2/255>nucs$N/2))
  
  exp.data<-c(all.nucs,HMGN5.nucs,GATA2.nucs,coexp)
}

# Calculate percentages and make a distribution histogram
co.results<-cbind (exp.results[,4]/exp.results[,2],
                   exp.results[,4]/exp.results[,3])

exp.results=100*exp.results/exp.results[,1]
co.results<-as.data.frame(cbind(exp.results,co.results*100))

colnames(co.results)<-c("all","hmgn5/all","gata2/all","co/all","co/hmgn5","co/gata2")

meta<-c(rep("all",breaks),rep("hmgn5/all",breaks),rep("gata2/all",breaks),rep("co/all",breaks),rep("co/hmgn5",breaks),rep("co/gata2",breaks))
plot.data<-data.frame(counts=c(as.matrix(co.results)),x=rep(1:breaks,6),meta=meta)

plotdata<-plot.data[breaks+1:(3*breaks),]
plotdata$facetmeta<-factor(plotdata$meta, levels=c("hmgn5/all","gata2/all","co/all"))

GH<-ggplot(plotdata, aes(x,counts, fill=meta, color=meta)) +
  geom_bar(stat = "identity",width = 0.5)+
  facet_grid(~ facetmeta)+
  theme_classic()+
  
  labs(title = "GATA2 and HMGN5 Co-expression",
       x="Distance from optic nerve (mm)",
       y="Positive cell count")+
  scale_y_continuous(expand=c(0,0),limits = c(0,50))+
  scale_x_continuous(breaks=seq(0, 25, by = 5),limits = c(0,20))

GH
