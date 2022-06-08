library(ggplot2)
library(viridis)
library(viridisLite)
library(dplyr)
library(tidyverse)

library(ggplot2)
#Get cumulative distributions for Kcats in Sce model 
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}

dataset <- read.csv(file = '../results/kineticData.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE, row.names = 'V1')
#kcats_sce <- kcats_all[grep('saccharomyces cerevisiae',kcats_all[,3]),]
#kcats_fungi <- kcats_all[grep(';fungi;',kcats_all[,3]),]
kingdom <- c('animals','plants','protists','fungi','bacteria','archaea')


data <- dataset
colnames(data) <- c('ecnumber','substrate','Organism','Kcat','sd','WC','metGroups')

idxs <- grep('Others',data$metGroups)
data <- data[-idxs,]

idxs <- grep('CEM',data$metGroups)
data$metGroups[idxs] <- 'CEM'
data$metGroups[-idxs] <- 'ALM & ISM'
idxs <- grep(';animals;',data$Organism)
data$Organism[idxs] <- 'animals'
idxs <- grep(';fungi;',data$Organism)
data$Organism[idxs] <- 'fungi'
idxs <- grep(';plants;',data$Organism)
data$Organism[idxs] <- 'plants'
idxs <- grep(';protists;',data$Organism)
data$Organism[idxs] <- 'protists'
idxs <- grep(';archaea;',data$Organism)
data$Organism[idxs] <- 'archaea'
idxs <- grep(';bacteria;',data$Organism)
data$Organism[idxs] <- 'bacteria'
idxs <- grep('//*//*',data$Organism)
data <- data[-idxs,]
pValues <- c()

for (cat in kingdom){
  print(cat)
  dist1 <- data$Kcat[intersect(grep(cat,data$Organism),grep('CEM',data$metGroups))]
  dist2 <- data$Kcat[intersect(grep(cat,data$Organism),grep('ALM & ISM',data$metGroups))]
  x <- ks.test(dist2,dist1,alternative = "greater")
  print(as.double(x$p.value[1]))
  pValues <- c(pValues,as.double(x$p.value[1]))
}
pValues
#metGroups <- data$metGroups
#data <- data %>% unite(class, Organism, metGroups)
#data$metGroups <- metGroups
p <- ggplot(data, aes(x=Organism, y=Kcat,fill=metGroups)) + 
  geom_violin() + scale_y_log10(limits=c(1E-6,1E6)) + 
  scale_fill_manual(values = cividis(2,direction=-1)) +
  #stat_summary(fun.y=median, geom="point", size=2, color="red")+
  #stat_summary(fun.data=mean_sdl, mult=1,geom="pointrange", color="red")+
  xlab('') + ylab('Kcats [1/s]') +
  labs(fill="Enzymes") +
  theme_bw(base_size = 28) 
  plotName <- paste('../results/Kcats_metGroups_allOrgs.eps',sep='')
  ggsave(p, file=plotName, device="eps",width = 13, height = 6)
  #png(plotName,width = 1400, height = 800)
  #plot(p)
#dev.off()
