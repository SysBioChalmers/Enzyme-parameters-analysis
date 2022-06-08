library(ggplot2)
library(viridis)
library(viridisLite)
#Get cumulative distributions for Kcats in Sce model 
if (exists("RStudio.Version")){
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  setwd(getSrcDirectory()[1])
}
kcats_new <- read.csv(file = '../results/newKcats.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
kcats_old <- read.csv(file = '../results/oldKcats.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
kcats_all <- read.csv(file = '../results/kineticData.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
kcats_sce <- kcats_all[grep('saccharomyces cerevisiae',kcats_all[,4]),]
kcats_fungi <- kcats_all[grep(';fungi;',kcats_all[,4]),]

kcats_sce <- as.matrix(kcats_sce[,5])
kcats_fungi <- as.matrix(kcats_fungi[,5])
#kcats_new <- as.matrix((kcats_new/3600))
#kcats_old <- as.matrix((kcats_old/3600))

kcats_new <- as.matrix(unique(kcats_new/3600))
kcats_old <- as.matrix(unique(kcats_old/3600))

df <- data.frame(x = c(kcats_sce,kcats_fungi,kcats_old,kcats_new),extension = factor(rep(c('Sce (BRENDA)','Fungi (BRENDA)','ecYeast7 (GECKO)','ecYeast7 (GECKO 2)'), c(length(kcats_sce),length(kcats_fungi),length(kcats_old),length(kcats_new)))))
df <- df[order(df$x), ]
df$x <- as.numeric(df$x)
df$ecdf <- ave(df$x, df$extension, FUN=function(x) seq_along(x)/length(x))
cdfP <- ggplot(df, aes(x, ecdf, color = extension)) + geom_line(size=1.5) + 
        theme_bw(base_size = 26)+scale_x_log10(limits = c(0.1,1000)) +
        labs(y = 'Relative frequency', x='Kcat [1/s]',color = '') + 
        geom_hline(data=df,aes(yintercept=0.5),linetype="dashed",color = "black", size=0.5) +
        scale_color_manual(values = cividis(4),name = "Distribution", labels =c('Fungi (BRENDA)','Sce (BRENDA)','ecYeast7 (GECKO)','ecYeast7 (GECKO 2)'))
plotName <- paste('../results/Kcats_cumDist_ecYeast.eps',sep='')
ggsave(cdfP, file=plotName, device="eps",width = 10, height = 6)
median(kcats_sce)
#Get spider plots comparing WCs between GECKO 1 and GECKO2
library(reshape)
library(ggbiplot)
library(ggplot2)
library(fmsb)
maxLim <- 3000
minLim <- 0
WCtable <- read.csv(file = '../results/WC_comparison.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
rownames(WCtable) <- c('WC0','WC1','WC2','WC3')
colnames(WCtable) <- c('GECKOv2','GECKO')
df <- t(WCtable)

df <- rbind(rep(minLim,ncol(df)),df)
df <- rbind(rep(maxLim,ncol(df)),df)
# Color vector
colors_border = c(rgb(0.8,0.6,0,0.8), rgb(0.1,0,0.8,0.8))
colors_in     = c(rgb(0.8,0.6,0,0.2), rgb(0.1,0,0.8,0.2))
newData <- as.data.frame(df)
newData <- (as.data.frame(df))
plotName <- paste('../results/wildCardsComparison.png',sep='')
png(plotName,width = 950, height = 900)
radarchart( newData  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="black", caxislabels=seq(minLim,maxLim,(maxLim-minLim)/4),cglwd=2,
            #custom labels
            vlcex=4, calcex = 3 
)
dev.off()
# Add a legend
legendStr <- c('GECKOV2','GECKO')
#legendStr <- c('PredictedEt','ExperimentalEt','PredictedGlc','ExperimentalGlc')
legend(x=1, y=1.3, legend =legendStr, bty = "n", pch=20 , col=colors_in , text.col = "black", cex=1.2, pt.cex=3)
WC_new <- read.csv(file = '../results/WC_GECKOv2.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)
WC_new$conditions <- c('1','2','3','4','5','6')
p<-ggplot(data=WC_new, aes(x=conditions,y=WC0)) +
  geom_bar(stat="identity",color="black", fill=cividis(1)) +
  theme_classic(base_size = 32) +labs(y = 'Number of matches',x='')
plotName <- paste('../results/Kcat_categories_ecYeast.eps',sep='')
ggsave(p, file=plotName, device="eps",width = 6, height = 4.5)