#!/usr/bin/env Rscript
# Functions for plotting the processed data in this study
# Ivan Domenzain

## @knitr getBarPlot
getBarPlot <- function(dataset,BinWidth,xLabel,yLabel,fontSize) {
  #Remove NaN and Inf values
  to.keep <- !(is.na(dataset[,1])|dataset==Inf)
  dataset <- as.data.frame(dataset[to.keep,])
  #plot data
  bars    <- ggplot (dataset, aes(x=dataset[,1])) + 
             geom_histogram(binwidth=BinWidth,color="black", fill=cividis(1)) + 
             labs(y = yLabel, x=xLabel) + theme_classic(base_size = fontSize)
  plot(bars)
}  

## @knitr getCumulativeDist
getCumulativeDist <- function(dataset,xLabel,yLabel,xLimit,fontSize,xIntercept) {
  nargin <- length(as.list(match.call())) -1
  #Remove NaN and Inf values
  to.keep     <- !(is.na(dataset[,1])|dataset==Inf)
  dataset     <- as.data.frame(dataset[to.keep,])
  cdfP        <- ggplot (data=dataset, aes(x= dataset[,1])) + 
                 stat_ecdf(geom = "step",color = cividis(1)) + 
                 labs(y = yLabel, x=xLabel) +xlim(0,xLimit) +
                 theme_classic(base_size = fontSize)
  if (nargin>5){
    cdfP <- cdfP + geom_vline(aes(xintercept=xIntercept),linetype='dotdash')
  }
  plot(cdfP)
}

## @knitr multipleCumDist
multipleCumDist <- function(vectorList,labelStrs,xLabel,yLabel,x_limits,classes,fontSize){
  V <- vectorList
  colores <- cividis(length(vectorList))
  if (length(vectorList)==6){
    df <- data.frame(x = c(V[[1]],V[[2]],V[[3]],V[[4]],V[[5]],V[[6]]),extension = (rep(colores, c(length(V[[1]]),length(V[[2]]),length(V[[3]]),length(V[[4]]),length(V[[5]]),length(V[[6]])))))
  }
  if (length(vectorList)==4){
    df <- data.frame(x = c(V[[1]],V[[2]],V[[3]],V[[4]]),extension = (rep(colores, c(length(V[[1]]),length(V[[2]]),length(V[[3]]),length(V[[4]])))))
  }
  df      <- df[order(df$x), ]
  df$ecdf <- ave(df$x, df$extension, FUN=function(x) seq_along(x)/length(x))
  cdfP    <- ggplot(df, aes(x, ecdf, colour = extension)) + geom_line() + 
             scale_colour_manual(name=classes,labels=labelStrs,values=colores) +
             labs(y = yLabel, x=xLabel) + theme_classic(base_size = 1.75*fontSize)+
             scale_x_log10(limits = x_limits)
  plot(cdfP)
}

## @knitr getPieChart
getPieChart <- function(dataset,column,fontSize,colors,titleStr){
  classes <- unique(dataset[,column])
  nargin  <- length(as.list(match.call())) -1
  if (nargin<4){colors <- factor(classes)}
  N_elements    <- c()
  prop <- c()
  for (element in classes){
    entries    <- length(which(dataset[,column]==element))
    N_elements <- c(N_elements,entries)
    percentage <- entries*100/nrow(dataset)
    prop       <- c(prop,percentage)
  }
  prop <- round(prop, digits = 2)
  df   <- data.frame(classes,N_elements,prop,stringsAsFactors = FALSE)
  # Add label position
  df <- df %>%arrange(desc(classes)) %>%mutate(lab.ypos = cumsum(prop) - 0.5*prop)
  #Create plot object
  colores <- cividis(length(classes))
  if (length(colores)>6){colores <- rev(colores)}
  p <- ggplot(df, aes(x = '', y = prop, fill = classes)) +
       geom_bar(stat = "identity", color = "white",width=0.2) +
       geom_text_repel(aes(y=lab.ypos,label = prop), color = "white",size=fontSize)+
       scale_fill_manual(values = colores) +
       theme_void(base_size = 2*fontSize)+coord_polar('y',start=0) #+xlim(0.5, 2.5) 
  if (nargin>4){
    p <- p + ggtitle(titleStr) 
  }
  return(p)
}

## @knitr plotVennDiagram
plotVennDiagram <- function(elementsList,categories,colorValues,intLabSize,scaleF){
# Function for getting Venn diagrams for 2, 3 and 4 lists based on VennDiagram library
#
# Inputs: 
#   - elementsList  List of lists with the different subsets to plot
#   - categories    Vector of character vectors indicating the name label for each subset
#   - colorValues   Vector of character vectors indicating the colors for each subset
#   - intLabSize    Label sizes for each subset in the plot (unique and intersected regions)
#   - scaleF        TRUE if subset ellipses must be scaled according to their number of elements
# 
# Outputs:
#   - overlap     List with the overlap between all plotted subsets
#
  ellipses <- length(elementsList)
  nargin   <- length(as.list(match.call())) -1
  if (nargin < 5){scaleF = FALSE}
  euler <- TRUE
  overlap <- list()  
  if (ellipses  >2){
    area   <- c()
    shared <- c()
    for (i in 1:3){
      if (i==3){j<-1}
      else {j<-i+1}
      area[i]      <- length(elementsList[[i]])
      overlap[[i]] <- intersect(elementsList[[i]],elementsList[[j]])
      shared[i]    <- length(overlap[[i]])
    }
    #Intersect 1,2,3
    overlap[[4]] <- intersect(intersect(elementsList[[1]],elementsList[[2]]),elementsList[[3]])
    shared[4] <- length(overlap[[4]])
    
    if (ellipses <4){
      venn.plot<-draw.triple.venn(area[1], area[2],area[3], 
                                  shared[1],shared[2],shared[3],shared[4],
                                  category = categories,scaled = scaleF,euler.d = euler,fill = colorValues, 
                                  lty = c(rep("solid",3)),lwd = c(rep(0.5,3)), cex = intLabSize, cat.cex = 0.5*min(intLabSize),
                                  ext.text=FALSE)}
    else{
      area[4] <- length(elementsList[[4]])
      #Intersect 2,4
      shared[5] <- length(intersect(elementsList[[2]],elementsList[[4]]))
      #Intersect 3,4
      shared[6] <- length(intersect(elementsList[[3]],elementsList[[4]]))
      #Intersect 4,1
      shared[7] <- length(intersect(elementsList[[1]],elementsList[[4]]))
      #Intersect 1,2,4
      shared[8] <- length(intersect(intersect(elementsList[[1]],elementsList[[2]]),elementsList[[4]]))
      #Intersect 1,3,4
      shared[9] <- length(intersect(intersect(elementsList[[1]],elementsList[[3]]),elementsList[[4]]))
      #Intersect 2,3,4
      shared[10] <- length(intersect(intersect(elementsList[[2]],elementsList[[3]]),elementsList[[4]]))
      #Intersect 1,2,3,4
      shared[11] <- length(intersect(intersect(elementsList[[1]],elementsList[[2]]),intersect(elementsList[[3]],elementsList[[4]]))) 
      
      venn.plot <- draw.quad.venn(area[1], area[2],area[3],area[4], 
                                  shared[1],shared[3],shared[7],shared[2],shared[5],
                                  shared[6],shared[4],shared[8],shared[9],shared[10],shared[11],
                                  #n12, n13,n14, n23, n24,
                                  #n34, n123, n124, n134, n234, n1234, 
                                  category = categories,scaled = scaleF,fill = colorValues, 
                                  lty = c(rep("solid",4)), lwd = c(rep(0.5,4)),cex = intLabSize, cat.cex = 0.5*min(intLabSize),euler.d = euler,ext.text=FALSE)
      
    }
  }
  if (ellipses  == 2){
    area   <- c()
    shared <- c()
    for (i in 1:2){
      area[i]      <- length(elementsList[[i]])
      overlap[[i]] <- elementsList[[i]]
    } 
    overlap[[3]] <- intersect(elementsList[[1]],elementsList[[2]])
    #Intersect 1,2
    shared  <- length(overlap[[3]])
    venn.plot<-draw.pairwise.venn(area[1], area[2],shared,
                                  category = categories,scaled = scaleF,fill = colorValues, 
                                  lty = c(rep("solid",2)), lwd = c(rep(0.5,2)),cex = intLabSize, cat.cex = 0.5*min(intLabSize),euler.d = euler,ext.text=FALSE)
    
    
  }
  
  if (ellipses  == 1){
    area   <- c()
    shared <- c()
    area[1]      <- length(elementsList[[1]])
    overlap[[1]] <- elementsList[[1]]
    #Intersect 1,2
    venn.plot<-draw.single.venn(area[1], 
                                  category = categories,scaled = scaleF,fill = colorValues, 
                                  lty = c(rep("solid",1)), lwd = c(rep(0.5,1)),cex = intLabSize, cat.cex = 0.5*min(intLabSize),euler.d = euler,ext.text=FALSE)
    
  }
  output <- list(venn.plot,overlap)
  return(output)
}

## @knitr multipleDensityPlot
multipleDensityPlot <- function(df,xLabel,yLabel,fontSize,x_Limits,medianLine,labelStr,colors) {
  #Default parameters if 
  nargin <- length(as.list(match.call())) -1
  if (nargin<8){
    colors <- cividis(length(unique(df$class)))
    if (nargin<7){
      labelStr <- factor(unique(df$class))
      if (nargin<6){
        medianLine <- FALSE
      }
    }
  }
  medianLine <- FALSE
  # Density plots with semi-transparent fill
  p <- ggplot(df, aes(x=values, fill=class)) + geom_density(alpha=.3) + 
       scale_fill_manual(values = colors,labels=labelStr) +
       labs(y = yLabel, x=xLabel) + theme_classic(base_size = fontSize)+
       scale_x_log10(limits = c(1E-6,max(df$values)))
  if (medianLine){
    cdat <- ddply(df, "class", summarise, values.median=median(values))
    p    <- p + geom_vline(data = cdat,
               aes(xintercept = values.median), na.rm = T,linetype ="longdash", size = .3) 
  }
  plot(p)
}

## @knitr multipleDensityPlot2
multipleDensityPlot2 <- function(df,xLabel,yLabel,fontSize,x_Limits,medianLine,labelStr,colors) {
  #Default parameters if 
  nargin <- length(as.list(match.call())) -1
  if (nargin<8){
    colors <- cividis(length(unique(df$class)))#factor(unique(df$class))
    if (nargin<7){
      labelStr <- factor(unique(df$class))
      if (nargin<6){
        medianLine <- FALSE
      }
    }
  }
  medianLine <- FALSE
  # Density plots with semi-transparent fill
  p <- ggplot(df, aes(x=values, fill=class)) + geom_density(alpha=.3) + 
    scale_fill_manual(values = colors,labels=labelStr) +
    labs(y = yLabel, x=xLabel) + theme_classic(base_size = fontSize)+
    scale_x_log10(limits = c(1E-6,max(df$values)))
  if (medianLine){
    cdat <- ddply(df, "class", summarise, values.median=median(values))
    p    <- p + geom_vline(data = cdat,
                           aes(xintercept = values.median), na.rm = T,linetype ="longdash", size = .3) 
  }
  return(p)
}


## @knitr barPlot_counts
barPlot_counts <- function(dat1,xLabel,yLabel,fontSize,colors){
  nargin <- length(as.list(match.call())) -1
  if (nargin<5){
    colors <- factor(dat1$metGroup)
  }
  p <- ggplot(data=dat1, aes(x=Kingdom,y=counts,fill=metGroup)) + geom_bar(stat='identity')+ 
       scale_fill_manual(values = colors) + theme_bw(base_size = 2*fontSize)
  plot(p)
}  




