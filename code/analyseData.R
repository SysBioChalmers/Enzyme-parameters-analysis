#!/usr/bin/env Rscript
# Functions for processing different datasets used in this study
# Ivan Domenzain. Last modified: 2019-07-17

## @knitr kcat_dispersion
#Function that gets the dispersion of Kcat value per enzyme class (EC#)
kcat_dispersion <- function(dataset){
  ECnumbers <- unique(dataset[,1])
  spanning  <- c()
  spreading <- c()
  N_entries <- c()
  EC_num    <- c()
  #Loop through each EC number and get the dispersion for their activity values
  for(ECnum in ECnumbers){
    indexes  <- which(dataset[,1]==ECnum)
    entries  <- length(indexes)
    if (entries >0){
      EC_num    <- c(EC_num,ECnum)
      N_entries <- c(N_entries,entries)
      #get distribution of values for the i-th EC#
      dist <- as.numeric(dataset[indexes,4])
      #Get the order of magnitude of the maximum and minimum values of dist
      maxOrder <- floor(log10(max(dist)))
      minOrder <- floor(log10(min(dist)))
      #Spanning is calculated as the amount of orders of magnitude that a given 
      #Kcat values distribution spans 
      range    <- maxOrder - minOrder
      spanning <- c(spanning,range)
      if (length(indexes)>2){
        #If there are more than 2 values reported for the i-th EC#
        #then a metric of "spreading" can be calculated as the ratio 
        #between the order of magnitude of the distribution's characteristic
        #value (median) and the spanning of the distribution
        maxOrder  <- log10(max(dist))
        minOrder  <- log10(min(dist))
        medOrder  <- log10(median(dist))
        range     <- abs(medOrder)/abs(maxOrder - minOrder)
        spreading <- c(spreading,range)
      }
    }
  }
  sorted      <- sort(N_entries,decreasing=TRUE,index.return=TRUE)
  percentages <- as.character(paste(round(as.numeric(sorted$x)*100/nrow(dataset),digits=2),' %'))
  N_entries   <- as.data.frame(cbind(EC_num[sorted$ix],as.numeric(sorted$x),percentages),stringsAsFactors = FALSE )
  colnames(N_entries) <- c('ECnumbers','Entries','% of total entries')
  
  plotVariables <- (list(N_entries,spanning,spreading))
  return(plotVariables)
} 

## @knitr data_orgs_composition
#Function that analyses the Kcats dataset composition in terms of studied
#organisms
data_orgs_composition <- function(dataset){
  organisms   <- unique(dataset[,3])
  org_entries <- c()
  org_class   <- c()
  #
  for (org in organisms){
    #Get the total number of entries reported for the i-th organism
    org_entries <- c(org_entries,length(which(dataset[,3]==org)))
    #Get organism class
    positions   <- gregexpr(pattern =';',org,fixed=TRUE)
    positions   <- positions[[1]]
    if (length(positions)>1){
      class             <- substr(org,positions[1]+1,positions[2]-1)
      substr(class,1,1) <- toupper(substr(class,1,1))
    } else {class <- 'non-KEGG'}
    org_class <- c(org_class,class)
  }
  #Sort of organisms by their number of entries in the dataset
  sorted      <- sort(org_entries,decreasing=TRUE,index.return=TRUE)
  percentages <- as.character(paste(round(as.numeric(sorted$x)*100/nrow(dataset),digits=2),' %'))
  appearences <- as.data.frame(cbind(organisms[sorted$ix],org_class[sorted$ix],as.numeric(sorted$x),percentages),stringsAsFactors = FALSE )
  colnames(appearences) <- c('Organism','Kingdom','Entries','% of total entries')
  #Capitalize first letter of organism names
  substr(appearences$Organism, 1, 1) <- toupper(substr(appearences$Organism, 1, 1))
  appearences$Entries                <- as.numeric(appearences$Entries)
  #dataset            <- cbind(dataset,org_class,stringsAsFactors = FALSE)
  #Remove tax. classification for appearences dataframe
  appearences[,1]    <- gsub("//.+$",'',appearences[,1],fixed=FALSE)
  return(list(appearences))
} 

## @knitr statTest_distributions
#Function that takes a list, with a distribution of values in each of its cells 
#and then computes a Kolmogorov-Smirnov pairwise statistical test to all the different
#combination of pairs of distributions
statTest_distributions <- function(distributions,tags){
  pVal_table <- c()
  for (dist_i in distributions){
    vector <- c()
    for (dist_j in distributions){
      pVal   <-ks.test(dist_i,dist_j)
      pVal   <- pVal$p.value
      vector <- c(vector,pVal)
    }
    pVal_table <- cbind(pVal_table,vector)
  } 
  df           <- data.frame(pVal_table,row.names = tags)
  colnames(df) <- tags
  return(df)
}

## @knitr add_WildCards
#Function that takes a dataset of enzyme parameters and adds the number of requested 
#wild cards to all of the EC numbers present in it. 
add_WildCards <- function(dataset,numWC){
  position <- 4-numWC
  for (i in 1:nrow(dataset)){
    ECnum    <- dataset[i,1]
    dots_Pos <- gregexpr(pattern='.',ECnum,fixed = TRUE)
    dots_Pos <- dots_Pos[[1]]
    if (length(dots_Pos==3)){
      if (position>=1){
        str          <- gsub(", ","",toString(rep('X.',numWC)))
        newECnum     <- paste(substr(ECnum,1,dots_Pos[position]),substr(str,1,(nchar(str)-1)),sep='')
        dataset[i,1] <- newECnum
      }
    }
  }
  return(dataset)
}

## @knitr get_WC1_distributions
#Function that takes a dataset of kinetic parameters where EC numbers have been modified 
#with one wild-card (df_WC1) and another one with two (df_WC2). The function gets the top represented 
#enzyme family (with two wild-cards) in df_WC2 and then gets a subset dataset for this enzyme family
#from df_WC1
get_WC1_distributions <- function(df_WC1,df_WC2,Ndists,rank){
  #Get the most represented families in the WC2 dataset
  df        <- as.data.frame(sort(table(df_WC2[,1]),decreasing=TRUE)[1:10])
  topEC_WC2 <- df$Var1[rank]
  #Get  a subset dataset from the the WC1 that matches with the top EC_WC2
  pos        <- gregexpr('.',topEC_WC2,fixed=TRUE)
  str        <- substr(topEC_WC2,1,pos[[1]][2])
  subset_WC1 <- df_WC1[grep(str,df_WC1[,1],fixed=TRUE),]
  classesWC1 <- unique(subset_WC1[,1])
  classesWC1 <- classesWC1[sample((1:length(classesWC1)),Ndists)]
  #Get distributions of values and strings for plotting
  colnames(df) <- c('Enzyme groups','Number of entries')
  output       <- list(subset_WC1,classesWC1,df)
  return(output)
}

## @knitr getDistributionsList
#Function that takes a dataset and list of classes for generating a list of distributions of values for each 
#unique class. The class comparison is done in the column provided as input.
getDistributionsList <- function(dataset,classes,column){
  average  <- c()
  V        <- list()
  labels   <- c()
  uniqueEC <- list()
  for (i in 1:length(classes)){
    class          <- classes[i]
    indexes        <- grep(class,dataset[,column])
    #indexes        <- gregexpr(class,dataset[,column],fixed=TRUE)
    #indexes       <- indexes[[1]]
    uniqueEC[[i]]  <- unique(dataset[indexes,1]) 
    #Get the average number of entries per unique classes
    average        <- c(average,floor(length(indexes)/length(uniqueEC[[i]])))
    V[[i]]         <- as.numeric(dataset[indexes,4])
    str            <- paste(class,'/',length(V[[i]]),'/',round(median(V[[i]]), digits = 2))
    labels         <- c(labels,str)
  }
  output <- list(V,labels,average,uniqueEC)
  return(output)
}

## @knitr DistList_to_DataFrame
#Function that gets a list with distributions of values and transforms it into a dataframe,
#adding additional information that allow the classification of such values.
DistList_to_DataFrame <- function(distList,classes,colors,tags){
  nargin <- length(as.list(match.call())) -1
  tag    <- c()
  class  <- c()
  values <- c()
  colours <- c()
  for (i in 1:length(distList)){
    nL     <- length(distList[[i]])
    values <- c(values,as.numeric(distList[[i]]))
    class  <- c(class,rep(classes[i],nL))
    if (nargin>2){
      colours <- c(colours,rep(colors[i],nL))
      if (nargin>3){
        tag <- c(tag,rep(tags[i],nL))
     }
    }
  }
  
  if (nargin>2) {
    df <- data.frame(class,values,colours,stringsAsFactors = TRUE)
    if (nargin>3) {
      df <- data.frame(tag,class,values,colours,stringsAsFactors = TRUE)
    }
  }else {
    df <- data.frame(class,values,stringsAsFactors = TRUE)
  }
}

## @knitr extract_data_subset
#Function that takes a dataset and extracts a subset of it 
extract_data_subset <- function(dataset,category,column,getUnique_EC){
  nargin <- length(as.list(match.call())) -1
  if (nargin<4){getUnique_EC <- FALSE}
  indexes  <- grep(category,dataset[,column],fixed=TRUE)
  subsetDF <- dataset[indexes,] 
  if (getUnique_EC){
    uniqueEC <- unique(subsetDF[,1])
    output   <- list(subsetDF,uniqueEC)
  } else{output <- subsetDF}
  return(output)
}

## @knitr get_orgs_metGroups_table
get_orgs_metGroups_table <- function(df,kingdoms,metGroups){
  countsTable <- c()
  counts   <- c()
  Kingdom  <- c()
  metGroup <- c()
  for (class in kingdoms){
    countVect <- c()
    indexesClass  <- grep(class,df[,3],fixed=TRUE)
    for (group in metGroups){
      indexesGroup  <- grep(group,df[,7],fixed=TRUE)
      indexes       <- intersect(indexesClass,indexesGroup)
      counts      <- c(counts,length(indexes))
      Kingdom     <- c(Kingdom,class)
      metGroup    <- c(metGroup,group)
    }
  }
  countsTable <- data.frame(counts,Kingdom,metGroup,stringsAsFactors = FALSE)
  return(countsTable)
}
