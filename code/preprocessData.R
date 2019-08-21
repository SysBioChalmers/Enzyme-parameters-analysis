#!/usr/bin/env Rscript
# Functions for processing different datasets used in this study
# Ivan Domenzain

#setwd('/Users/ivand/Documents/GitHub/Kcat_analysis/code')

## @knitr convertSAdata
# Kcats are also equivalent to SA*MW so more data entries can be 
#obtained from the loaded BRENDA datasets
convertSAdata <- function(Sp_Act,MWs) {
  numel <- nrow(Sp_Act)
  #Loop through each SA entry and find if there exists a reported 
  #MW for the same EC# and organism
  SA_kcats   <- data.frame()
  for(i in 1:numel) {
    ECnum    <- Sp_Act[i,1]
    organism <- Sp_Act[i,3]
    indexes  <- which(MWs[,1]==ECnum & MWs[,3]==organism)
    if (length(indexes)>0){
      MWeigth  <- max(MWs[indexes,4])
      newKcat  <- Sp_Act[i,4]*MWeigth/(60*1000)
      newRow   <- cbind(ECnum,'*',organism,newKcat,'*')
      SA_kcats <- rbind(SA_kcats,newRow)
    }
  }
  SA_kcats           <- data.frame(SA_kcats,stringsAsFactors = FALSE)
  colnames(SA_kcats) <- colnames(kcats)
  return(SA_kcats)
}

## @knitr extendKcatData
# Complement kcat data with those couples EC#-organism that are 
# present in SA_kcats but not in the original dataset
extendKcatData <- function(kcats,SA_kcats) {
  ext_Kcats <- data.frame()
  #Extract unique EC numbers
  SA_EC      <- as.character(unique(SA_kcats[,1]))
  ECnumbers  <- unique(c(unique(kcats[,1]),SA_EC))

  #Loop through each unique enzyme class
  for(i in 1:length(ECnumbers)) {
    ECnum      <- ECnumbers[i]
    #First add all the entries from kcat file (with substrate info)
    indexes <- which(kcats[,1]==ECnum)
    entries <- length(indexes)
    if (entries>0){
      ext_Kcats <- rbind(ext_Kcats,kcats[indexes,])
      orgsKcats <- kcats[indexes,3]
      #Search for additional entries in SA_kcats for the same EC#
      indexes   <- which(SA_kcats[,1]==ECnum)
      orgsSAct  <- SA_kcats[indexes,3]
      indexes   <- indexes[is.na(match(orgsSAct,orgsKcats))]
      ext_Kcats <- rbind(ext_Kcats,SA_kcats[indexes,])
     } else{
      #search ec number in SA_kcats
      indexes   <- which(SA_kcats[,1]==ECnum)
      ext_Kcats <- rbind(ext_Kcats,SA_kcats[indexes,],stringsAsFactors = FALSE)
     }
  }
  #save enzymes class information
  enzClasses <- c('EC1.X.X.X','EC2.X.X.X','EC3.X.X.X','EC4.X.X.X','EC5.X.X.X','EC6.X.X.X')
  ECclass    <- rep('',nrow(ext_Kcats))
  for (class in enzClasses){
    ECgroup <- substr(as.character(class),1,4)
    indexes <- grep(ECgroup,ext_Kcats[,1],value=FALSE,fixed=TRUE)
    ECclass[indexes] <- class
  }
  ext_Kcats <- cbind(ext_Kcats,ECclass,stringsAsFactors = FALSE)
  return(ext_Kcats)
}

## @knitr processKEGG_enzymes
#Get related pathways and genes for each EC number available in the KEGG ftp database
processKEGG_enzymes <- function(enzData) {
  delimiters <- grep('ENTRY       EC',enzData[,1])
  EC_info    <- c()
  for (i in 1:(length(delimiters)-1)){
    #Get data for the i-th enzyme
    ECdata           <- c()
    ECdata           <- as.data.frame(enzData[((delimiters[i]):(delimiters[i+1]-1)),1],stringsAsFactors = FALSE)
    colnames(ECdata) <- 'column'
    #Get EC number
    ECpos  <- gregexpr(pattern='EC ',ECdata$column[1],fixed = TRUE)
    ECpos  <- ECpos[[1]]
    if (length(ECpos)==1){
      ECpos <- ECpos[1]
      dataStr <- as.character(ECdata$column[1])
      ECnum   <- as.character(substr(dataStr,ECpos,nchar(dataStr)))
      blanks  <- gregexpr(pattern=' ',ECnum,fixed = TRUE)
      obsEnz  <- gregexpr(pattern='Obsolete  Enzyme',ECnum,fixed = TRUE)
      #discard empty and obsolote EC numbers
      if (nchar(ECnum)>3 & sum(obsEnz[[1]])<1){
        ECnum     <- substr(ECnum,1,blanks[[1]][2]-1)
        subst_pos <- as.numeric(grep('SUBSTRATE  ',ECdata$column,fixed = TRUE))
        produ_pos <- as.numeric(grep('PRODUCT    ',ECdata$column,fixed = TRUE))
        #Add default values
        substrates <- '//'
        pathways   <- '//'
        #Get substrate information
        if (sum(subst_pos)>0 & sum(produ_pos)>0){
          subData    <- as.character(ECdata$column[((subst_pos):(produ_pos-1))])
          substrates <- c()
          for (rowStr in subData){
            #Extract the substrate name from each substrate row
            bracket    <- gregexpr(pattern='[CPD:',rowStr,fixed = TRUE)
            rowStr     <- substr(rowStr,1,bracket[[1]]-2)
            rowStr     <- gsub('SUBSTRATE','',rowStr)
            rowStr     <- trimws(rowStr,which="left")
            substrates <- tolower(c(substrates,rowStr))
          }
          substrates <- paste(substrates, collapse="//") 
        }
        #Get pathway information
        path_pos <- grep('PATHWAY    ',ECdata$column,fixed = TRUE)
        if (length(path_pos)>=1){
          pathways <- c()
          orth_pos <- grep('ORTHOLOGY  ',ECdata$column,fixed = TRUE)
          gene_pos <- grep('GENES      ',ECdata$column,fixed = TRUE)
          link_pos <- grep('DBLINKS    ',ECdata$column,fixed = TRUE)
          hist_pos <- grep('HISTORY    ',ECdata$column,fixed = TRUE)
          comm_pos <- grep('COMMENT    ',ECdata$column,fixed = TRUE)
          jour_pos <- grep('JOURNAL    ',ECdata$column,fixed = TRUE)
          #Pathway information lie between PATHWAY and the next category
          indxVect <- c(orth_pos,gene_pos,link_pos,hist_pos,comm_pos)
          indxVect <- indxVect[which(indxVect>path_pos)]
          index    <- min(indxVect)
          if (length(index)>=1){
            subData  <- ECdata$column[((path_pos):(index-1))]
            subData  <- gsub('PATHWAY','',subData)
            subData  <- trimws(subData,which="left")
            pathways <- paste(subData, collapse="//")
          }
        }
        EC_info <- rbind(EC_info,cbind(ECnum,substrates,pathways))
      }
    }
    #Update the position in the enzyme file
    rowIndx  <- delimiters[i]
  }
  EC_info[,1] <- gsub(' ','',EC_info[,1])
  EC_info[,3] <- gsub('  ',' ',EC_info[,3])
  indexes     <- grep('EC7.',EC_info[,1],fixed = TRUE,invert=TRUE)
  EC_info     <- EC_info[indexes,]
  EC_info     <- data.frame(EC_info,stringsAsFactors = FALSE)
  return(EC_info)
} 


## @knitr classify_pathways
#Get related pathways and genes for each EC number available in the KEGG ftp database
classify_pathways <- function(EC_info) {
  #Get KEGG pathways classification file
  pathways   <- read.csv(file = '../data/pathway.txt', sep = '\n', header = FALSE, stringsAsFactors = FALSE)
  delimiters <- grep('##',pathways[,1])
  groups     <- gsub('##','',pathways[delimiters,1],fixed = TRUE)
  grouping   <- list(c('Carbohydrate metabolism','Energy metabolism'),
                     c('Lipid metabolism','Nucleotide metabolism','Amino acid metabolism'),
                     c('Metabolism of other amino acids','Glycan biosynthesis and metabolism','Metabolism of cofactors and vitamins'))
  groupNames <- c('CEM','ALM','ISM')
  #Get pathGroups: a list with all related pathways for the indicated groups above
  pathGroups <- data.frame()
  i <- 1
  for (group in grouping){
    vector <- c()
    for (pathway in group){
      index    <- grep(pathway,groups,fixed=TRUE)
      nextIndx <- delimiters[index+1]
      index    <- delimiters[index]
      vector   <- c(vector,pathways[(index+1):(nextIndx-1),1])
    }
    vector     <- gsub('\t',' ',vector)
    newRow     <- cbind(rep('ec',length(vector)),vector,rep(groupNames[i],length(vector)))
    pathGroups <- rbind(pathGroups,newRow)
    i <- i+1
  }
  EC_info[,'metGroup'] <- NA
  pathGroups           <- unite(pathGroups,'pathways','V1','vector', sep = "", remove = TRUE)
  colnames(pathGroups) <- c('pathways','group')
  pathGroups[,1]       <- trimws(pathGroups[,1],which="both")
  pathGroups           <- data.frame(pathGroups,stringsAsFactors = FALSE)
  #Match metabolic groups to pathway information for each EC number entry
  for (i in 1:nrow(EC_info)){
    paths  <- strsplit(EC_info$pathways[i], '//', fixed = TRUE)
    groups <- c()
    for (path in paths[[1]]){
      #Search for the i-th pathway in the pathGroups dataframe
      if (nchar(path)>0){
        index <- grep(path,pathGroups$pathways,fixed = TRUE)
        if (sum(index)>=1){
          groups <- c(groups,as.character(pathGroups[index,2]))
        }
      }
    }
    if(length(groups)>=1){groups <- paste(unique(groups), collapse="//") }
    else {groups <- 'Others'}
    #Add metGroup information to the EC_info kegg dataframe
    EC_info$metGroup[i] <- groups
  }
  return(list(EC_info,pathGroups))
}

## @knitr add_MetGroup
#Add metabolic group information to each enzyme entry in the studied dataset
add_MetGroup <- function(dataset,EC_info){
  dataset[,'metGroup'] <- NA
  for (i in 1:nrow(dataset)){
    ecNum <- dataset[i,1]
    #Search EC number entry in KEGG EC dataset
    index <- grep(ecNum,EC_info$ECnum,fixed = TRUE)
    if (sum(index)>=1){
      metGroup            <- as.character(EC_info$metGroup[index])
      dataset$metGroup[i] <- metGroup
    }
  }
  dataset <- data.frame(dataset,stringsAsFactors = FALSE)
  colnames(dataset) <- c('ECnumber','Substrate','Organism','Kcat','StdDev','enz_Family','metGroups')
  return(dataset)
}

