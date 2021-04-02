#!/usr/bin/env Rscript
# Functions for loading and minor formatting of all different datasets used in this study
# Ivan Domenzain

## @knitr loadKcats
kcats <- read.csv(file = '../data/max_KCAT.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE)

## @knitr load_Sp_Act
Sp_Act <- read.csv(file = '../data/max_SA.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE)

## @knitr loadMWs
MWs <- read.csv(file = '../data/max_MW.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE)

## @knitr loadKMs
KMs <- read.csv(file = '../data/min_KM.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE)

## @knitr loadKEGG_enzData
enzData <- read.csv(file = '../data/enzyme.txt', sep = '\t', header = FALSE, stringsAsFactors = FALSE)