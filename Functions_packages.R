#Obtaining datasets and packages needed

library(devtools)
load_all("~/git/chmiddbb/")
dir <- chmi.phendir()
chmi.required_packages(
  locale = 'es_ES.UTF-8',
  update_packages = FALSE)

dat_file <- list.files(
  dir, pattern = '.csv',
  full.names = TRUE,
  recursive = TRUE)


### read data in 'dat.file'
dat_imp <- lapply(
  dat_file,
  fread,
  sep = ';',
  dec = '.',
  na.strings = c('', 'NA'))


# update 'dat_imp'
dat_imp <- chmi.update_phenos.dat_imp(dat_imp)

#Libraries pf packages needed

library(SummarizedExperiment)
library(DESeq2)
library(tweeDEseq)
library(stats)
library(pheatmap)
library(ggplot2)
library(ggfortify)
library(edgeR)
library(org.Hs.eg.db)
library(clusterProfiler)
library(sva)
library(EnsDb.Hsapiens.v79)
library(tidyverse)
library(mixOmics)
library(limma)
library(GSA)
library(VennDiagram)
library(FactoMineR)
library(annotate)


