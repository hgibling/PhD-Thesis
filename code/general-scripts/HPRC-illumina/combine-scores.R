#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]

library(tidyr)
library(dplyr)
library(ggplot2)

base.path <- "/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/validation"
illumina.path <- paste(base.path, "PanGenomeProject/illumina-genotyping", sep="/")
count.path <- paste(illumina.path, "gbkc-calls", sep="/")
distance.path <- paste(illumina.path, "gbkc-distance-calls", sep="/")

# Files
count.data <- read.table(paste(count.path, 
                               paste0(sample, ".filtered-all-methods.tsv"), sep="/"), 
                               col.names=c("CountMethod", "k", 
                               "Genotype", "CountScore"))
distance.data <- read.table(paste(distance.path, 
                            paste0(sample, ".filtered-paired-all-methods.tsv"), sep="/"),
                            col.names=c("DistanceMethod", "k", 
                            "Genotype", "DistanceScore"))

# Data prep
counts <- count.data %>%
  group_by(CountMethod, k) %>%
  mutate(CountRank=rank(-CountScore, ties.method="min")) %>% 
  ungroup()

distances <- distance.data %>%
  group_by(DistanceMethod, k) %>%
  mutate(DistanceRank=rank(-DistanceScore, ties.method="min")) %>% 
  ungroup()

combined <- counts %>% 
  full_join(distances) %>% 
  mutate(CombinedMethod=paste(CountMethod, DistanceMethod, sep="_"),
         CombinedScore=CountScore+DistanceScore) %>%
  select(Genotype, k, CombinedMethod, CombinedScore) %>%
  group_by(k, CombinedMethod) %>% 
  mutate(CombinedRank=rank(-CombinedScore, ties.method="min")) %>% 
  filter(CombinedRank<=100) %>%
  arrange(k, CombinedMethod, CombinedRank)
  
write.table(combined, paste0(illumina.path, "/combined-calls/", sample,
                             "-filtered-combined-calls-top100.tsv"),
            sep="\t", quote=F, col.names=T, row.names=F)
