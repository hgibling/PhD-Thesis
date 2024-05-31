#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
sample <- args[1]

library(tidyr)
library(dplyr)
library(ggplot2)

base.path <- "/.mounts/labs/simpsonlab/users/hgibling/OHS/illumina"
count.path <- paste(base.path, "gbkc-count-calls-corrected-36", sep="/")
distance.path <- paste(base.path, "gbkc-distance-calls-corrected-36", sep="/")

# Files
count.data <- read.table(paste(count.path, 
                               paste0(sample, ".filtered-paired-corrected-all-methods.tsv"), sep="/"), 
                               col.names=c("CountMethod", "k", 
                               "Genotype", "CountScore"))
distance.data <- read.table(paste(distance.path, 
                            paste0(sample, ".filtered-paired-corrected-all-methods.tsv"), sep="/"),
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

write.table(combined, paste0(base.path, "/combined-calls-corrected-36/", sample,
                             "-filtered-combined-calls-corrected-top100.tsv"),
            sep="\t", quote=F, col.names=T, row.names=F)
