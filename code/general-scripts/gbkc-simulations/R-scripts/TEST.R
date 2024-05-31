library(tidyr)
library(ggplot2)
library(dplyr)

args <- commandArgs(trailingOnly=T)


setwd("/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/HMM/simulations")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))


haploid <- read.csv("distances/distance-flank10k-pr.csv", header=F, stringsAsFactors=F,
                     col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                 "Iteration", "ReadLength", "FragmentLength", "Method", "Flank")) %>%
    mutate(Method=paste0(Method, "-Haploid")) %>%
    select(-Flank) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    filter(kmer <= ReadLength)

haploid.each <- read.csv("distances/distance-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                 col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", 
                             "kmer", "Iteration", "ReadLength", "FragmentLength", "Method", "Flank")) %>%
    mutate(Method=paste0(Method, "-Haploid")) %>%
    select(-Flank) %>%
    group_by(ReadLength, FragmentLength, Allele, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

diploid <- read.csv("diploid-distances/distance-flank10k-pr.csv", header=F, stringsAsFactors=F,
                    col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                "Iteration", "ReadLength", "FragmentLength", "Method")) %>%
    mutate(Method=paste0(Method, "-Diploid")) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    filter(kmer <= ReadLength)

diploid.each <- read.csv("diploid-distances/distance-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                         col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", 
                                     "kmer", "Iteration", "ReadLength", "FragmentLength", "Method")) %>%
    mutate(Method=paste0(Method, "-Diploid")) %>%
    group_by(ReadLength, FragmentLength, Allele, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

counts <- read.csv("counts/counts-flank10k-pr.csv", header=F, stringsAsFactors=F,
                       col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                   "Iteration", "ReadLength", "FragmentLength", "Method", "Flank"))  %>% 
    mutate(Method=paste0(Method, "-Haploid")) %>%
    select(-Flank) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

counts.each <- read.csv("counts/counts-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                   col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", "kmer", 
                               "Iteration", "ReadLength", "FragmentLength", "Method", "Flank"))  %>% 
    mutate(Method=paste0(Method, "-Haploid")) %>%
    select(-Flank) %>%
    group_by(Allele, ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

distances.counts <- bind_rows(haploid, diploid, counts)
distances.counts$Error <- factor(distances.counts$Error, levels=unique(distances.counts$Error))
distances.counts$Coverage <- factor(distances.counts$Coverage, levels=unique(distances.counts$Coverage))
distances.counts$Method <- factor(distances.counts$Method, levels=unique(distances.counts$Method))


distances.counts.each <- bind_rows(haploid.each, diploid.each, counts.each)
distances.counts.each$Method <- factor(distances.counts.each$Method, levels=unique(distances.counts.each$Method))
distances.counts.each$Allele <- factor(distances.counts.each$Allele, levels=unique(distances.counts.each$Allele))

# Plots

cols <- c("red", "blue", "pink", "lightblue", "black")

distances.counts.each <- distances.counts.each %>%
    ungroup() %>%
    mutate(Allele=gsub("/.*", "", Allele))

x <- ggplot(distances.counts.each %>% filter(Allele=="A", Method=="max-Haploid", Error==0.01), aes(kmer, Coverage)) +
    theme_bw() +
    geom_contour(aes(z=FScore))

ggsave("test.pdf", x)
