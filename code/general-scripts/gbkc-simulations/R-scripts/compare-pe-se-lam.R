library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(dplyr)


setwd("~/awadalla/PRDM9-Project/HMM/sim-reads")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

se.combined <- read.csv("se-counts/se-counts-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="Single") %>%
    arrange(desc(Error))

se.250.combined <- read.csv("se-counts/se-counts-250bp-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="Single250") %>%
    arrange(desc(Error))

pe.combined <- read.csv("pe-counts/pe-counts-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="Paired") %>%
    arrange(desc(Error))

combined <- bind_rows(pe.combined, se.combined, se.250.combined)


se.each <- read.csv("se-counts/se-counts-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="Single") %>%
    arrange(desc(Error))

se.250.each <- read.csv("se-counts/se-counts-250bp-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="Single250") %>%
    arrange(desc(Error))

pe.each <- read.csv("pe-counts/pe-counts-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="Paired") %>%
    arrange(desc(Error))

each <- bind_rows(pe.each, se.each, se.250.each)

combined$Error <- factor(combined$Error, levels=c(0.005, 0.002, 0.001, 0))
each$Error <- factor(each$Error, levels=c(0.005, 0.002, 0.001, 0))
each$Allele <- factor(each$Allele, levels=allele.order)


# Plots

ggplot(combined, aes(kmer, FScore)) +
    theme_bw() +
    geom_point(aes(color=Reads)) +
    geom_line(aes(color=Reads)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(combined %>% filter(Reads!="Single250"), aes(kmer, FScore)) +
    theme_bw() +
    geom_point(aes(color=Reads)) +
    geom_line(aes(color=Reads)) +
    scale_color_manual(values=c("purple", "black")) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(each, aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=Reads)) +
    geom_line(aes(color=Reads)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(each %>% filter(kmer>=75), aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=Reads)) +
    geom_line(aes(color=Reads)) +
    scale_x_continuous(breaks=seq(75, 99, 4)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")
