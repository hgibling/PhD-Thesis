library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(dplyr)


setwd("~/awadalla/PRDM9-Project/HMM/sim-reads")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

se.combined <- read.csv("se-counts-half/se-half-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="Single") %>%
    arrange(desc(Error))

pef.combined <- read.csv("pe-counts-half/pe-forward-half-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="PairedForward") %>%
    arrange(desc(Error))

per.combined <- read.csv("pe-counts-half2/pe-reverse-half-pr.csv", header=F, stringsAsFactors=F,
                         col.names=c("Precision", "Recall", "Coverage", 
                                     "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="PairedReverse") %>%
    arrange(desc(Error))

combined <- bind_rows(se.combined, pef.combined, per.combined)


se.each <- read.csv("se-counts-half/se-half-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="Single") %>%
    arrange(desc(Error))

pef.each <- read.csv("pe-counts-half/pe-forward-half-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="PairedForward") %>%
    arrange(desc(Error))

per.each <- read.csv("pe-counts-half2/pe-reverse-half-pr-each.csv", header=F, stringsAsFactors=F,
                     col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                 "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Reads="PairedReverse") %>%
    arrange(desc(Error))

each <- bind_rows(se.each, pef.each, per.each)

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
    ggtitle("50X coverage, 0% error rate")

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
    ggtitle("50X coverage, 0% error rate")

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
    ggtitle("50X coverage, 0% error rate")
