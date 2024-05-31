library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(dplyr)


setwd("~/awadalla/sim-reads/flank10k")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))


combined <- read.csv("flank10-pr.csv", header=F, stringsAsFactors=F,
                     col.names=c("Type", "Precision", "Recall", "Coverage", 
                                 "Error", "kmer", "Iteration")) %>%
    filter(Coverage==100, Error==0) %>%
    group_by(Type, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))

each <- read.csv("flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                 col.names=c("Type", "Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                             "Error", "kmer", "Iteration")) %>%
    filter(Coverage==100, Error==0) %>%
    group_by(Type, Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))


combined$Type <- factor(combined$Type, levels=unique(combined$Type)
each$Type <- factor(each$Type, levels=unique(each$Type)
each$Allele <- factor(each$Allele, levels=allele.order)


# Plots

ggplot(combined, aes(kmer, FScore)) +
    theme_bw() +
    geom_point(aes(color=Type)) +
    geom_line(aes(color=Type)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("100X coverage, 0% error rate, 10k flank")


ggplot(each, aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=Type)) +
    geom_line(aes(color=Type)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F1 Score") +
    ggtitle("100X coverage, 0% error rate, 10k flank")

ggplot(each %>% filter(kmer>=75), aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=Type)) +
    geom_line(aes(color=Type)) +
    scale_x_continuous(breaks=seq(75, 99, 4)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F1 Score") +
    ggtitle("100X coverage, 0% error rate, 10k flank")


