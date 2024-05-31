library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)


setwd("~/awadalla/test-sims-fragments")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))


combined <- read.csv("frag-pr.csv", header=F, stringsAsFactors=F,
                     col.names=c("Precision", "Recall", "Coverage", "Error", 
                                 "kmer", "Iteration", "ReadLength", "FragmentLength")) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    filter(kmer <= ReadLength)

each <- read.csv("frag-pr-each.csv", header=F, stringsAsFactors=F,
                 col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                             "Error", "kmer", "Iteration", "ReadLength", "FragmentLength")) %>%
    group_by(ReadLength, FragmentLength, Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))


combined$ReadLength <- factor(combined$ReadLength, levels=unique(combined$ReadLength))
combined$FragmentLength <- factor(combined$FragmentLength, levels=unique(combined$FragmentLength))

each$Allele <- factor(each$Allele, levels=allele.order)
each$ReadLength <- factor(each$ReadLength, levels=unique(each$ReadLength))
each$FragmentLength <- factor(each$FragmentLength, levels=unique(each$FragmentLength))


# Plots

ggplot(combined %>% filter(ReadLength==100), aes(kmer, FScore)) +
    theme_bw() +
    geom_point(aes(color=FragmentLength)) +
    geom_line(aes(color=FragmentLength)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("100X coverage, 0% error rate, 10k flank")

ggplot(combined %>% filter(FragmentLength==250, kmer>13), aes(kmer, FScore)) +
    theme_bw() +
    geom_point(aes(color=ReadLength)) +
    geom_line(aes(color=ReadLength)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("100X coverage, 0% error rate, 10k flank")


ggplot(each %>% filter(ReadLength==100, kmer > 13), aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=FragmentLength)) +
    geom_line(aes(color=FragmentLength)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F1 Score") +
    ggtitle("100X coverage, 0% error rate, 10k flank")

ggplot(each %>% filter(FragmentLength==250, kmer > 13), aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=ReadLength)) +
    geom_line(aes(color=ReadLength)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F1 Score") +
    ggtitle("100X coverage, 0% error rate, 10k flank")


