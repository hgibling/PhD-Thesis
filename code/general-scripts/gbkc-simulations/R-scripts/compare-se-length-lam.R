library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(dplyr)


setwd("~/awadalla/PRDM9-Project/HMM/sim-reads/count-tests/se-counts-length")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

combined <- read.csv("se-counts-lengths-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration", "ReadLength")) %>%
    group_by(Coverage, Error, kmer, ReadLength) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))

each <- read.csv("se-counts-lengths-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration", "ReadLength")) %>%
    group_by(Allele, Coverage, Error, kmer, ReadLength) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))


combined$ReadLength <- factor(combined$ReadLength, levels=seq(25, 500, 25))
each$ReadLength <- factor(each$ReadLength, levels=seq(25, 500, 25))
each$Allele <- factor(each$Allele, levels=allele.order)


# Plots

ggplot(combined %>% filter(ReadLength %in% seq(50, 500, 50)), aes(kmer, FScore)) +
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(each, aes(kmer, FScore)) +
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(each %>% filter(kmer>=75), aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=ReadLength)) +
    geom_line(aes(color=ReadLength)) +
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

# see where pe reads fall

pe.combined <- read.csv("~/awadalla/PRDM9-Project/HMM/sim-reads/pe-counts/pe-counts-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), ReadLength="Paired250") %>%
    arrange(desc(Error)) %>%
    dplyr::select(Coverage, Error, kmer, ReadLength, MeanPrecision, MeanRecall, FScore)

subset.combined <- combined %>%
    filter(as.numeric(ReadLength)==200 | as.numeric(ReadLength)==250) %>%
    mutate(ReadLength=paste0("Single", ReadLength)) 

se.pe <- bind_rows(pe.combined, subset.combined)

se.pe$ReadLength <- factor(se.pe$ReadLength, levels=c("Single200", "Single250", "Paired250"))


# Plots

ggplot(se.pe, aes(kmer, FScore)) +
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")
