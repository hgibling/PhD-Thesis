library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

setwd("~/prdm9")

diploid.all <- read.csv("counts-flank10k-pr.csv", header=F, stringsAsFactors=F,
                           col.names=c("kmer", "Genotype", "Precision", "Recall", "F1Score", "Coverage", "Error",
                                       "Iteration", "ReadLength", "FragmentLength"))

diploid.average <- diploid.all %>%
    filter(Genotype == "Average") %>%
    group_by(kmer, ReadLength, FragmentLength, Coverage, Error) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

diploid.each <- diploid.all %>%
    filter(Genotype != "Average") %>%
    group_by(Genotype, kmer, ReadLength, FragmentLength, Coverage, Error) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

ggplot(diploid.average, aes(kmer, MeanF1Score)) +
    facet_grid(Error ~ Coverage) +
    theme_bw() +
    geom_point() +
    geom_line() +
    theme(axis.text.x=element_text(size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("")

sub <- diploid.average
sub$Coverage <- factor(sub$Coverage, levels=c(20, 100))

ggplot(sub, aes(kmer, MeanF1Score)) +
    facet_grid(. ~ Error) +
    theme_bw() +
    geom_point(aes(color=Coverage)) +
    geom_line(aes(color=Coverage)) +
    theme(axis.text.x=element_text(size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("")

ggplot(diploid.each %>% filter(Coverage==100, Error==0.001), aes(Genotype, MeanF1Score)) +
    facet_wrap(~ kmer) +
    theme_bw() +
    geom_point() +
    geom_line() +
    theme(axis.text.x=element_text(size=15, angle=45, hjust=1),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("100X coverage, 0% error rate")


## compare 100x 0e flank10k and flank1k and python flank1k

library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

setwd("~/prdm9")


flank10k.all <- read.csv("counts-flank10k-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("kmer", "Genotype", "Precision", "Recall", "F1Score", "Coverage", "Error",
                                    "Iteration", "ReadLength", "FragmentLength")) %>%
    filter(Coverage==100, Error==0, Iteration<=25, kmer==51) %>%
    mutate(Type="gbkc-10k")

flank1k.all <- read.csv("counts-flank1k-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("kmer", "Genotype", "Precision", "Recall", "F1Score", "Coverage", "Error",
                                    "Iteration", "ReadLength", "FragmentLength")) %>%
    filter(Coverage==100, Error==0, Iteration<=25, kmer==51) %>%
    mutate(Type="gbkc-1k")

flank1k.py.all <- read.csv("counts-flank1k-python-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("kmer", "Genotype", "Precision", "Recall", "F1Score", "Coverage", "Error",
                                    "Iteration", "ReadLength", "FragmentLength")) %>%
    filter(Coverage==100, Error==0, Iteration<=25) %>%
    mutate(Type="Python-1k")

diploid.all <- flank10k.all %>%
    bind_rows(flank1k.all) %>%
    bind_rows(flank1k.py.all)

diploid.average <- diploid.all %>%
    filter(Genotype == "Average") %>%
    group_by(kmer, ReadLength, FragmentLength, Coverage, Error, Type) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

diploid.each <- diploid.all %>%
    filter(Genotype != "Average") %>%
    group_by(Genotype, kmer, ReadLength, FragmentLength, Coverage, Error, Type) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

ggplot(diploid.average, aes(kmer, MeanF1Score)) +
    facet_grid(Error ~ Coverage) +
    theme_bw() +
    geom_point(aes(color=Type)) +
    geom_line(aes(color=Type)) +
    theme(axis.text.x=element_text(size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("")

ggplot(diploid.each, aes(Genotype, MeanF1Score)) +
    facet_wrap(~ kmer) +
    theme_bw() +
    geom_point(aes(color=Flank), alpha=0.5) +
    #geom_line(aes(color=Flank)) +
    theme(axis.text.x=element_text(size=15, angle=45, hjust=1),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("")
