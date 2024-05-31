library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(dplyr)


setwd("~/awadalla/PRDM9-Project/HMM/sim-reads/pe-frag/")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

combined <- read.csv("pe-counts-frag-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration", "FragmentLength")) %>%
    group_by(Coverage, Error, kmer, FragmentLength) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))

each <- read.csv("pe-counts-frag-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration", "FragmentLength")) %>%
    group_by(Allele, Coverage, Error, kmer, FragmentLength) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))


combined$FragmentLength <- factor(combined$FragmentLength, levels=seq(200, 500, 50))
each$FragmentLength <- factor(each$FragmentLength, levels=seq(200, 500, 50))
each$Allele <- factor(each$Allele, levels=allele.order)


# Plots

ggplot(combined, aes(kmer, FScore)) +
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(each, aes(kmer, FScore)) +
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(each %>% filter(kmer>=75), aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=FragmentLength)) +
    geom_line(aes(color=FragmentLength)) +
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

# check kmer distributions

k1 <- read.csv("kmer-count-check.csv", header=F, stringsAsFactors=F,
              col.names=c("Allele", "kmer", "Count", "FragmentLength")) %>%
    mutate(Sample=1)
k2 <- read.csv("kmer-count-check2.csv", header=F, stringsAsFactors=F,
               col.names=c("Allele", "kmer", "Count", "FragmentLength")) %>%
    mutate(Sample=2)
k3 <- read.csv("kmer-count-check3.csv", header=F, stringsAsFactors=F,
               col.names=c("Allele", "kmer", "Count", "FragmentLength")) %>%
    mutate(Sample=3)
k4 <- read.csv("kmer-count-check4.csv", header=F, stringsAsFactors=F,
               col.names=c("Allele", "kmer", "Count", "FragmentLength")) %>%
    mutate(Sample=4)

k <- bind_rows(k1, k2, k3, k4)

k$FragmentLength <- factor(k$FragmentLength, levels=seq(200, 500, 50))


ggplot(k %>% filter(FragmentLength==200 | FragmentLength==500), aes(Count)) +
    theme_bw() +
    facet_wrap(~ Sample) +
    geom_histogram(aes(fill=FragmentLength), position="identity", alpha=0.3,
                   binwidth=1) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    xlab("Number of kmers") +
    ylab("kmer Count") +
    ggtitle("100X coverage, 0% error rate")


kf <- read.csv("check-first-mate-frag.csv", header=F, stringsAsFactors=F,
               col.names=c("Allele", "kmer", "Count", "FragmentLength", "Iteration")) %>%
    mutate(Iteration=case_when(
        Iteration == 50 ~ 2,
        Iteration == 100 ~ 3,
        Iteration == 13 ~ 4,
        T ~ 1))
kr <- read.csv("check-second-mate-frag.csv", header=F, stringsAsFactors=F,
               col.names=c("Allele", "kmer", "Count", "FragmentLength", "Iteration")) %>%
    mutate(Iteration=case_when(
        Iteration == 50 ~ 2,
        Iteration == 100 ~ 3,
        Iteration == 13 ~ 4,
        T ~ 1))

kf$FragmentLength <- factor(kf$FragmentLength, levels=seq(200, 500, 50))
kr$FragmentLength <- factor(kr$FragmentLength, levels=seq(200, 500, 50))


ggplot(kf, aes(Count)) +
    theme_bw() +
    geom_histogram(aes(fill=FragmentLength), position="identity", alpha=0.3,
                   binwidth=1) +
    facet_wrap(~ Iteration) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    xlab("Number of kmers") +
    ylab("kmer Count") +
    ggtitle("50X coverage, first mates, 0% error rate")

ggplot(kr, aes(Count)) +
    theme_bw() +
    geom_histogram(aes(fill=FragmentLength), position="identity", alpha=0.3,
                   binwidth=1) +
    facet_wrap(~ Iteration) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    xlab("Number of kmers") +
    ylab("kmer Count") +
    ggtitle("50X coverage, second mates, 0% error rate")


# merge?

kf <- kf %>%
    rename(FirstCount=Count)
kr <- kr %>%
    rename(SecondCount=Count)

kfr <- full_join(kf, kr, by=c("Allele", "kmer", "FragmentLength", "Iteration")) %>%
    mutate(FirstCount=ifelse(is.na(FirstCount), 0, FirstCount),
           SecondCount=ifelse(is.na(SecondCount), 0, SecondCount)) %>%
    mutate(TotalCount=FirstCount+SecondCount)

kfr$FragmentLength <- factor(kfr$FragmentLength, levels=seq(200, 500, 50))

ggplot(kfr, aes(TotalCount)) +
    theme_bw() +
    geom_histogram(aes(fill=FragmentLength), position="identity", alpha=0.3,
                   binwidth=1) +
    facet_wrap(~ Iteration) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    xlab("Number of kmers") +
    ylab("kmer Count") +
    ggtitle("50X coverage, second mates, 0% error rate")

# same as pe