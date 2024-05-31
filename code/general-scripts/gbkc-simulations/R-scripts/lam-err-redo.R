library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)


setwd("~/awadalla/PRDM9-Project/HMM/sim-reads/counts/pr-meta")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

poisson <- read.csv("lam-err-redo-meta-batch.csv", header=F,
                    stringsAsFactors=F,
                    col.names=c("Count", "Coverage", "Error", "kmer", "Correct", "Incorrect", "Tied"))

poisson.pr.base <- read.csv("counts-pr.csv", header=F, stringsAsFactors=F,
                       col.names=c("Precision", "Recall", "Coverage", 
                                   "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))

poisson.pr <- poisson.pr.base
poisson.pr$Error <- factor(poisson.pr$Error, levels=c(0.005, 0.002, 0.001, 0))

greens <- brewer.pal(9, "Greens")[c(9, 7, 5, 4)]

ggplot(poisson.pr, aes(Coverage, FScore)) +
    theme_bw() +
    facet_wrap(~ kmer) +
    scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
                       name="Substitution\nerror",
                       labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    geom_line(aes(color=Error)) +
    geom_point(aes(color=Error)) +
    scale_x_continuous(breaks=seq(40, 200, 40)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F Score")

ggplot(poisson.pr, aes(kmer, FScore)) +
    theme_bw() +
    facet_wrap(~ Error) +
    #scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
    #                   name="Substitution\nerror",
    #                   labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    #geom_line(aes(color=kmer)) +
    geom_point(aes(color=kmer)) +
    #scale_x_continuous(breaks=seq(40, 200, 40)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F Score")



ggplot(poisson.pr %>% filter(kmer %in% c(1:7, 47:53, 93:99)), aes(Coverage, FScore)) +
    theme_bw() +
    facet_wrap(~ kmer) +
    scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
                       name="Substitution\nerror",
                       labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    geom_line(aes(color=Error)) +
    geom_point(aes(color=Error)) +
    scale_x_continuous(breaks=seq(20, 200, 20)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F Score")


ggplot(poisson.pr %>% filter(kmer %in% c(43:59)), aes(Coverage, FScore)) +
    theme_bw() +
    facet_wrap(~ kmer) +
    scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
                       name="Substitution\nerror",
                       labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    geom_line(aes(color=Error)) +
    geom_point(aes(color=Error)) +
    scale_x_continuous(breaks=seq(20, 200, 20)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F Score")



ggplot(poisson.pr, aes(kmer, FScore)) +
    theme_bw() +
    facet_wrap(~ Coverage, nrow=2) +
    scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
                       name="Substitution\nerror",
                       labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    geom_line(aes(color=Error)) +
    geom_point(aes(color=Error)) +
    scale_x_continuous(breaks=seq(11, 99, 10)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F Score")

poisson.pr.max <- poisson.pr %>% 
    group_by(Coverage, Error) %>% 
    top_n(1, FScore) %>%
    arrange(Coverage, Error, kmer) %>%
    top_n(1, -kmer)

poisson.pr.max$Error <- factor(poisson.pr.max$Error, levels=c(0, 0.001, 0.002, 0.005))


ggplot(poisson.pr %>% filter(kmer > 9), aes(kmer, FScore)) +
    theme_bw() +
    facet_wrap(~ Coverage, nrow=2) +
    scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
                       name="Substitution\nerror",
                       labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    geom_line(aes(color=Error)) +
    geom_point(data=poisson.pr.max, aes(color=Error)) +
    scale_x_continuous(breaks=seq(9, 99, 10)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F Score")





ggplot(poisson.pr.max, aes(kmer, FScore)) +
    theme_bw() +
    facet_wrap(~ Error, nrow=2) +
    #scale_color_manual(#values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
                       #name="Substitution\nerror",
                       #labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    scale_color_discrete(name="Coverage") +
    #geom_line(aes(color=as.factor(Coverage))) +
    geom_path(alpha=0.6) +
    geom_point(size=3, alpha=0.6, aes(color=as.factor(Coverage))) +
    scale_x_continuous(breaks=seq(11, 99, 10)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F Score")

poisson.pr.sub <- poisson.pr.base %>%
    filter(kmer>=43 & kmer<=81) %>%
    arrange(desc(Error))

poisson.pr.sub$Error <- factor(poisson.pr.sub$Error, levels=c(0.005, 0.002, 0.001, 0))

ggplot(poisson.pr.sub, aes(Coverage, FScore)) +
    theme_bw() +
    facet_wrap(~ kmer) +
    geom_point(aes(color=Error)) +
    geom_line(aes(color=Error)) +
    scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
                       name="Substitution\nerror",
                       labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    scale_x_continuous(breaks=seq(20, 200, 20)) +
    scale_y_continuous(breaks=seq(0.94, 1, 0.01)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12),
          panel.grid.minor=element_blank()) +
    ylab("Average F Score")


poisson.pr.each <- read.csv("counts-pr-each.csv", header=F, stringsAsFactors=F,
                            col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                        "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))

poisson.pr.each.sub <- poisson.pr.each %>%
    filter(kmer==95) %>%
    arrange(desc(Error))

poisson.pr.each.sub$Error <- factor(poisson.pr.each.sub$Error, levels=c(0.005, 0.002, 0.001, 0))
poisson.pr.each.sub$Allele <- factor(poisson.pr.each.sub$Allele, levels=allele.order)


poisson.pr.each$Error <- factor(poisson.pr.each$Error, levels=c(0.005, 0.002, 0.001, 0))
poisson.pr.each$Allele <- factor(poisson.pr.each$Allele, levels=allele.order)



ggplot(poisson.pr.each, aes(kmer, FScore)) +
    theme_bw() +
    facet_wrap(~ Allele) +
    geom_point(aes(color=Error)) +
    geom_line(aes(color=Error)) +
    #scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
    #                   name="Substitution\nerror",
    #                   labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    #scale_x_continuous(breaks=seq(20, 200, 20)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12),
          panel.grid.minor=element_blank()) +
    ylab("Average F Score")

###############
# Heatmaps
###############
poisson.max <- read.csv("called/all-called-lam-err-redo.csv", stringsAsFactors=F,
                        col.names=c("Simulated", "Tested", "Score", 
                                    "Correct", "Duplicated", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Simulated)