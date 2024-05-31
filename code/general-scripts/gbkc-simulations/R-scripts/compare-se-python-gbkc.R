library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(dplyr)


setwd("~/awadalla/sim-reads")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

py.combined <- read.csv("se-counts/se-counts-pr.csv", header=F, stringsAsFactors=F,
                     col.names=c("Precision", "Recall", "Coverage", 
                                 "Error", "kmer", "Iteration")) %>%
    filter(Coverage==100, Error==0) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="python")

py.each <- read.csv("se-counts/se-counts-pr-each.csv", header=F, stringsAsFactors=F,
                 col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                             "Error", "kmer", "Iteration")) %>%
    filter(Coverage==100, Error==0) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="python")

gbkc.combined <- read.csv("gbkc/gbkc-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="gbkc")

gbkc.each <- read.csv("gbkc/gbkc-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="gbkc")

combined <- bind_rows(py.combined, gbkc.combined)
combined$Type <- factor(combined$Type, levels=c("gbkc", "python"))

each <-bind_rows(py.each, gbkc.each)
each$Type <- factor(each$Type, levels=c("gbkc", "python"))
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")


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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")


####################3
# comapre wgsim and vg sim

vg.combined <- gbkc.combined %>%
    mutate(Type="vg")

vg.each <- gbkc.each %>%
    mutate(Type="vg")


wgsim.combined <- read.csv("gbkc-wgsim/gbkc-pr.csv", header=F, stringsAsFactors=F,
                          col.names=c("Precision", "Recall", "Coverage", 
                                      "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="wgsim")

wgsim.each <- read.csv("gbkc-wgsim/gbkc-pr-each.csv", header=F, stringsAsFactors=F,
                      col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                  "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="wgsim")

combined <- bind_rows(vg.combined, wgsim.combined)
combined$Type <- factor(combined$Type, levels=c("vg", "wgsim"))

each <-bind_rows(vg.each, wgsim.each)
each$Type <- factor(each$Type, levels=c("vg", "wgsim"))
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")


##########
# vg forward only vs vg both

forward.combined <- read.csv("gbkc/gbkc-pr.csv", header=F, stringsAsFactors=F,
                          col.names=c("Precision", "Recall", "Coverage", 
                                      "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="Forward")

forward.each <- read.csv("gbkc/gbkc-pr-each.csv", header=F, stringsAsFactors=F,
                      col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                  "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="Forward")

both.combined <- read.csv("gbkc-vg/vg-gbkc-pr.csv", header=F, stringsAsFactors=F,
                             col.names=c("Precision", "Recall", "Coverage", 
                                         "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="Both")

both.each <- read.csv("gbkc-vg/vg-gbkc-pr-each.csv", header=F, stringsAsFactors=F,
                         col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                     "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="Both")


combined <- bind_rows(forward.combined, both.combined)
combined$Type <- factor(combined$Type, levels=c("Forward", "Both"))

each <-bind_rows(forward.each, both.each)
each$Type <- factor(each$Type, levels=c("Forward", "Both"))
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")



##########
# vg forward only vs vg both vs wgsim outer vs wgsim inner

forward.combined <- read.csv("gbkc/gbkc-pr.csv", header=F, stringsAsFactors=F,
                             col.names=c("Precision", "Recall", "Coverage", 
                                         "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="vgForwardSE")

forward.each <- read.csv("gbkc/gbkc-pr-each.csv", header=F, stringsAsFactors=F,
                         col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                     "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="vgForwardSE")

both.combined <- read.csv("gbkc-vg/vg-gbkc-pr.csv", header=F, stringsAsFactors=F,
                          col.names=c("Precision", "Recall", "Coverage", 
                                      "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="vgBothSE")

both.each <- read.csv("gbkc-vg/vg-gbkc-pr-each.csv", header=F, stringsAsFactors=F,
                      col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                  "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="vgBothSE")

outer.combined <- read.csv("gbkc-wgsim/gbkc-pr.csv", header=F, stringsAsFactors=F,
                             col.names=c("Precision", "Recall", "Coverage", 
                                         "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="wgsimPE")

outer.each <- read.csv("gbkc-wgsim/gbkc-pr-each.csv", header=F, stringsAsFactors=F,
                         col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                     "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="wgsimPE")

se.combined <- read.csv("gbkc-se/se-wgsim-pr.csv", header=F, stringsAsFactors=F,
                           col.names=c("Precision", "Recall", "Coverage", 
                                       "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="wgsimSE")

se.each <- read.csv("gbkc-se/se-wgsim-pr-each.csv", header=F, stringsAsFactors=F,
                       col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                   "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="wgsimSE")

vgpe.combined <- read.csv("gbkc-vg-pe/vg-pe-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="vgPE")

vgpe.each <- read.csv("gbkc-vg-pe/vg-pe-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="vgPE")

my.combined <- read.csv("gbkc-my/my-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="mySE")

my.each <- read.csv("gbkc-my/my-pr-each.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error)) %>%
    mutate(Type="mySE")




combined <- bind_rows(forward.combined, both.combined, vgpe.combined,
                      outer.combined, my.combined, se.combined)
combined$Type <- factor(combined$Type, levels=c("vgPE", "vgForwardSE", "vgBothSE",
                                                "wgsimPE", "wgsimSE", "mySE"))

each <- bind_rows(forward.each, both.each, vgpe.each,
                      outer.each, my.each, se.each)
each$Type <- factor(each$Type, levels=c("vgPE", "vgForwardSE", "vgBothSE",
                                                "wgsimPE", "wgsimSE", "mySE"))
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")


ggplot(combined %>% filter(!(Type %in% c("wgsimPE", "vgPE"))) %>%
           filter(kmer < 99), aes(kmer, FScore)) +
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(each %>% filter(!(Type %in% c("wgsimPE", "vgPE"))) %>%
           filter(kmer < 99), aes(kmer, FScore)) +
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(each %>% filter(!(Type %in% c("wgsimPE", "vgPE"))) %>%
           filter(Allele %in% c("A", "C", "L2", "L7", "L11", "L14", "L16", "L36")) %>%
           filter(kmer < 99), 
       aes(kmer, FScore)) +
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
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")
