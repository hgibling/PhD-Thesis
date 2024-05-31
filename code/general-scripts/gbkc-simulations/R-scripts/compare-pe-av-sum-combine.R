library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)
library(dplyr)


setwd("~/awadalla/PRDM9-Project/HMM/sim-reads/combine")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

av.count.combined <- read.csv("av-pe-pr-combined.csv", header=F, stringsAsFactors=F,
                              col.names=c("Precision", "Recall", "Coverage", 
                                          "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="AverageCount") %>%
    arrange(desc(Error))

sum.count.combined <- read.csv("sum-pe-pr-combined.csv", header=F, stringsAsFactors=F,
                               col.names=c("Precision", "Recall", "Coverage", 
                                           "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="SumCount") %>%
    arrange(desc(Error))

av.combined <- read.csv("~/awadalla/PRDM9-Project/HMM/sim-reads/pe-lam/av-pe-pr.csv", header=F, stringsAsFactors=F,
                        col.names=c("Precision", "Recall", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="Average") %>%
    arrange(desc(Error))

sum.combined <- read.csv("~/awadalla/PRDM9-Project/HMM/sim-reads/pe-lam/sum-pe-pr.csv", header=F, stringsAsFactors=F,
                         col.names=c("Precision", "Recall", "Coverage", 
                                     "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="Sum") %>%
    arrange(desc(Error))

counts.combined <- read.csv("~/awadalla/PRDM9-Project/HMM/sim-reads/se-counts/se-counts-pr.csv", 
                            header=F, stringsAsFactors=F,
                            col.names=c("Precision", "Recall", "Coverage", 
                                        "Error", "kmer", "Iteration")) %>%
    filter(Error==0 & Coverage==100) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="Count") %>%
    arrange(desc(Error))

combined <- bind_rows(av.count.combined, sum.count.combined, av.combined, sum.combined, counts.combined)

## 
av.count.each <- read.csv("av-pe-pr-each-combined.csv", header=F, stringsAsFactors=F,
                    col.names=c("Allele", "Precision", "Recall", "Coverage", 
                                "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="AverageCount") %>%
    arrange(desc(Error))

sum.count.each <- read.csv("sum-pe-pr-each-combined.csv", header=F, stringsAsFactors=F,
                     col.names=c("Allele", "Precision", "Recall", "Coverage", 
                                 "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="SumCount") %>%
    arrange(desc(Error))

av.each <- read.csv("~/awadalla/PRDM9-Project/HMM/sim-reads/pe-lam/av-pe-pr-each.csv", header=F, stringsAsFactors=F,
                          col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                      "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="Average") %>%
    arrange(desc(Error))

sum.each <- read.csv("~/awadalla/PRDM9-Project/HMM/sim-reads/pe-lam/sum-pe-pr-each.csv", header=F, stringsAsFactors=F,
                           col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                       "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="Sum") %>%
    arrange(desc(Error))

counts.each <- read.csv("~/awadalla/PRDM9-Project/HMM/sim-reads/se-counts/se-counts-pr-each.csv", 
                        header=F, stringsAsFactors=F,
                        col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                    "Error", "kmer", "Iteration")) %>%
    filter(Error==0 & Coverage==100) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore), Scoring="Count") %>%
    arrange(desc(Error))

each <- bind_rows(av.count.each, sum.count.each, av.each, sum.each, counts.each)

combined$Error <- factor(combined$Error, levels=c(0.005, 0.002, 0.001, 0))
combined$Scoring <- factor(combined$Scoring, levels=c("AverageCount", "SumCount", "Average", "Sum", "Count"))
each$Error <- factor(each$Error, levels=c(0.005, 0.002, 0.001, 0))
each$Scoring <- factor(each$Scoring, levels=c("AverageCount", "SumCount", "Average", "Sum", "Count"))
each$Allele <- factor(each$Allele, levels=allele.order)


# Plots

gg.color.hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

cols <- c(gg.color.hue(4), "black")

ggplot(combined, aes(kmer, FScore)) +
    theme_bw() +
    geom_point(aes(color=Scoring)) +
    geom_line(aes(color=Scoring)) +
    scale_color_manual(values=cols) +
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
    geom_point(aes(color=Scoring)) +
    geom_line(aes(color=Scoring)) +
    scale_color_manual(values=cols) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")

ggplot(each %>% filter(kmer>=87), aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=Scoring)) +
    geom_line(aes(color=Scoring)) +
    scale_x_continuous(breaks=seq(75, 99, 4)) +
    scale_color_manual(values=cols) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F Score") +
    ggtitle("100X coverage, 0% error rate")



### Confusion matrix

pe.max <- read.csv("meta/meta-flank1000-100-bp-100-x-e-0-F-pe-99-k.csv", stringsAsFactors=F,
                   header=F, col.names=c("Simulated", "Tested", "Score", 
                                         "Correct", "Duplicated", "Iteration", "Scoring")) 

pe.max.count <- pe.max %>%
    group_by(Simulated, Iteration, Scoring) %>%
    count()

pe.max.sum <- pe.max %>%
    full_join(pe.max.count, 
              by=c("Simulated", "Iteration", "Scoring")) %>%
    mutate(Count=1/n) %>%
    ungroup() %>%
    group_by(Simulated, Tested, Scoring) %>%
    summarize(Sum=sum(Count))

av.pe.max.sum.spread <- pe.max.sum %>%
    filter(Scoring=="Average") %>%
    ungroup() %>%
    mutate(Simulated=factor(Simulated, levels=allele.order),
           Tested=factor(Tested, levels=allele.order)) %>%
    arrange(Simulated, Tested) %>%
    spread(Tested, Sum, fill=0) %>%
    select(-Simulated)

av.missing.alleles <- allele.order[which(allele.order %in% colnames(av.pe.max.sum.spread)==F)]
# missing columns because never called

av.missing.alleles.dataframe <- data.frame(matrix(ncol=length(av.missing.alleles), nrow=length(allele.order), 0))
colnames(av.missing.alleles.dataframe) <- av.missing.alleles

av.pe.max.sum.spread.complete <- av.pe.max.sum.spread %>%
    bind_cols(av.missing.alleles.dataframe) %>%
    select(allele.order)

av.pe.max.sum.spread.mat <- data.matrix(av.pe.max.sum.spread.complete)

rownames(av.pe.max.sum.spread.mat) <- colnames(av.pe.max.sum.spread.mat)

levelplot(t(av.pe.max.sum.spread.mat), 
          scale=list(x=list(rot=45)),
          col.regions=colorRampPalette(c("darkgreen","white","#D55E00"))(1e2),
          ylab="Simulated",
          xlab="Called")


sum.pe.max.sum.spread <- pe.max.sum %>%
    filter(Scoring=="Sum") %>%
    ungroup() %>%
    mutate(Simulated=factor(Simulated, levels=allele.order),
           Tested=factor(Tested, levels=allele.order)) %>%
    arrange(Simulated, Tested) %>%
    spread(Tested, Sum, fill=0) %>%
    select(-Simulated)

sum.missing.alleles <- allele.order[which(allele.order %in% colnames(sum.pe.max.sum.spread)==F)]
# missing columns because never called

sum.missing.alleles.dataframe <- data.frame(matrix(ncol=length(sum.missing.alleles), nrow=length(allele.order), 0))
colnames(sum.missing.alleles.dataframe) <- sum.missing.alleles

sum.pe.max.sum.spread.complete <- sum.pe.max.sum.spread %>%
    bind_cols(sum.missing.alleles.dataframe) %>%
    select(allele.order)

sum.pe.max.sum.spread.mat <- data.matrix(sum.pe.max.sum.spread.complete)

rownames(sum.pe.max.sum.spread.mat) <- colnames(sum.pe.max.sum.spread.mat)

levelplot(t(sum.pe.max.sum.spread.mat), 
          scale=list(x=list(rot=45)),
          col.regions=colorRampPalette(c("darkgreen","white","#D55E00"))(1e2),
          ylab="Simulated",
          xlab="Called")
