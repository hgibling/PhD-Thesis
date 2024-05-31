library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(lattice)


setwd("~/awadalla/PRDM9-Project/HMM/sim-reads/counts")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

pe.combined <- read.csv("100/meta-counts-pr.csv", header=F, stringsAsFactors=F,
                            col.names=c("Precision", "Recall", "Coverage", 
                                        "Error", "kmer", "Iteration")) %>%
    group_by(Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))

pe.each <- read.csv("100/meta-counts-pr-each.csv", header=F, stringsAsFactors=F,
                            col.names=c("Allele", "Precision", "Recall", "Accuracy", "Coverage", 
                                        "Error", "kmer", "Iteration")) %>%
    group_by(Allele, Coverage, Error, kmer) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    arrange(desc(Error))

se.combined

se.each


pe.combined$Error <- factor(pe.combined$Error, levels=c(0.005, 0.002, 0.001, 0))
pe.each$Error <- factor(pe.each$Error, levels=c(0.005, 0.002, 0.001, 0))
pe.each$Allele <- factor(pe.each$Allele, levels=allele.order)


greens <- brewer.pal(9, "Greens")[c(9, 7, 5, 4)]

ggplot(pe.combined, aes(Coverage, FScore)) +
    theme_bw() +
    facet_wrap(~ kmer) +
    #scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
    #                   name="Substitution\nerror",
    #                   labels=c("0%", "0.1%", "0.2%", "0.5%")) +
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

ggplot(pe.each, aes(Allele, FScore)) +
    theme_bw() +
    #facet_wrap(~ Allele) +
    geom_point(aes(color=kmer)) +
    #geom_line(aes(color=Error)) +
    #scale_color_manual(values=rev(greens), breaks=c(0, 0.001, 0.002, 0.005),
                       #name="Substitution\nerror",
                       #labels=c("0%", "0.1%", "0.2%", "0.5%")) +
    #scale_x_continuous(breaks=seq(20, 200, 20)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12),
          panel.grid.minor=element_blank()) +
    guides(color=F) +
    ylab("Average F Score")


### Heatmaps

pe.max <- read.csv("meta/meta-flank1000-100-bp-100-x-e-0-F-pe-51-k.csv", stringsAsFactors=F,
                   header=F, col.names=c("Simulated", "Tested", "Score", 
                               "Correct", "Duplicated", "Iteration")) 

pe.max.count <- pe.max %>%
    group_by(Simulated, Iteration) %>%
    count()

pe.max.sum <- pe.max %>%
    full_join(pe.max.count, 
              by=c("Simulated", "Iteration")) %>%
    mutate(Count=1/n) %>%
    ungroup() %>%
    group_by(Simulated, Tested) %>%
    summarize(Sum=sum(Count))

pe.max.sum.spread <- pe.max.sum %>%
    ungroup() %>%
    mutate(Simulated=factor(Simulated, levels=allele.order),
           Tested=factor(Tested, levels=allele.order)) %>%
    arrange(Simulated, Tested) %>%
    spread(Tested, Sum, fill=0) %>%
    select(-Simulated)

missing.alleles <- allele.order[which(allele.order %in% colnames(pe.max.sum.spread)==F)]
# missing columns because never called

missing.alleles.dataframe <- data.frame(matrix(ncol=length(missing.alleles), nrow=length(allele.order), 0))
colnames(missing.alleles.dataframe) <- missing.alleles

pe.max.sum.spread.complete <- pe.max.sum.spread %>%
    bind_cols(missing.alleles.dataframe) %>%
    select(allele.order)

pe.max.sum.spread.mat <- data.matrix(pe.max.sum.spread.complete)

rownames(pe.max.sum.spread.mat) <- colnames(pe.max.sum.spread.mat)

levelplot(t(pe.max.sum.spread.mat), 
          scale=list(x=list(rot=45)),
          col.regions=colorRampPalette(c("darkgreen","white","#D55E00"))(1e2),
          ylab="Simulated (20X, 0% error)",
          xlab="Called")
