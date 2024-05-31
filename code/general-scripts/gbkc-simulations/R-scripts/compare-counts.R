library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lattice)


setwd("~/prdm9/simulations/counts")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))


combined <- read.csv("counts-flank10k-pr.csv", header=F, stringsAsFactors=F,
                     col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                 "Iteration", "ReadLength", "FragmentLength", "Method", "Flank")) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    filter(kmer <= ReadLength) %>%
    arrange(Coverage)

each <- read.csv("counts-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                 col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", "kmer", 
                             "Iteration", "ReadLength", "FragmentLength", "Method", "Flank")) %>%
    group_by(ReadLength, FragmentLength, Allele, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))


combined$Coverage <- factor(combined$Coverage, levels=seq(20, 100, 20))
combined$Error <- factor(combined$Error, levels=c(0, 0.001, 0.01))

each$Allele <- factor(each$Allele, levels=allele.order)
each$Coverage <- factor(each$Coverage, levels=seq(20, 100, 20))
each$Error <- factor(each$Error, levels=c(0.01, 0.001, 0))


# Plots

greens.error <- brewer.pal(9, "Greens")[c(9, 7, 5)]
greens.coverage <- brewer.pal(9, "Greens")[5:9]

set1.coverage <- brewer.pal(9, "Set1")[c(2, 5, 4, 3, 1)]
dark2.coverage <- brewer.pal(8, "Dark2")[c(1:4, 6)]


ggplot(combined, aes(kmer, FScore)) +
    theme_bw() +
    facet_wrap(~ Coverage) +
    #scale_color_manual(values=rev(greens.error), breaks=c(0, 0.001, 0.01),
    #                   name="Substitution\nerror",
    #                   labels=c("0%", "0.1%", "1%")) +
    geom_point(aes(color=Error)) +
    geom_line(aes(color=Error)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("Count scores for 10k flank")

error.lab <- c("0% error", "0.1% error", "1% error")
names(error.lab) <- c(0, 0.001, 0.01)

ggplot(combined, aes(kmer, FScore)) +
    theme_bw() +
    facet_wrap(~ Error, labeller=labeller(Error=error.lab)) +
    scale_color_manual(values=(set1.coverage), breaks=rev(seq(20, 100, 20)), drop=F) +
    geom_point(aes(color=Coverage), size=2) +
    geom_line(aes(color=Coverage)) +
    theme(axis.text.x=element_text(size=20),
          axis.title.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.title.y=element_text(size=20),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20),
          strip.text.x=element_text(size=20)) +
    ylab("") +
    xlab("") +
    ggtitle("")

ggplot(combined, aes(kmer, FScore)) +
    theme_bw() +
    facet_wrap(~ Error, labeller=labeller(Error=error.lab)) +
    scale_color_manual(values=(set1.coverage), breaks=rev(seq(20, 100, 20)), drop=F) +
    geom_point(aes(color=Coverage), size=2, alpha=0) +
    #geom_line(aes(color=Coverage)) +
    theme(axis.text.x=element_text(size=20),
          axis.title.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.title.y=element_text(size=20),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20),
          strip.text.x=element_text(size=20)) +
    ylab("") +
    xlab("") +
    ggtitle("")

ggplot(combined, aes(kmer, FScore)) +
    theme_bw() +
    facet_wrap(~ Error, labeller=labeller(Error=error.lab)) +
    scale_color_manual(values=(set1.coverage), breaks=rev(seq(20, 100, 20)), drop=F) +
    geom_point(aes(color=Coverage), size=2) +
    geom_line(aes(color=Coverage)) +
    theme(axis.text.x=element_text(size=20),
          axis.title.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.title.y=element_text(size=20),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20),
          strip.text.x=element_text(size=20)) +
    coord_cartesian(ylim=c(0.8, 1)) +
    ylab("") +
    xlab("") +
    ggtitle("")



ggplot(combined, aes(as.numeric(levels(Coverage))[Coverage], FScore)) +
    theme_bw() +
    facet_wrap(~ kmer) +
    scale_color_manual(values=rev(greens.error), breaks=c(0, 0.001, 0.01),
                       name="Substitution\nerror",
                       labels=c("0%", "0.1%", "1%")) +
    geom_point(aes(color=Error)) +
    geom_line(aes(color=Error)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    xlab("Coverage") +
    ggtitle("Count scores for 10k flank")

ggplot(each %>% filter(Coverage==100, Error==0), aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
    theme_bw() +
    geom_point(aes(color=Method)) +
    geom_line(aes(color=Method)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F1 Score") +
    ggtitle("100X coverage, 0% error rate, 10k flank")


########## heatmaps

counts.max <- read.csv("counts-max.csv", stringsAsFactors=F,
                   header=F, col.names=c("Simulated", "Tested", "Score", "Correct", "Duplicated", "Coverage", 
                                          "Error", "kmer", "Iteration", "ReadLength", "FragmentLength", "Method")) 

counts.max.count <- counts.max %>%
    group_by(Simulated, Iteration, Coverage, Error, kmer, Method) %>%
    count()

counts.max.sum <- counts.max %>%
    full_join(counts.max.count, 
              by=c("Simulated", "Iteration", "Coverage", "Error", "kmer", "Method")) %>%
    mutate(Count=1/n) %>%
    ungroup() %>%
    group_by(Simulated, Tested, Coverage, Error, kmer, Method) %>%
    summarize(Sum=sum(Count))

av.counts.max.sum.spread <- counts.max.sum %>%
    filter(Coverage==20, Error==0.1, kmer==51) %>%
    ungroup() %>%
    mutate(Simulated=factor(Simulated, levels=allele.order),
           Tested=factor(Tested, levels=allele.order)) %>%
    arrange(Simulated, Tested) %>%
    spread(Tested, Sum, fill=0) %>%
    select(-Simulated)

av.missing.alleles <- allele.order[which(allele.order %in% colnames(av.counts.max.sum.spread)==F)]
# missing columns because never called

av.missing.alleles.dataframe <- data.frame(matrix(ncol=length(av.missing.alleles), nrow=length(allele.order), 0))
colnames(av.missing.alleles.dataframe) <- av.missing.alleles

av.counts.max.sum.spread.complete <- av.counts.max.sum.spread %>%
    bind_cols(av.missing.alleles.dataframe) %>%
    select(allele.order)

av.counts.max.sum.spread.mat <- data.matrix(av.counts.max.sum.spread.complete)

rownames(av.counts.max.sum.spread.mat) <- colnames(av.counts.max.sum.spread.mat)

levelplot(t(av.counts.max.sum.spread.mat), 
          scale=list(x=list(rot=45, cex=1), y=list(cex=1)),
          col.regions=colorRampPalette(c("darkgreen","white","#D55E00"))(100),
          colorkey=list(at=seq(0,100,10), labels=list(cex=1)),
          ylab="Simulated",
          xlab="Called")
