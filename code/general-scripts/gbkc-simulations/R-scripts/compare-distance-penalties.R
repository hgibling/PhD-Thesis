library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lattice)



setwd("~/prdm9/simulations")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))


combined <- read.csv("penalties/penalties-pr.csv", header=F, stringsAsFactors=F,
                     col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                 "Iteration", "ReadLength", "FragmentLength", "Method", "Penalty")) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method, Penalty) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    filter(kmer <= ReadLength)

each <- read.csv("penalties/penalties-pr-each.csv", header=F, stringsAsFactors=F,
                 col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", 
                             "kmer", "Iteration", "ReadLength", "FragmentLength", "Method", "Penalty")) %>%
    group_by(ReadLength, FragmentLength, Allele, Coverage, Error, kmer, Method, Penalty) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) 

counts <- read.csv("counts/counts-flank10k-pr.csv", header=F, stringsAsFactors=F,
                       col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                   "Iteration", "ReadLength", "FragmentLength", "Method", "Flank"))  %>% 
    select(-Flank) %>%
    filter(Coverage==20) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    mutate(Penalty1=1, Penalty10=10, Penalty100=100, Penalty1000=1000,
           Penalty10000=10000, Penalty100000=100000, Penalty1000000=1000000, Penalty10000000=10000000) %>%
    gather("Penalty", "Value", 10:17) %>%
    mutate(Penalty=Value) %>%
    select(-Value)
  

counts.each <- read.csv("counts/counts-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                   col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", "kmer", 
                               "Iteration", "ReadLength", "FragmentLength", "Method", "Flank"))  %>% 
    select(-Flank) %>%
    filter(Coverage==20) %>%
    group_by(Allele, ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    mutate(Penalty1=1, Penalty10=10, Penalty100=100, Penalty1000=1000,
           Penalty10000=10000, Penalty100000=100000, Penalty1000000=1000000, Penalty10000000=10000000) %>%
    gather("Penalty", "Value", 11:18) %>%
    mutate(Penalty=Value) %>%
    select(-Value)

distances.counts <- bind_rows(combined, counts)
distances.counts$Method <- factor(distances.counts$Method, levels=unique(distances.counts$Method))

distances.counts.each <- bind_rows(each, counts.each)
distances.counts.each$Method <- factor(distances.counts.each$Method, levels=unique(distances.counts.each$Method))

# Plots

cols <- c("red", "blue", "black")

ggplot(distances.counts, aes(kmer, FScore)) +
    facet_grid(Error ~ Penalty) +
    theme_bw() +
    geom_point(aes(color=Method)) +
    geom_line(aes(color=Method)) +
    scale_color_manual(values=cols) +
    theme(axis.text.x=element_text(size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("20X coverage")

ggplot(distances.counts, aes(kmer, FScore)) +
    facet_grid(Error ~ Penalty) +
    theme_bw() +
    geom_point(aes(color=Method)) +
    geom_line(aes(color=Method)) +
    scale_color_manual(values=cols) +
    coord_cartesian(ylim=c(0.5, 1)) +
    theme(axis.text.x=element_text(size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("20X coverage")

ggplot(distances.counts %>% filter(Penalty %in% c(1, 10)) %>% mutate(Penalty=as.factor(Penalty)), aes(kmer, FScore)) +
    facet_grid(Error ~ Method) +
    theme_bw() +
    geom_point(aes(color=Penalty)) +
    geom_line(aes(color=Penalty)) +
    #scale_color_manual(values=cols) +
    coord_cartesian(ylim=c(0.5, 1)) +
    theme(axis.text.x=element_text(size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("20X coverage")

ggplot(distances.counts.each %>% filter(Allele=="A"), aes(kmer, FScore)) +
    facet_grid(Error ~ Penalty) +
    theme_bw() +
    geom_point(aes(color=Method)) +
    geom_line(aes(color=Method)) +
    scale_color_manual(values=cols) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F1 Score") +
    ggtitle("Allele A, 20X coverage")

ggplot(distances.counts.each %>% filter(Allele=="A") %>% mutate(Penalty=as.factor(Penalty)), aes(kmer, FScore)) +
    facet_grid(Error ~ Method) +
    theme_bw() +
    geom_point(aes(color=Penalty)) +
    geom_line(aes(color=Penalty)) +
    #scale_color_manual(values=cols) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F1 Score") +
    ggtitle("Allele A, 20X coverage")

ggplot(distances.counts.each %>% filter(Allele=="A", Penalty %in% c(1, 10, 1000, 10000)) %>% mutate(Penalty=as.factor(Penalty)), aes(kmer, FScore)) +
    facet_grid(Error ~ Method) +
    theme_bw() +
    geom_point(aes(color=Penalty)) +
    geom_line(aes(color=Penalty)) +
    #scale_color_manual(values=cols) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12),
          axis.title.x=element_text(size=12),
          axis.text.y=element_text(size=12),
          axis.title.y=element_text(size=12),
          legend.text=element_text(size=12),
          legend.title=element_text(size=12),
          strip.text.x=element_text(size=12)) +
    ylab("Average F1 Score") +
    ggtitle("Allele A, 20X coverage")



#######################
#
# confusion matrices
#
#######################

setwd("~/prdm9/simulations/")

counts.max <- read.csv("counts/counts-max.csv", header=F, stringsAsFactors=F, 
                       col.names=c("Simulated", "Tested", "Score", "Correct", "Duplicated", "Coverage", "Error",
                                   "kmer", "Iteration", "ReadLength", "FragmentLength", "Method")) %>%
    mutate(Method="kmer-counts")
distances.max <- read.csv("distances/distances-max.csv", header=F, stringsAsFactors=F, 
                          col.names=c("Simulated", "Tested", "Score", "Correct", "Duplicated", "Coverage", "Error",
                                      "kmer", "Iteration", "ReadLength", "FragmentLength", "Method")) %>%
    mutate(Method=ifelse(Method=="max", "distance-max", "distance-sum"))
counts.distances.max <- bind_rows(counts.max, distances.max)

cm.counts <- counts.distances.max %>%
    select(-ReadLength, -FragmentLength) %>%
    group_by(Simulated, Iteration, Method, Coverage, Error, kmer) %>%
    count()

cm.max.sum <- counts.distances.max %>%
    full_join(cm.counts, by=c("Simulated", "Iteration",  "Method", "Coverage", "Error", "kmer")) %>%
    mutate(Count=1/n) %>%
    ungroup() %>%
    group_by(Simulated, Tested, Method, Coverage, Error, kmer) %>%
    summarize(Sum=sum(Count)) %>%
    ungroup()

make.cm <- function(dat, cov, e, m) {
    data.spread <- dat %>%
        filter(Method==UQ(m), Error==UQ(e), Coverage==UQ(cov)) %>%
        mutate(Simulated=factor(Simulated, levels=allele.order),
               Tested=factor(Tested, levels=allele.order)) %>%
        arrange(Simulated, Tested) %>%
        spread(Tested, Sum, fill=0) %>%
        select(-Simulated)

    data.missing.alleles <- allele.order[which(allele.order %in% colnames(data.spread)==F)]
    # missing columns because never called
    
    data.missing.alleles.dataframe <- data.frame(matrix(ncol=length(data.missing.alleles), nrow=length(allele.order), 0))
    colnames(data.missing.alleles.dataframe) <- data.missing.alleles
    
    data.spread.complete <- data.spread %>%
        bind_cols(data.missing.alleles.dataframe) %>%
        select(allele.order)
    
    data.spread.mat <- data.matrix(data.spread.complete)
    
    rownames(data.spread.mat) <- colnames(data.spread.mat)
    
    #pdf(paste(cov, "x", e, "e", m, "m.pdf", sep="-"))
    print(levelplot(t(data.spread.mat), 
              scale=list(x=list(rot=45)),
              col.regions=colorRampPalette(c("darkgreen","white","#D55E00"))(1e2),
              ylab="Simulated",
              xlab="Called",
              main=paste0("Confusion matrix for ", cov, "X, ", e, " error, k=51 using ", m)))
    #dev.off()
}

pdf("confusion-matrices.pdf")
for (coverage in unique(cm.max.sum$Coverage)){
    for (error in unique(cm.max.sum$Error)) {
        for (method in unique(cm.max.sum$Method)) {
            make.cm(cm.max.sum, coverage, error, method)
        }
    }
}
dev.off()
