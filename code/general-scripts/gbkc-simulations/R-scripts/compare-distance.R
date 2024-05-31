library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lattice)



setwd("~/prdm9/simulations")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))


combined <- read.csv("distances/distance-flank10k-pr.csv", header=F, stringsAsFactors=F,
                     col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                 "Iteration", "ReadLength", "FragmentLength", "Method", "Flank")) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method, Flank) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    filter(kmer <= ReadLength)

each <- read.csv("distances/distance-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                 col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", 
                             "kmer", "Iteration", "ReadLength", "FragmentLength", "Method", "Flank")) %>%
    group_by(ReadLength, FragmentLength, Allele, Coverage, Error, kmer, Method, Flank) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

counts <- read.csv("counts/counts-flank10k-pr.csv", header=F, stringsAsFactors=F,
                       col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                   "Iteration", "ReadLength", "FragmentLength", "Method", "Flank"))  %>% 
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method, Flank) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

counts.each <- read.csv("counts/counts-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                   col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", "kmer", 
                               "Iteration", "ReadLength", "FragmentLength", "Method", "Flank"))  %>% 
    group_by(Allele, ReadLength, FragmentLength, Coverage, Error, kmer, Method, Flank) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

distances.counts <- bind_rows(combined, counts) %>%
    ungroup() %>%
    mutate(Method=case_when(
        Method == "counts" ~ "Counts",
        Method == "max" ~ "Distance Max",
        Method == "sum" ~ "Distance Sum"
    ))
distances.counts$Method <- factor(distances.counts$Method, levels=c("Counts", "Distance Max", "Distance Sum"))

distances.counts.each <- bind_rows(each, counts.each)
distances.counts.each$Method <- factor(distances.counts.each$Method, levels=unique(distances.counts.each$Method))

# Plots

cols <- c("black", "#4daf4a", "#984ea3")

error.lab <- c("0% error", "1% error")
names(error.lab) <- c(0, 0.01)
cov.lab <- c("20X coverage", "100X coverage")
names(cov.lab) <- c(20, 100)

ggplot(distances.counts %>% filter(Coverage %in% c(20, 100), Error != 0.001), aes(kmer, FScore)) +
    facet_grid(Coverage ~ Error, labeller=labeller(Error=error.lab, Coverage=cov.lab)) +
    theme_bw() +
    geom_point(aes(color=Method), size=3) +
    geom_line(aes(color=Method)) +
    scale_color_manual(values=cols) +
    theme(axis.text.x=element_text(size=20),
          axis.title.x=element_text(size=20),
          axis.text.y=element_text(size=20),
          axis.title.y=element_text(size=20),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20),
          strip.text.x=element_text(size=20),
          strip.text.y=element_text(size=20)) +
    ylab("Average F1 Score") +
    ggtitle("")

ggplot(distances.counts.each %>% filter(Coverage==40, Error==0), aes(kmer, FScore)) +
    facet_wrap(~ Allele) +
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
    ggtitle("40X coverage, 0% error rate")


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
