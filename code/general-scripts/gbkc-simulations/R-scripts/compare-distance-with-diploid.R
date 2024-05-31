library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(lattice)



setwd("~/prdm9/simulations")

allele.order <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))


haploid <- read.csv("distances/distance-flank10k-pr.csv", header=F, stringsAsFactors=F,
                     col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                 "Iteration", "ReadLength", "FragmentLength", "Method", "Flank")) %>%
    mutate(Method=paste0(Method, "-Haploid")) %>%
    select(-Flank) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    filter(kmer <= ReadLength)

haploid.each <- read.csv("distances/distance-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                 col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", 
                             "kmer", "Iteration", "ReadLength", "FragmentLength", "Method", "Flank")) %>%
    mutate(Method=paste0(Method, "-Haploid")) %>%
    select(-Flank) %>%
    group_by(ReadLength, FragmentLength, Allele, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

diploid <- read.csv("diploid-distances/distance-flank10k-pr.csv", header=F, stringsAsFactors=F,
                    col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                "Iteration", "ReadLength", "FragmentLength", "Method")) %>%
    mutate(Method=paste0(Method, "-Diploid")) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore)) %>%
    filter(kmer <= ReadLength)

diploid.each <- read.csv("diploid-distances/distance-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                         col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", 
                                     "kmer", "Iteration", "ReadLength", "FragmentLength", "Method")) %>%
    mutate(Method=paste0(Method, "-Diploid")) %>%
    group_by(ReadLength, FragmentLength, Allele, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

counts <- read.csv("counts/counts-flank10k-pr.csv", header=F, stringsAsFactors=F,
                       col.names=c("Precision", "Recall", "Coverage", "Error", "kmer", 
                                   "Iteration", "ReadLength", "FragmentLength", "Method", "Flank"))  %>% 
    mutate(Method=paste0(Method, "-Haploid")) %>%
    select(-Flank) %>%
    group_by(ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

counts.each <- read.csv("counts/counts-flank10k-pr-each.csv", header=F, stringsAsFactors=F,
                   col.names=c("Allele", "Precision", "Recall", "Coverage", "Error", "kmer", 
                               "Iteration", "ReadLength", "FragmentLength", "Method", "Flank"))  %>% 
    mutate(Method=paste0(Method, "-Haploid")) %>%
    select(-Flank) %>%
    group_by(Allele, ReadLength, FragmentLength, Coverage, Error, kmer, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall)) %>%
    mutate(FScore=2 * ((MeanPrecision * MeanRecall) / (MeanPrecision + MeanRecall))) %>%
    mutate(FScore=ifelse(is.nan(FScore), 0, FScore))

distances.counts <- bind_rows(haploid, diploid, counts)
distances.counts$Error <- factor(distances.counts$Error, levels=unique(distances.counts$Error))
distances.counts$Coverage <- factor(distances.counts$Coverage, levels=unique(distances.counts$Coverage))
distances.counts$Method <- factor(distances.counts$Method, levels=unique(distances.counts$Method))


distances.counts.each <- bind_rows(haploid.each, diploid.each, counts.each)
distances.counts.each$Method <- factor(distances.counts.each$Method, levels=unique(distances.counts.each$Method))
distances.counts.each$Allele <- factor(distances.counts.each$Allele, levels=unique(distances.counts.each$Allele))

# Plots

cols <- c("red", "blue", "pink", "lightblue", "black")

ggplot(distances.counts, aes(kmer, FScore)) +
    facet_grid(Error ~ Coverage) +
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
    ggtitle("")

ggplot(distances.counts, aes(kmer, FScore)) +
    facet_grid(Error ~ Coverage) +
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
    coord_cartesian(ylim=c(0.75, 1)) +
    ylab("Average F1 Score") +
    ggtitle("")

distances.counts.each <- distances.counts.each %>%
    ungroup() %>%
    mutate(Allele=gsub("/.*", "", Allele))

ggplot(distances.counts.each %>% filter(Coverage==80, Error==0), aes(kmer, FScore)) +
    facet_wrap( ~ Allele) +
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
    ggtitle("80X coverage, 0% error rate")


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



##

sub <- distances.counts %>% filter(Error==0)

> plot_ly(z = ~volcano) %>% add_surface()


plot_ly(sub, x=~Coverage, y=~kmer, z=~FScore, color=~Method) %>% add_markers()



hc <- sub %>% filter(Method=="counts-Haploid") %>%
    ungroup() %>%
    select(Coverage, kmer, FScore) %>%
    spread(kmer, FScore) %>%
    as.data.frame()

hc <- hc[,-1]
colnames(hc) <- NULL
mhc <- as.matrix(hc)

hmax <- sub %>% filter(Method=="max-Haploid") %>%
    ungroup() %>%
    select(Coverage, kmer, FScore) %>%
    spread(kmer, FScore) %>%
    as.data.frame()

hmax <- hmax[,-1]
colnames(hmax) <- NULL
mhmax <- as.matrix(hmax)

hsum <- sub %>% filter(Method=="sum-Haploid") %>%
    ungroup() %>%
    select(Coverage, kmer, FScore) %>%
    spread(kmer, FScore) %>%
    as.data.frame()

hsum <- hsum[,-1]
colnames(hsum) <- NULL
mhsum <- as.matrix(hsum)

dmax <- sub %>% filter(Method=="max-Diploid") %>%
    ungroup() %>%
    select(Coverage, kmer, FScore) %>%
    spread(kmer, FScore) %>%
    as.data.frame()

dmax <- dmax[,-1]
colnames(dmax) <- NULL
mdmax <- as.matrix(dmax)

dsum <- sub %>% filter(Method=="sum-Diploid") %>%
    ungroup() %>%
    select(Coverage, kmer, FScore) %>%
    spread(kmer, FScore) %>%
    as.data.frame()

dsum <- dsum[,-1]
colnames(dsum) <- NULL
mdsum <- as.matrix(dsum)


plot_ly(showscale=F) %>% 
    add_surface(z = ~mhc, opacity= 0.5, colorscale=list(c(0,1),c("#36333d", "#635d70"))) %>%        #black
    add_surface(z = ~mhmax, opacity= 0.5, colorscale=list(c(0,1),c("#5b9eeb","#3676bf"))) %>%       #yellow
    add_surface(z = ~mhsum, opacity= 0.5, colorscale=list(c(0,1),c("#ed535d","#ab242d"))) %>%       #red
    add_surface(z = ~mdmax, opacity= 0.5, colorscale=list(c(0,1),c("826de8", "#4d39ad"))) %>%       #purple
    add_surface(z = ~mdsum, opacity= 0.5, colorscale=list(c(0,1),c("#517143","#7FB069")))           #green
