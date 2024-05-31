library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

setwd("~/prdm9/simulations/diploid-all")

diploid.all.20 <- read.csv("diploid-distances-20.csv", header=F, stringsAsFactors=F,
                    col.names=c("kmer", "Genotype", "Precision", "Recall", "F1Score", "Coverage", "Error",
                                "Iteration", "ReadLength", "FragmentLength", "Method"))

diploid.average.20 <- diploid.all.20 %>%
    filter(Genotype == "Average") %>%
    group_by(kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

diploid.each.20 <- diploid.all.20 %>%
    filter(Genotype != "Average") %>%
    group_by(Genotype, kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

write.table(diploid.average.20, "diploid-20-average.csv", quote=F, row.names=F, col.names=F, sep=",")
write.table(diploid.each.20, "diploid-20-each.csv", quote=F, row.names=F, col.names=F, sep=",")


diploid.all.40 <- read.csv("diploid-distances-40.csv", header=F, stringsAsFactors=F,
                           col.names=c("kmer", "Genotype", "Precision", "Recall", "F1Score", "Coverage", "Error",
                                       "Iteration", "ReadLength", "FragmentLength", "Method"))

diploid.average.40 <- diploid.all.40 %>%
    filter(Genotype == "Average") %>%
    group_by(kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

diploid.each.40 <- diploid.all.40 %>%
    filter(Genotype != "Average") %>%
    group_by(Genotype, kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

write.table(diploid.average.40, "diploid-40-average.csv", quote=F, row.names=F, col.names=F, sep=",")
write.table(diploid.each.40, "diploid-40-each.csv", quote=F, row.names=F, col.names=F, sep=",")


diploid.all.60 <- read.csv("diploid-distances-60.csv", header=F, stringsAsFactors=F,
                           col.names=c("kmer", "Genotype", "Precision", "Recall", "F1Score", "Coverage", "Error",
                                       "Iteration", "ReadLength", "FragmentLength", "Method"))

diploid.average.60 <- diploid.all.60 %>%
    filter(Genotype == "Average") %>%
    group_by(kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

diploid.each.60 <- diploid.all.60 %>%
    filter(Genotype != "Average") %>%
    group_by(Genotype, kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

write.table(diploid.average.60, "diploid-60-average.csv", quote=F, row.names=F, col.names=F, sep=",")
write.table(diploid.each.60, "diploid-60-each.csv", quote=F, row.names=F, col.names=F, sep=",")


diploid.all.80 <- read.csv("diploid-distances-80.csv", header=F, stringsAsFactors=F,
                           col.names=c("kmer", "Genotype", "Precision", "Recall", "F1Score", "Coverage", "Error",
                                       "Iteration", "ReadLength", "FragmentLength", "Method"))

diploid.average.80 <- diploid.all.80 %>%
    filter(Genotype == "Average") %>%
    group_by(kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

diploid.each.80 <- diploid.all.80 %>%
    filter(Genotype != "Average") %>%
    group_by(Genotype, kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
    summarize(MeanPrecision=mean(Precision), MeanRecall=mean(Recall), MeanF1Score=mean(F1Score))

write.table(diploid.average.80, "diploid-80-average.csv", quote=F, row.names=F, col.names=F, sep=",")
write.table(diploid.each.80, "diploid-80-each.csv", quote=F, row.names=F, col.names=F, sep=",")


# cat together in bash


### clear workspace

average <- read.csv("diploid-all-average.csv", stringsAsFactors=F, 
                    col.names=c("kmer", "ReadLength", "FragmentLength", "Coverage", "Error", "Method",
                                "MeanPrecision", "MeanRecall", "MeanF1Score"))

ggplot(average, aes(kmer, MeanF1Score)) +
    facet_grid(Error ~ Coverage) +
    theme_bw() +
    geom_point(aes(color=Method)) +
    geom_line(aes(color=Method)) +
    theme(axis.text.x=element_text(size=15),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("")



### each

each <- read.csv("diploid-all-each.csv", stringsAsFactors=F, 
                    col.names=c("Genotype", "kmer", "ReadLength", "FragmentLength", "Coverage", "Error", "Method",
                                "MeanPrecision", "MeanRecall", "MeanF1Score"))

ggplot(each %>% filter(Coverage==80, Error==0), aes(Genotype, MeanF1Score)) +
    facet_wrap(~ kmer) +
    theme_bw() +
    geom_point(aes(color=Method)) +
    geom_line(aes(color=Method)) +
    theme(axis.text.x=element_text(size=15, angle=45, hjust=1),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("80X coverage, 0% error rate")

ggplot(each %>% filter(Coverage==80, Error==0, MeanF1Score>0.5, kmer==99), aes(Genotype, MeanF1Score)) +
    facet_wrap(~ kmer) +
    theme_bw() +
    geom_point(aes(color=Method)) +
    geom_line(aes(color=Method)) +
    theme(axis.text.x=element_text(size=15, angle=45, hjust=1),
          axis.title.x=element_text(size=15),
          axis.text.y=element_text(size=15),
          axis.title.y=element_text(size=15),
          legend.text=element_text(size=15),
          legend.title=element_text(size=15),
          strip.text.x=element_text(size=15)) +
    ylab("Average F1 Score") +
    ggtitle("80X coverage, 0% error rate")

