library(dplyr)
library(ggplot2)

setwd("~/prdm9")

compare <- read.csv("compare-profiles.csv", col.names=c("First", "Second", "kmer", "Difference", "k"), stringsAsFactors=F) 

compare.num.diff.kmers <- compare %>%
    group_by(First, Second, k) %>%
    count(Second) %>%
    arrange(n)

compare.abs.count.diff.kmers <- compare %>%
    group_by(First, Second, k) %>%
    summarize(Sum=sum(abs(Difference))) %>%
    arrange(Sum)

# duplicate data so each allele appears first and second
compare.equal.labels <- compare.abs.count.diff.kmers %>%
    rename(tempfirst=First, tempsecond=Second) %>%
    rename(First=tempsecond, Second=tempfirst) %>%
    bind_rows(compare.abs.count.diff.kmers) 

compare.means <- compare.equal.labels %>%
    group_by(First, k) %>%
    summarize(Mean=mean(Sum), StandardDev=sd(Sum)) %>%
    arrange(k, Mean)

ggplot(compare.means, aes(First, Mean)) +
    facet_wrap(~ k) +
    geom_point()

ggplot(compare.means, aes(k, Mean)) +
    geom_point(aes(color=First)) +
    geom_line(aes(color=First)) +
    geom_text(data=compare.means %>% filter(k==99), aes(label = First, colour = First), hjust = - 2) +
    scale_colour_discrete(guide='none')



####### PCA

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("~/prdm9")

profiles <- read.csv("all-alleles-flank10k-profiles-51k.csv", col.names=c("k", "Allele", "kmer", "Count"), 
                     stringsAsFactors=F) %>%
    spread(kmer, Count, fill=0)


row.names(profiles) <- profiles$Allele

profiles.mat <- profiles %>%
    select(-k, -Allele)

# remove kmers with identical counts across all alleles
keep.cols <- c()
for (i in 1:ncol(profiles.mat)) {
    if (length(unique(profiles.mat[,i])) != 1) {
        keep.cols <- c(keep.cols, i)
    }
}

profiles.mat.diff <- profiles.mat[,keep.cols]

profiles.pca <- prcomp(profiles.mat.diff)

# get znf content info
znf <- read.csv("allele-znf-content.csv", stringsAsFactors=F, col.names=c("Allele", "Znf"), header=F)

for (i in letters) {
    znf <- znf %>%
        mutate(!!i := as.character(nchar(gsub(paste0("[^", !!i, "]"), "", Znf))))
    if (length(unique(znf[,i])) == 1) {
        znf <- znf %>%
            select(-!!i)
    }
}


# combine into single dataframe

profiles.znf <- data.frame(profiles.pca$x) %>%
    add_rownames(var="Allele") %>%
    full_join(znf)

ggplot(profiles.znf, aes(PC1, PC2)) +
    theme_bw() +
    geom_point(aes(color=k), alpha=0.5) +
    geom_text_repel(label=rownames(profiles.mat.diff)) +
    xlab(paste0("PC1 (", summary(profiles.pca)$importance[2,1]*100, "%)")) +
    ylab(paste0("PC2 (", summary(profiles.pca)$importance[2,2]*100, "%)")) +
    ggtitle("PCA for k-mer count profiles at k=51")
    
    
ggplot(profiles.pca$x, aes(PC1, PC2)) +
    theme_bw() +
    geom_point(aes(alpha=0.5)) +
    geom_text_repel(label=rownames(profiles.mat.diff)) +
    xlab(paste0("PC1 (", summary(profiles.pca)$importance[2,1]*100, "%)")) +
    ylab(paste0("PC2 (", summary(profiles.pca)$importance[2,2]*100, "%)")) +
    ggtitle("PCA for k-mer count profiles at k=51")

profiles.pca.cs <- prcomp(profiles.mat.diff, center=T, scale=T)

ggplot(profiles.pca.cs$x, aes(PC1, PC2)) +
    theme_bw() +
    geom_point(aes(alpha=0.5)) +
    geom_text_repel(label=rownames(profiles.mat.diff)) +
    xlab(paste0("PC1 (", summary(profiles.pca)$importance[2,1]*100, "%)")) +
    ylab(paste0("PC2 (", summary(profiles.pca)$importance[2,2]*100, "%)")) +
    ggtitle("PCA for k-mer count profiles at k=51 (centred and scaled)")



####### PCA diploid

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)

setwd("~/prdm9")

profiles <- read.csv("all-genotypes-flank10k-profiles-51k.csv", col.names=c("k", "Genotype", "kmer", "Count"), 
                     stringsAsFactors=F) %>%
    spread(kmer, Count, fill=0)


row.names(profiles) <- profiles$Genotype

profiles.mat <- profiles %>%
    select(-k, -Genotype)

# remove kmers with identical counts across all alleles
keep.cols <- c()
for (i in 1:ncol(profiles.mat)) {
    if (length(unique(profiles.mat[,i])) != 1) {
        keep.cols <- c(keep.cols, i)
    }
}

profiles.mat.diff <- profiles.mat[,keep.cols]

profiles.pca <- prcomp(profiles.mat.diff)

ggplot(data.frame(profiles.pca$x), aes(PC1, PC2)) +
    theme_bw() +
    geom_point(alpha=0.5) +
    #geom_text_repel(label=rownames(profiles.mat.diff)) +
    xlab(paste0("PC1 (", summary(profiles.pca)$importance[2,1]*100, "%)")) +
    ylab(paste0("PC2 (", summary(profiles.pca)$importance[2,2]*100, "%)")) +
    ggtitle("PCA for k-mer count profiles at k=51")

pca3d(profiles.pca$x)


# get znf content info
znf <- read.csv("allele-znf-content.csv", stringsAsFactors=F, col.names=c("Allele", "Znf"), header=F)

for (i in letters) {
    znf <- znf %>%
        mutate(!!i := as.character(nchar(gsub(paste0("[^", !!i, "]"), "", Znf))))
    if (length(unique(znf[,i])) == 1) {
        znf <- znf %>%
            select(-!!i)
    }
}

geno <- data.frame(Genotype=rownames(profiles)) %>%
    separate(Genotype, c("Allele1", "Allele2"), remove=F) %>%
    gather("Position", "Allele", 2:3) %>%
    arrange(Genotype) %>%
    full_join(znf) %>%
    select(-c(Position, Allele, Znf)) %>%
    gather("Finger", "Count", -1) %>%
    mutate(Count=as.numeric(Count)) %>%
    group_by(Genotype, Finger) %>%
    summarize(Count=sum(Count)) %>%
    ungroup() %>%
    mutate(Count=as.character(Count)) %>%
    spread(Finger, Count) %>%
    arrange(Genotype)

# combine

profiles.znf <- data.frame(profiles.pca$x[,1:3]) %>%
    add_rownames(var="Genotype") %>%
    full_join(geno)

ggplot(profiles.znf, aes(PC1, PC2)) +
    theme_bw() +
    geom_point(aes(color=k), alpha=0.5) +
    #geom_text_repel(label=rownames(profiles.mat.diff)) +
    xlab(paste0("PC1 (", summary(profiles.pca)$importance[2,1]*100, "%)")) +
    ylab(paste0("PC2 (", summary(profiles.pca)$importance[2,2]*100, "%)")) +
    ggtitle("PCA for k-mer count profiles at k=51")

