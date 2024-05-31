library(dplyr)

setwd("~/awadalla/PRDM9-Project/HMM/kmers/300/")

all.300 <- read.csv("all-alleles-300-k.kmercounts", header=F, 
                    col.names=c("Allele", "kmer", "Count")) %>%
    spread(Allele, Count, fill=0)

compare.profiles <- function(df) {
    pairs <- combn(36, 2, simplify=F)
    res <- rep(NA, length(pairs))
    for (i in 1:length(pairs)) {
        res[i] <- identical(df[,pairs[[i]][1]], df[,pairs[[i]][2]])
    }
    same <- pairs[which(res==TRUE)]
    return(length(same))
}

s <- seq(1, 99, 2)
compare <- data.frame(kmer=s, same=NA)
for (i in 1:50) {
    compare[i, 2] <- compare.profiles(profiles, (i*2)-1)
}

compare[which(compare$same>0),]