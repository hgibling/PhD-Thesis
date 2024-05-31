library(dplyr)
library(tidyr)
library(ggplot2)
library(lattice)

setwd("~/Documents/PRDM9-Allele-Calling")

all.probs.N <- read.csv("all-perfect-updated-probs-cov40.csv", header=F, 
                      col.names=c("Simulated", "Tested", "Probability"))

all.probs.N <- group_by(all.probs.N, Simulated) %>%
    arrange(desc(Probability)) %>%
    arrange(Simulated) %>%
    mutate(Correct=ifelse(Simulated==Tested, T, F)) %>%
    group_by(Simulated) %>%
    mutate(Unique=ifelse(Probability %in% Probability[duplicated(Probability)], F, T))

all.probs <- filter(all.probs.N, Simulated!="N", Tested!="N")

order <-c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))



# test

ggplot(filter(all.probs.N, Simulated=="N"), 
       aes(reorder(Tested, -Probability), Probability)) +
    geom_point(size=5, aes(color=Correct, shape=Unique)) +
    scale_shape_manual(values=c(15,16)) +
    theme_bw() +
    xlab("Allele Tested") +
    ylab("Log Probability")


# plot all

for (i in unique(all.probs$Simulated)) {
    ggplot(filter(all.probs, Simulated==i), 
       aes(reorder(Tested, -Probability), Probability)) +
    geom_point(size=5, aes(color=Correct, shape=Unique)) +
    scale_shape_manual(values=c(15, 16)) +
    theme_bw() +
    xlab("Allele Tested") +
    ylab("Log Probability")
    ggsave(paste0("plots/perfect-updated/cov20/", i, "-probs.pdf"), width=15, height=5)
}


all.probs.max <- filter(all.probs, Probability==max(Probability))
# duplicates because some alleles had more than one tested allele with 
    # highest probability

all.probs.max.true <- filter(all.probs.max, Correct==T)



# heatmap

all.probs.wide <- all.probs %>%
    select(Simulated, Tested, Probability) %>%
    spread(Tested, Probability) %>%
    arrange(factor(Simulated, levels=rev(order))) %>%
    ungroup() %>%
    select(order)
rownames(all.probs.wide) <- rev(colnames(all.probs.wide))

levelplot(t(as.matrix(all.probs.wide)), scale=list(x=list(rot=45)),
          col.regions=colorRampPalette(c("#004949","#009292","#ff6db6"))(1e4),
          aspect="fill", 
          xlab='Tested Allele \n', 
          ylab='Simulated Allele',
          colorkey=list(space='bottom'),
          main="Scores of tested alleles (log probability)"
          )

probs.binary <- t(apply(all.probs.wide, 1, function(z){
    1 * (z==max(z))
}))

levelplot(t(as.matrix(probs.binary)), scale=list(x=list(rot=45)),
          col.regions=colorRampPalette(c("white", "black"))(2),
          aspect="fill", 
          xlab="Tested Allele \n",
          ylab="Simulated Allele",
          colorkey=list(space='bottom'),
          main="Tested alleles with highest score"
)
