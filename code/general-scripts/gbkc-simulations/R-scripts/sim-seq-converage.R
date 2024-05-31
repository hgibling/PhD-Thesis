library(dplyr)
library(ggplot2)

setwd("~/prdm9/spotcheck/")

all.pos <- read.csv("all-pos.csv", stringsAsFactors=F, col.names=c("Iteration", "Error", "Coverage", "Position"))
all.pos$Iteration <- factor(all.pos$Iteration, levels=unique(all.pos$Iteration))
all.pos$Error <- factor(all.pos$Error, levels=unique(all.pos$Error))

ggplot(all.pos, aes(Iteration, Position)) +
    facet_grid(Error ~ .) +
    geom_violin(aes(color=Iteration))

ggplot(all.pos %>% filter(Error==0), aes(Coverage, Position)) +
    facet_wrap(~ Iteration) +
    geom_violin(aes(color=Iteration)) +
    geom_hline(yintercept=10000) +
    geom_hline(yintercept=11092)

ggplot(all.pos %>% filter(Error==0), aes(Position)) +
    geom_density(aes(color=Iteration)) +
    geom_vline(xintercept=10000) +
    geom_vline(xintercept=11092)

ggplot(all.pos %>% filter(Error==0), aes(Position)) +
    geom_histogram(binwidth=10, aes(fill=Iteration)) +
    coord_cartesian(xlim=c(9000, 12100)) +
    geom_vline(xintercept=10000) +
    geom_vline(xintercept=11092)

ggplot(all.pos %>% filter(Error==0), aes(Position)) +
    geom_histogram(binwidth=100, aes(fill=Iteration)) +
    geom_vline(xintercept=10000) +
    geom_vline(xintercept=11092)

ggplot(all.pos %>% filter(Error==0), aes(Position)) +
    facet_wrap(~ Iteration) +
    geom_histogram(binwidth=100, aes(fill=Iteration)) +
    geom_vline(xintercept=10000) +
    geom_vline(xintercept=11092)

ggplot(all.pos %>% filter(Error==0, Iteration==1), aes(Position)) +
    geom_histogram(binwidth=10, aes(fill=Iteration)) +
    geom_vline(xintercept=10000) +
    geom_vline(xintercept=11092)
