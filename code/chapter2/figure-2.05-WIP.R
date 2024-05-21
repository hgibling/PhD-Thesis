library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)


setwd("~/cluster/scratch/test-sims-fragments/new-sims")

reads.dat <- read.table("all-reads-average-F1-noflank.tsv", header=T)
frags.dat <- read.table("all-frags-average-F1-noflank.tsv", header=T)

reads.flank.dat <- read.table("all-reads-average-F1-flank.tsv", header=T)
frags.flank.dat <- read.table("all-frags-average-F1-flank.tsv", header=T)

reads <- reads.dat %>% 
  mutate(Flank="no flank") %>% 
  bind_rows(reads.flank.dat %>%
              mutate(Flank="10k flank"))

frags <- frags.dat %>% 
  mutate(Flank="no flank") %>% 
  bind_rows(frags.flank.dat %>%
              mutate(Flank="10k flank"))

reads$ReadLength <- factor(reads$ReadLength, levels=as.character(seq(50, 300, 25)))
reads$Flank <- factor(reads$Flank, levels=c("no flank", "10k flank"))
frags$FragmentLength <- factor(frags$FragmentLength, levels=as.character(seq(200, 500, 50)))
frags$Flank <- factor(frags$Flank, levels=c("no flank", "10k flank"))


rds <- ggplot(reads, aes(k, AverageF1)) +
  theme_bw() +
  facet_grid(~Flank) +
  geom_point(aes(color=ReadLength)) +
  geom_line(aes(color=ReadLength)) +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylab("Average F1 Score") +
  xlab("*k*") +
  guides(color=guide_legend(title="Read length")) 


frg <- ggplot(frags, aes(k, AverageF1)) +
  theme_bw() +
  facet_grid(~Flank) +
  geom_point(aes(color=FragmentLength)) +
  geom_line(aes(color=FragmentLength)) +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylab("Average F1 Score") +
  xlab("*k*") +
  guides(color=guide_legend(title="Fragment length")) 

plot_grid(rds, NULL, frg, ncol=1, align='v', rel_heights=c(1,0.1,1), 
          labels=c("A", "", "B"))  
