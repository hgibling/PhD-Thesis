library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)

ill.dat <- read.table("data/chapter2/HG003-illumina-jcount.tsv", header=T) %>%
  mutate(type=case_when(
    type=="znf" ~ "zinc finger",
    type=="left" ~ "left flank",
    type=="right" ~ "right flank"
  ))
pac.dat <- read.table("data/chapter2/HG003-pacbio-jcount.tsv", header=T) %>%
  mutate(type=case_when(
    type=="znf" ~ "zinc finger",
    type=="left" ~ "left flank",
    type=="right" ~ "right flank"
  ))

ill.dat$type <- factor(ill.dat$type, levels=c("left flank", "zinc finger", "right flank"))
pac.dat$type <- factor(ill.dat$type, levels=c("left flank", "zinc finger", "right flank"))



pacplot <- ggplot(pac.dat, aes(position, (exact_kmer_count/base_depth))) +
  theme_bw() +
  geom_line(aes(color=type)) +
  scale_color_manual(values=c(palette()[2:4])) +
  ggtitle("PacBio HiFi") +
  guides(color=guide_legend("Region")) +
  theme(axis.title.x=ggtext::element_markdown(),
        axis.title.y=ggtext::element_markdown()) +
  ylab("*k*-mer count : per-base coverage ratio") +
  xlab("Position of *k*-mer (GRCh38)")

illplot <- ggplot(ill.dat, aes(position, (exact_kmer_count/base_depth))) +
  theme_bw() +
  geom_line(aes(color=type)) +
  scale_color_manual(values=c(palette()[2:4])) +
  ggtitle("Illumina 2x250bp") +
  guides(color=guide_legend("Region")) +
  theme(axis.title.x=ggtext::element_markdown(),
        axis.title.y=ggtext::element_markdown()) +
  ylab("*k*-mer count : per-base coverage ratio") +
  xlab("Position of *k*-mer (GRCh38)")

plot_grid(illplot, NULL, pacplot, ncol=1, labels=c("A","",  "B"), rel_heights=c(1,0.1,1))
