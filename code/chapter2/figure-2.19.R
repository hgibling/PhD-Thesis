library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(scales)

HG003.raw.dat <- read.table("data/chapter2/HG003-raw.k-71.tsv", header=T)
HG003.spades.dat <- read.table("data/chapter2/HG003-spades.k-71.tsv", header=T)
HG003.graphaligner.dat <- read.table("data/chapter2/HG003-graphaligner.k-71.tsv", header=T)
HG003.aligncorrect.dat <- read.table("data/chapter2/HG003-aligncorrect.k-71.tsv", header=T)

HG003.raw <- HG003.raw.dat %>%
  mutate(Type=ifelse(Type=="ZnF", "Zinc finger", Type))
HG003.spades <- HG003.spades.dat %>%
  mutate(Type=ifelse(Type=="ZnF", "Zinc finger", Type))
HG003.graphaligner <- HG003.graphaligner.dat %>%
  mutate(Type=ifelse(Type=="ZnF", "Zinc finger", Type))
HG003.aligncorrect <- HG003.aligncorrect.dat %>%
  mutate(Type=ifelse(Type=="ZnF", "Zinc finger", Type))


raw <- ggplot(HG003.raw, aes(Position, ReadCount)) +
  theme_bw() +
  geom_bar(stat="identity", aes(fill=Type)) +
  scale_fill_manual(values=c(palette()[2:4])) +
  ggtitle("Raw reads") +
  guides(fill=guide_legend("Region")) +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylim(0,70) +
  ylab("Number of reads") +
  xlab("Position of *k*-mer")

spades <- ggplot(HG003.spades, aes(Position, ReadCount)) +
  theme_bw() +
  geom_bar(stat="identity", aes(fill=Type)) +
  scale_fill_manual(values=c(palette()[2:4])) +
  ggtitle("SPAdes corrected") +
  guides(fill=guide_legend("Region")) +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylim(0,70) +
  ylab("Number of reads") +
  xlab("Position of *k*-mer")

ga <- ggplot(HG003.graphaligner, aes(Position, ReadCount)) +
  theme_bw() +
  geom_bar(stat="identity", aes(fill=Type)) +
  scale_fill_manual(values=c(palette()[2:4])) +
  ggtitle("GraphAligner corrected") +
  guides(fill=guide_legend("Region")) +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylim(0,70) +
  ylab("Number of reads") +
  xlab("Position of *k*-mer")

ac <- ggplot(HG003.aligncorrect, aes(Position, ReadCount)) +
  theme_bw() +
  geom_bar(stat="identity", aes(fill=Type)) +
  scale_fill_manual(values=c(palette()[2:4])) +
  ggtitle("AlignCorrect corrected") +
  guides(fill=guide_legend("Region")) +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylim(0,70) +
  ylab("Number of reads") +
  xlab("Position of *k*-mer")

plot_grid(raw,spades,ga,ac, ncol=1, labels=LETTERS[1:4])