library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(scales)
library(stringr)

HG003.raw.dat <- read.table("data/chapter2/HG003-raw.k-71.tsv", header=T)

trio.kmercounts.dat <- read.csv("data/chapter2/ash-trio-illumina-raw-noncanon.k-71.jf.flank.csv", 
                            stringsAsFactors=F, header=F,
                            col.names=c("Sample", "kmer", "Count")) 

gc.bins <- read.table("data/chapter2/B-flank10k.k-71-GC.tsv", header=T)

HG003.raw <- HG003.raw.dat %>%
  mutate(Type=ifelse(Type=="ZnF", "Zinc finger", Type))

bins.pos <- gc.bins %>%
  mutate(Bin=Bin*250-250)
num.bins.pos <- nrow(bins.pos) -1

bin.dat <- bins.pos[1:num.bins.pos,]
maxval <- max(HG003.raw$ReadCount)
# 61

gc.hg003 <- ggplot(HG003.raw %>% mutate(Sample="HG003"), aes(Position, ReadCount)) +
  theme_bw() +
  facet_grid(~Sample) +
  scale_y_continuous(sec.axis=sec_axis(trans=~./maxval, name="*k*-mer GC%",
                                       label=percent_format(scale=100))) +
  scale_fill_manual(values=c(palette()[2:4])) +
  coord_cartesian(ylim=c(0, maxval)) +
  geom_bar(stat="identity", aes(fill=Type)) +
  geom_point(data=bin.dat, aes(x=Bin, y=GC*maxval), alpha=0.5) +
  geom_line(data=bin.dat, aes(x=Bin, y=GC*maxval), alpha=0.5) +
  guides(fill=guide_legend("Region")) +
  theme(axis.title.x=ggtext::element_markdown(),
        axis.title.y.right=ggtext::element_markdown()) +
  ylab("Number of reads") +
  xlab("Position of *k*-mer")


 trio.kmercounts <- trio.kmercounts.dat %>%
  rowwise() %>%
  mutate(PercentGC=(str_count(kmer, "G|C")/71) * 100) %>%
  select(-kmer) 

trio.10 <- trio.kmercounts %>%
  mutate(GCBin=case_when(
    PercentGC < 10 ~ "0 - 10",
    PercentGC >= 10 & PercentGC < 20 ~ "10 - 20",
    PercentGC >=20 & PercentGC < 30 ~ "20 - 30",
    PercentGC >=30 & PercentGC < 40 ~ "30 - 40",
    PercentGC >=40 & PercentGC < 50 ~ "40 - 50",
    PercentGC >=50 & PercentGC < 60 ~ "50 - 60",
    PercentGC >=60 & PercentGC < 70 ~ "60 - 70",
    PercentGC >=70 & PercentGC < 80 ~ "70 - 80",
    PercentGC >=80 & PercentGC < 90 ~ "80 - 90",
    PercentGC >= 100 ~ "90 - 100")) %>%
  select(-PercentGC) %>%
  group_by(Sample, Count, GCBin) %>%
  summarize(NumKmer=n())

trio10plot <- ggplot(trio.10,  aes(Count, NumKmer)) +
  theme_bw() +
  facet_grid(Sample ~ .) +
  geom_bar(stat="identity", aes(fill=GCBin)) +
  scale_fill_manual(values=c(palette()[2:9])) +
  xlab("*k*-mer count") +
  ylab("Number of distinct *k*-mers") +
  guides(fill=guide_legend("GC% bin")) +
  theme(axis.title.x=ggtext::element_markdown(),
        axis.title.y=ggtext::element_markdown())

plot_grid(gc.hg003, NULL, trio10plot, ncol=1, labels=c("A", "", "B"),
          align="v", axis="lr", rel_heights=c(1,0.1,2))
