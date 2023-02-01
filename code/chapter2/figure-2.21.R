library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(scales)

dat <- read.csv("data/chapter2/ash-trio-count-results-corrected.csv",
                header=F, col.names=c("Sample", "Data", "Method",
                                      "k", "Tested", "Score"))
all.samples <- dat %>% 
  group_by(Sample, Data, Method, k) %>%
  mutate(Rank=rank(-Score, ties.method="min"),
         Genotype=case_when(
           Sample=="HG004" ~ "A/L37",
           T ~ "A/A"
         )) %>%
  mutate(Called=case_when(
    Rank==1 & Tested==Genotype ~ "True",
    T ~ "False"),
    Data=case_when(
      Data=="300x" ~ "2x150bp 300X",
      Data=="2x250" ~ "2x250bp 60X"),
    Method=case_when(
      Method!="coverage" ~ paste("flank", Method, sep=" "),
      T ~ Method
    )) %>%
  filter(Tested==Genotype)

all.samples$Called <- factor(all.samples$Called, levels=c("True", "False"))

short <- ggplot(all.samples %>% filter(Data=="2x150bp 300X"), aes(k, Rank)) +
  theme_bw() +
  facet_grid(Sample ~ Method) +
  geom_point(aes(color=Called, shape=Data), alpha=0.5) +
  geom_line(aes(linetype=Data), color=palette()[2], alpha=0.25) +
  scale_shape_manual(values=c(19)) +
  scale_color_manual(values=palette()[c(4,2)]) +
  xlab("*k*") +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylab("Rank of truth genotype") +
  guides(shape=guide_legend(override.aes=list(color=palette()[1])),
         color=guide_legend(title="Called correctly")) +
  coord_cartesian(ylim=c(0,200))

long <- ggplot(all.samples %>% filter(Data=="2x250bp 60X"), aes(k, Rank)) +
  theme_bw() +
  facet_grid(Sample ~ Method) +
  geom_point(aes(color=Called, shape=Data), alpha=0.5) +
  geom_line(aes(linetype=Data), color=palette()[2], alpha=0.25) +
  scale_shape_manual(values=c(1)) +
  scale_color_manual(values=palette()[c(4,2)]) +
  xlab("*k*") +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylab("Rank of truth genotype") +
  guides(shape=guide_legend(override.aes=list(color=palette()[1])),
         color=guide_legend(title="Called correctly")) +
  coord_cartesian(ylim=c(0,200))

plot_grid(short, NULL, long, ncol=1, labels=c("A", "", "B"),
          align="v", axis="lr", rel_heights=c(1,0.1,1))
