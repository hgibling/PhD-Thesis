library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(scales)

dat <- read.csv("data/chapter2/ash-trio-all-count-results.csv",
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

ggplot(all.samples, aes(k, Rank)) +
  theme_bw() +
  facet_grid(Sample ~ Method) +
  geom_point(aes(color=Called, shape=Data), alpha=0.5) +
  geom_line(aes(linetype=Data), color=palette()[2], alpha=0.25) +
  scale_shape_manual(values=c(19,1)) +
  scale_color_manual(values=palette()[c(2)]) +
  xlab("*k*") +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylab("Rank of truth genotype") +
  guides(shape=guide_legend(override.aes=list(color=palette()[1])),
         color=guide_legend(title="Called correctly"))

