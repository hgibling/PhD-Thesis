### plot for figure 2.7
### comparison of F1 scores for modified cortex (gbkc), cortex, pearson, and HMM methods

# load libaries
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)

# load data files for modified cortex/gbkc, cortex, and pearson 20X and 100X
cortex.100x.data <- read.table("data/chapter2/average-cortex-F1-100x.tsv", header=T)
cortex.20x.data <- read.table("data/chapter2/average-cortex-F1-20x.tsv", header=T)

pearson.100x.data <- read.table("data/chapter2/average-pearson-F1-100x.tsv", header=T)
pearson.20x.data <- read.table("data/chapter2/average-pearson-F1-20x.tsv", header=T)

gbkc.100x.data <- read.table("data/chapter2/average-gbkc-F1-100x.tsv", header=T)
gbkc.20x.data <- read.table("data/chapter2/average-gbkc-F1-20x.tsv", header=T)

# combine into single data set and give names for scoring methods
all.methods <- gbkc.100x.data %>% 
  # subset to count-coverage model and remove extra columns
  filter(Method=="coverage") %>%
  select(Error, k, Coverage, AverageF1) %>%
  mutate(Scoring="Modified Cortex") %>%
  bind_rows(gbkc.20x.data %>% 
              filter(Method=="coverage") %>%
              select(Error, k, Coverage, AverageF1) %>%
              mutate(Scoring="Modified Cortex")) %>%
  bind_rows(cortex.100x.data %>%
              mutate(Scoring="Cortex", Coverage=100)) %>%
  bind_rows(cortex.20x.data %>%
              mutate(Scoring="Cortex", Coverage=20)) %>%
  bind_rows(pearson.100x.data %>%
              mutate(Scoring="Pearson", Coverage=100)) %>%
  bind_rows(pearson.20x.data %>%
              mutate(Scoring="Pearson", Coverage=20)) %>%
  # clean up values for error and coverage
  mutate(Error=case_when(
    Error==0 ~ "0%",
    Error==0.001 ~ "0.1%",
    Error==0.01 ~ "1%"),
    Coverage=paste0(as.character(Coverage), "X"))

# order factor levels
all.methods$Error <- factor(all.methods$Error, levels=c("0%", "0.1%", "1%"))
all.methods$Scoring <- factor(all.methods$Scoring, levels=c("Cortex", "Modified Cortex", "Pearson"))
all.methods$Coverage <- factor(all.methods$Coverage, levels=c("20X", "100X"))

# plot comparisons
kmer.methods.plot <- ggplot(all.methods, aes(k, AverageF1)) +
  theme_bw() +
  facet_grid(Coverage~Error) +
  geom_point(aes(color=Scoring), alpha=0.5) +
  geom_line(aes(color=Scoring), alpha=0.5) +
  theme(axis.title.x=ggtext::element_markdown()) +
  scale_color_manual(values=c(palette()[c(4,2,7)])) +
  guides(color=guide_legend(title="Method")) +
  ylab("Average F1 score") +
  xlab("*k*") +
  coord_cartesian(ylim=c(0.8,1))

# load data files for HMM and modified cortex/gbkc results (20X error-freee; did not use simulations with 10kb flanks)
hmm.data <- read.csv("data/chapter2/all-HMM-scores-20x-noflank.csv", 
                     header=F,
                     col.names=c("Simulated", "Tested", "Score", "Iteration"))
gbkc.data <- read.csv("data/chapter2/all-gbkc-scores-20x-noflank.csv", 
                      header=F,
                      col.names=c("Simulated", "k", "Tested", "Score", "Iteration"))

# define allele names in order
alleles.list <- c(LETTERS[1:5], paste0("L", seq(1,24)), paste0("L", seq(32,38)))

# calculate per-iteration F1 scores for HMM data
hmm.iteration <- hmm.data %>% 
  group_by(Simulated, Iteration) %>%
  mutate(Rank=rank(-Score, ties.method="min")) %>%
  # define true and false positives and negatives
  mutate(Type=case_when(
    Simulated==Tested & Rank==1 ~ "TP",
    Simulated==Tested & Rank!=1 ~ "FN",
    Simulated!=Tested & Rank==1 ~ "FP",
    Simulated!=Tested & Rank!=1 ~ "TN"
  )) %>% 
  mutate(Type=factor(Type, levels=c("TP", "FN", "FP", "TN"))) %>% 
  ungroup() %>%
  group_by(Simulated, Iteration, Type, .drop=F) %>% 
  summarize(TypeCount=n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from=Type, values_from="TypeCount") %>% 
  group_by(Simulated, Iteration) %>% 
  mutate(IterPrecision=TP/(TP+FP),
         IterRecall=TP/(TP+FN)) %>% 
  replace_na(list(IterPrecision=0, IterRecall=0)) %>% 
  mutate(IterF1=2*((IterPrecision*IterRecall)/(IterPrecision+IterRecall))) %>%
  replace_na(list(IterF1=0)) %>%
  ungroup()

# calculate per-allele F1 scores for HMM data
hmm.alleles <- hmm.iteration %>% 
  group_by(Simulated) %>%
  summarize(AlleleF1=mean(IterF1)) %>%
  ungroup()

# calculate average F1 score for HMM data
hmm.average <- hmm.alleles %>% 
  summarize(AverageF1=mean(AlleleF1))

# calculate per-iteration F1 scores for gbkc data
gbkc.iteration <- gbkc.data %>% 
  group_by(Simulated, k, Iteration) %>%
  mutate(Rank=rank(-Score, ties.method="min")) %>%
  mutate(Type=case_when(
    Simulated==Tested & Rank==1 ~ "TP",
    Simulated==Tested & Rank!=1 ~ "FN",
    Simulated!=Tested & Rank==1 ~ "FP",
    Simulated!=Tested & Rank!=1 ~ "TN"
  )) %>% 
  mutate(Type=factor(Type, levels=c("TP", "FN", "FP", "TN"))) %>% 
  ungroup() %>%
  group_by(Simulated, k, Iteration, Type, .drop=F) %>% 
  summarize(TypeCount=n()) %>% 
  ungroup() %>% 
  pivot_wider(names_from=Type, values_from="TypeCount") %>% 
  group_by(Simulated, k, Iteration) %>% 
  mutate(IterPrecision=TP/(TP+FP),
         IterRecall=TP/(TP+FN)) %>% 
  replace_na(list(IterPrecision=0, IterRecall=0)) %>% 
  mutate(IterF1=2*((IterPrecision*IterRecall)/(IterPrecision+IterRecall))) %>%
  replace_na(list(IterF1=0)) %>%
  ungroup()

# calculate per-allele F1 scores for gbkc data
gbkc.alleles <- gbkc.iteration %>% 
  group_by(Simulated, k) %>%
  summarize(AlleleF1=mean(IterF1)) %>%
  ungroup()

# calculate average F1 scores for gbkc data
gbkc.average <- gbkc.alleles %>% 
  group_by(k) %>%
  summarize(AverageF1=mean(AlleleF1)) %>%
  mutate(Type="All alleles")

# order factor levels for allele names
hmm.alleles$Simulated <- factor(hmm.alleles$Simulated, levels=alleles.list)
gbkc.alleles$Simulated <- factor(gbkc.alleles$Simulated, levels=alleles.list)

# plot average F1 scores
hmm.average.plot <- ggplot(gbkc.average, aes(k, AverageF1)) +
  theme_bw() +
  facet_wrap(~Type) +
  geom_point(alpha=0.5, color=palette()[2]) +
  geom_line(alpha=0, aes(linetype="Modified Cortex"), color=palette()[2]) +
  geom_line(alpha=0.5, color=palette()[2]) +
  geom_hline(data=hmm.average, aes(yintercept=AverageF1, linetype="HMM"), 
             alpha=0.75, size=1) +
  ylab("Average F1 score") +
  xlab("*k*") +
  theme(axis.title.x=ggtext::element_markdown()) +
  ylim(0,1) +
  theme(legend.position="none")

# plot per-allele F1 scores
hmm.allele.plot <- ggplot(gbkc.alleles, aes(k, AlleleF1)) +
  theme_bw() +
  facet_wrap(~ Simulated) +
  geom_point(alpha=0.5, color=palette()[2]) +
  geom_line(alpha=0, aes(linetype="Modified Cortex"), color=palette()[2]) +
  geom_line(alpha=0.5, color=palette()[2]) +
  geom_hline(data=hmm.alleles, aes(yintercept=AlleleF1, linetype="HMM"),  
             alpha=0.75, size=1) +
  ylab("Average F1 score") +
  xlab("*k*") +
  scale_x_continuous(breaks=c(0,50,100)) +
  scale_y_continuous(breaks=c(0,0.5,1)) +
  theme(axis.title.x=ggtext::element_markdown()) +
  guides(linetype=guide_legend(title="Method", 
                               override.aes=list(color=palette()[c(1,2)], 
                                                 linetype=c(1,1))))

# combine HMM vs gbkc plots
bottom <- plot_grid(hmm.average.plot, NULL, hmm.allele.plot,
                    ncol=3, rel_widths=c(1,0.05,3.1),
                    labels=c("B","","C"))

# combine all 3 plots
plot_grid(kmer.methods.plot, NULL, bottom,
          nrow=3, rel_heights=c(1,0.05,1.5), labels=c("A","",""))
