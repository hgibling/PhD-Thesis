library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(scales)

cortex.dat <- read.table("data/chapter2/average-cortex-F1-100x.tsv", header=T)
cortex.20x.dat <- read.table("data/chapter2/average-cortex-F1-20x.tsv", header=T)

pearson.dat <- read.table("data/chapter2/average-pearson-F1-100x.tsv", header=T)
pearson.20x.dat <- read.table("data/chapter2/average-pearson-F1-20x.tsv", header=T)

gbkc.dat <- read.table("data/chapter2/average-gblr-F1-100x.tsv", header=T)
gbkc.20x.dat <- read.table("data/chapter2/average-gblr-F1-20x.tsv", header=T)

both <- gbkc.dat %>% 
  filter(Method=="coverage") %>%
  select(Error, k, Coverage, AverageF1) %>%
  mutate(Scoring="Modified Cortex") %>%
  bind_rows(gbkc.20x.dat %>% 
              filter(Method=="coverage") %>%
              select(Error, k, Coverage, AverageF1) %>%
              mutate(Scoring="Modified Cortex")) %>%
  bind_rows(cortex.dat %>%
              mutate(Scoring="Cortex", Coverage=100)) %>%
  bind_rows(cortex.20x.dat %>%
              mutate(Scoring="Cortex", Coverage=20)) %>%
  bind_rows(pearson.dat %>%
              mutate(Scoring="Pearson", Coverage=100)) %>%
  bind_rows(pearson.20x.dat %>%
              mutate(Scoring="Pearson", Coverage=20)) %>%
  mutate(Error=case_when(
    Error==0 ~ "0%",
    Error==0.001 ~ "0.1%",
    Error==0.01 ~ "1%"),
    Coverage=paste0(as.character(Coverage), "X"))

both$Error <- factor(both$Error, levels=c("0%", "0.1%", "1%"))
both$Scoring <- factor(both$Scoring, levels=c("Cortex", "Modified Cortex", "Pearson"))
both$Coverage <- factor(both$Coverage, levels=c("20X", "100X"))


cortex.plot <- ggplot(both, aes(k, AverageF1)) +
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




hmm.data <- read.csv("data/chapter2/all-HMM-scores-20x-noflank.csv", 
                     header=F,
                     col.names=c("Simulated", "Tested", "Score", "Iteration"))
gbkc.data <- read.csv("data/chapter2/all-gbkc-scores-20x-noflank.csv", 
                      header=F,
                      col.names=c("Simulated", "k", "Tested", "Score", "Iteration"))

alleles.list <- c(LETTERS[1:5], paste0("L", seq(1,24)), paste0("L", seq(32,38)))

hmm.iteration <- hmm.data %>% 
  group_by(Simulated, Iteration) %>%
  mutate(Rank=rank(-Score, ties.method="min")) %>%
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

hmm.alleles <- hmm.iteration %>% 
  group_by(Simulated) %>%
  summarize(AlleleF1=mean(IterF1)) %>%
  ungroup()

hmm.average <- hmm.alleles %>% 
  summarize(AverageF1=mean(AlleleF1))


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

gbkc.alleles <- gbkc.iteration %>% 
  group_by(Simulated, k) %>%
  summarize(AlleleF1=mean(IterF1)) %>%
  ungroup()

gbkc.average <- gbkc.alleles %>% 
  group_by(k) %>%
  summarize(AverageF1=mean(AlleleF1)) %>%
  mutate(Type="All alleles")


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

hmm.alleles$Simulated <- factor(hmm.alleles$Simulated, levels=alleles.list)
gbkc.alleles$Simulated <- factor(gbkc.alleles$Simulated, levels=alleles.list)

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

bottom <- plot_grid(hmm.average.plot, NULL, hmm.allele.plot,
                    ncol=3, rel_widths=c(1,0.05,3.1),
                    labels=c("B","","C"))

plot_grid(cortex.plot, NULL, bottom,
          nrow=3, rel_heights=c(1,0.05,1.5), labels=c("A","",""))
