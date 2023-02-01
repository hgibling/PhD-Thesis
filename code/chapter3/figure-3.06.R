library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)
library(cowplot)
library(ggtext)

hprc <- read.table("data/chapter3/gbkc-36-vs-106-HPRC.tsv", header=T, sep="\t")
ohs <- read.table("data/chapter3/gbkc-36-vs-106-OHS.tsv", header=T, sep="\t")

both <- hprc %>% 
  mutate(Dataset="HPRC++") %>% 
  bind_rows(ohs %>% mutate(Dataset="OHS")) %>% 
  rename(Comparison=Truth, Reads=Corrected) %>% 
  mutate(Comparison=ifelse(Comparison==T, "Truth genotypes", 
                           "Realignment model calls"),
         Reads=ifelse(Reads==T, "Corrected", "Raw"),
         NumSamples=case_when(
           Dataset=="HPRC++" & PRDM9List==36 & Comparison=="Truth genotypes" ~ 45,
           Dataset=="HPRC++" & PRDM9List==36 & Comparison=="Realignment model calls" ~ 52,
           Dataset=="OHS" & PRDM9List==36 & Comparison=="Truth genotypes" ~ 46,
           Dataset=="OHS" & PRDM9List==36 & Comparison=="Realignment model calls" ~ 52,
           Dataset=="HPRC++" & PRDM9List==106 ~ 52,
           Dataset=="OHS" & PRDM9List==106 ~ 49,
         )) %>% 
  mutate(TwoCorrectPerc=TwoCorrect/NumSamples,
         OneCorrectPerc=OneCorrect/NumSamples,
         Top10Perc=Top10/NumSamples)

both$PRDM9List <- factor(both$PRDM9List, levels=c(36, 106))
both$Model <- factor(both$Model, levels=c("count", "distance",  "combined"))
both$Comparison <- factor(both$Comparison, levels=c("Truth genotypes", 
                                                    "Realignment model calls"))

truth <- both %>% filter(Comparison=="Truth genotypes")

truth.long <- truth %>%
  pivot_longer(cols=c(TwoCorrectPerc, Top10Perc),
               names_to="Perc", values_to="Value")

truth <- truth %>%
  mutate(Reads=ifelse(Reads=="Corrected", "corrected", "raw"))

truth$Reads <- factor(truth$Reads, levels=c("raw", "corrected"))


hh <- ggplot(truth %>% filter(Dataset=="HPRC++"),
             aes(x=Method, group=PRDM9List)) + 
  theme_bw() +
  facet_grid(Reads ~ Model, scales="free_x", space="free_x") +
  geom_col(aes(y=Top10Perc,fill=PRDM9List), 
           position="dodge", alpha=0.3, color="black") +
  geom_col(aes(y=TwoCorrectPerc,fill=PRDM9List), 
           position="dodge", color="black") +
  scale_fill_manual(values=c(palette()[c(3:4)]),
                    labels=c("*PRDM9*-36","*PRDM9*-106")) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        plot.title = element_text(hjust = 0.5),
        legend.text=ggtext::element_markdown()) +
  guides(fill = guide_legend(override.aes = list(alpha = 0.3))) +
  labs(fill="Genotype in top 10") +
  xlab("") + 
  ylab("Proportion of samples called") +
  ggtitle("HPRC++") +
  ylim(0,1)

oo <- ggplot(truth %>% filter(Dataset=="OHS"),
             aes(x=Method, group=PRDM9List)) + 
  theme_bw() +
  facet_grid(Reads ~ Model, scales="free_x", space="free_x") +
  geom_col(aes(y=Top10Perc,fill=PRDM9List), 
           position="dodge", alpha=0.3, color="black") +
  geom_col(aes(y=TwoCorrectPerc,fill=PRDM9List), 
           position="dodge", color="black") +
  scale_fill_manual(values=c(palette()[c(3:4)]),
                    labels=c("*PRDM9*-36","*PRDM9*-106")) +
  theme(axis.text.x=element_text(angle=30, hjust=1),
        plot.title = element_text(hjust = 0.5),
        legend.text=ggtext::element_markdown()) +
  guides(fill = guide_legend(override.aes = list(alpha = 1))) +
  labs(fill="Genotype correct") +
  xlab("") + 
  ylab("Proportion of samples called") +
  ggtitle("OHS") +
  ylim(0,1)

plot_grid(hh, NULL, oo, ncol=1, align="v", axis="lr", 
          rel_heights=c(1,0.01,1), labels=c("A","","B"))