library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)

all.alleles.data <- read.table("data/chapter2/all-allele-hap-F1.tsv", header=T)

alleles <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

all.alleles <- all.alleles.data %>% 
  mutate(Model=gsub("s$", "", Model),
         Coverage=paste0(Coverage, "X"),
         Error=case_when(
           Error==0 ~ "0%",
           Error==0.001 ~ "0.1%",
           Error==0.01 ~ "1%"
         )) %>% 
  mutate(Coverage=factor(Coverage, levels=c("20X", "40X", "60X", "80X", "100X")),
         Error=factor(Error, levels=c("0%", "0.1%", "1%")),
         Simulated=factor(Simulated, levels=alleles))


### HAPLOID COUNT COVERAGE INDIVIDUAL ALLELE

ggplot(all.alleles %>% filter(Model=="count", Method=="coverage", Error=="0%"), aes(k, AlleleF1)) +
  theme_bw() +
  facet_wrap(~Simulated) +
  geom_point(aes(color=Coverage), alpha=0.5) +
  geom_line(aes(color=Coverage), alpha=0.5) +
  scale_color_manual(values=c(palette()[c(1:4,6)])) +
  xlab("*k*") +
  ylab("Average F1 score") +
  theme(axis.title.x = ggtext::element_markdown()) +
  coord_cartesian(ylim=c(0.75,1))
