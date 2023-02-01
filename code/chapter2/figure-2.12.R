library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggtext)

all.average.data <- read.table("data/chapter2/all-average-hap-F1.tsv", header=T)

all.average <- all.average.data %>% 
  mutate(Model=gsub("s$", "", Model),
         Coverage=paste0(Coverage, "X"),
         Error=case_when(
           Error==0 ~ "0%",
           Error==0.001 ~ "0.1%",
           Error==0.01 ~ "1%"
         )) %>% 
  mutate(Coverage=factor(Coverage, levels=c("20X", "40X", "60X", "80X", "100X")),
         Error=factor(Error, levels=c("0%", "0.1%", "1%")))

## COMBINED

all.average.combined <- all.average %>% 
  filter(Model=="combined") %>% 
  mutate(CombinedMethod=Method) %>% 
  separate(Method, into=c("CountMethod", "DistanceMethod")) %>% 
  mutate(CountMethod=case_when(
    CountMethod=="mean" ~ "flank mean",
    CountMethod=="median" ~ "flank median",
    CountMethod=="coverage" ~ "coverage"
  ))

all.comb <- ggplot(all.average.combined %>% filter(Error=="0%", Coverage=="20X"), aes(k, AverageF1)) +
  theme_bw() +
  facet_grid(Coverage ~ Error) +
  geom_point(aes(color=DistanceMethod), alpha=0.5) +
  geom_line(aes(color=DistanceMethod, linetype=CountMethod), alpha=0.5) +
  scale_color_manual(values=c(palette()[c(1:4)])) +
  scale_linetype_manual(values=c(1:3)) +
  xlab("*k*") +
  ylab("Average F1 score") +
  labs(color="Distance method", linetype="Count method") +
  theme(axis.title.x = ggtext::element_markdown()) +
  coord_cartesian(ylim=c(0.9,1))

all.average.subset <- all.average %>% 
  filter(Model=="count", Method=="coverage") %>% 
  bind_rows(all.average %>% 
              filter(Model=="distance", Method=="max")) %>% 
  bind_rows(all.average %>% 
              filter(Model=="combined", Method=="coverage_max"))

best <- ggplot(all.average.subset, aes(k, AverageF1)) +
  theme_bw() +
  facet_grid(Coverage ~ Error) +
  geom_point(aes(color=Method), alpha=0.5) +
  geom_line(aes(color=Method), alpha=0.5) +
  scale_color_manual(values=c(palette()[c(1,2,4)]),
                     labels=c("count-coverage", 
                              "combined count-coverage \n& distance-max",
                              "distance-max")) +
  xlab("*k*") +
  ylab("Average F1 score") +
  theme(axis.title.x = ggtext::element_markdown()) +
  coord_cartesian(ylim=c(0.90,1))

plot_grid(all.comb, best, labels=LETTERS[1:2], ncol=1, 
          align=T, axis="lr", rel_heights=c(1,2))
