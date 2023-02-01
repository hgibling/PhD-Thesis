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

dplot <- ggplot(all.average %>% filter(Model=="distance"), aes(k, AverageF1)) +
  theme_bw() +
  facet_grid(Coverage ~ Error) +
  geom_point(aes(color=Method), alpha=0.5) +
  geom_line(aes(color=Method), alpha=0.5) +
  scale_color_manual(values=c(palette()[c(1, 3, 4, 6)])) +
  xlab("*k*") +
  ylab("Average F1 score") +
  theme(axis.title.x = ggtext::element_markdown())

dplot.zoom <- ggplot(all.average %>% filter(Model=="distance", Method=="max"), aes(k, AverageF1)) +
  theme_bw() +
  facet_grid(Method ~ Error) +
  geom_point(aes(color=Coverage), alpha=0.5) +
  geom_line(aes(color=Coverage), alpha=0.5) +
  scale_color_manual(values=c(palette()[c(1:4, 6)])) +
  xlab("*k*") +
  ylab("Average F1 score") +
  theme(axis.title.x = ggtext::element_markdown()) +
  coord_cartesian(ylim=c(0.75, 1))


plot_grid(dplot, dplot.zoom, labels=LETTERS[1:2], ncol=1, 
          align=T, axis="lr", rel_heights=c(2,1))