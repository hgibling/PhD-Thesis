library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(ggtext)
library(scales)

F1.data <- read.table("data/chapter4/all-graphs-kmer-F1", header=T)

F1 <- F1.data %>%
  filter(!is.na(graph)) %>%
  mutate_all(~replace(., is.na(.), 0)) %>%
  mutate(graph=gsub("-no", "", graph)) %>%
  mutate(graph=case_when(
    graph == "allele-msa-dag" ~ "allele-msa-acyclic",
    graph == "stacked" ~ "allele-stacked",
    graph == "znf-loop" ~ "zinc-finger-loop",
    graph == "znf-msa" ~ "zinc-finger-msa",
    T ~ graph))

F1$graph <- factor(F1$graph, levels=c("GRCh38", "allele-msa", 
                                      "allele-msa-acyclic", "allele-stacked", 
                                      "zinc-finger-loop", "zinc-finger-msa"))
nums <- F1 %>%
  select(graph, k, graph.num) %>%
  bind_rows(F1 %>%
              filter(graph=="allele-stacked") %>%
              select(k, allele.num) %>%
              mutate(graph="all alleles", .before=1) %>%
              rename(graph.num=allele.num)) 

nums$graph <- factor(nums$graph, levels=c("GRCh38", "allele-msa", 
                                          "allele-msa-acyclic", "allele-stacked", 
                                          "zinc-finger-loop", "zinc-finger-msa",
                                          "all alleles"))

nums.sum <- nums %>%
  group_by(graph) %>%
  summarize(tot=sum(graph.num))

ga2 <- ggplot(F1, aes(k, F1)) +
  theme_bw() + 
  geom_line(aes(color=graph, linetype=graph), size=1.5) +
  scale_color_manual(values=c(palette()[c(1,5,4,6,7,3)])) +
  scale_linetype_manual(values=c("longdash", "dashed", "twodash",
                                 "solid", "dotted", "dotdash")) +
  theme(axis.title.x=ggtext::element_markdown(),
        legend.position="none") +
  labs(x="*k*", y="Pseudo-F1 score")

gb <- ggplot(F1, aes(precision, recall)) +
  theme_bw() + 
  facet_wrap(~graph) + 
  geom_point(aes(color=k), size=2) +
  scale_colour_gradient(low = "#D0ECE7",high ="#0B5345") +
  theme(legend.title=ggtext::element_markdown(),
        axis.text.x = element_text(angle=45, hjust=1),
        plot.margin = margin(l=0, t=0, b=0, r=70, 'pt')) +
  labs(x="Precision", y="Recall", color="*k*")

gc <- ggplot(nums, aes(k, graph.num)) +
  theme_bw() +
  geom_line(aes(color=graph, linetype=graph), size=1.5) +
  scale_y_log10(labels=comma) +
  scale_color_manual(values=c(palette()[c(8,1,5,4,6,7,3)]),
                     limits=c("all alleles", "GRCh38", "allele-msa", 
                              "allele-msa-acyclic", "allele-stacked", 
                              "zinc-finger-loop", "zinc-finger-msa")) +
  scale_linetype_manual(values=c("dashed", "longdash", "dashed", "twodash",
                                 "solid", "dotted", "dotdash"),
                        limits=c("all alleles", "GRCh38", "allele-msa", 
                                 "allele-msa-acyclic", "allele-stacked", 
                                 "zinc-finger-loop", "zinc-finger-msa")) +
  #scale_y_continuous(labels = label_comma()) +
  theme(axis.title.x=ggtext::element_markdown(),
        axis.title.y=ggtext::element_markdown()) +
  guides(linetype=F) +
  labs(x="*k*", y="Number of distinct *k*-mers", color="Reference")

top <- plot_grid(ga2,NULL,gc, nrow=1, rel_widths=c(1,0.1,1.5), labels=c("A","", "B"))
plot_grid(top, NULL, gb, nrow=3, rel_heights=c(1,0.1,1), labels=c("","", "C"))