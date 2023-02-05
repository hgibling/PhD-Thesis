library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

corrected.stats.data <-read.csv("data/chapter4/all-corrected-stats.csv", header=T) 

allele.names <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

graph.stats <- graph.stats.data %>%
  mutate(PropReadsNoMapRef=NumReadsNoMapRef/NumReads,
         PropReadsNoMapGraph=NumReadsNoMapGraph/NumReads,
         Aligner=case_when(
           Aligner=="ga" ~ "GraphAligner",
           Aligner=="mg" ~ "minigraph",
           Aligner=="vg" ~ "vg"
         ),
         Error=case_when(
           Error==0 ~ "0%",
           Error==0.001 ~ "0.1%",
           Error==0.01 ~ "1%"
         ),
         Graph=case_when(
           Graph=="allele-msa-dag" ~ "allele-msa-acyclic",
           Graph=="znf-loop" ~ "zinc-finger-loop",
           Graph=="znf-msa" ~ "zinc-finger-msa",
           Graph=="stacked" ~ "allele-stacked",
           T ~ Graph
         ),
         Path=case_when(
           Path=="path" ~ "yes",
           Path=="nopath" ~ "no"
         )) %>% 
  mutate(across(matches("Num|Mean|Prop"), ~
                  ifelse(Aligner=="vg" & Graph=="allele-msa", NA, .)))
# alter diff scores for vg allele-msa to NA, since couldn't index

graph.stats$Allele <- factor(graph.stats$Allele, levels=allele.names)
graph.stats$Error <- factor(graph.stats$Error, levels=c("0%", "0.1%", "1%"))
graph.stats$Path <- factor(graph.stats$Path, levels=c("yes", "no"))

graph.stats.av <- graph.stats %>%
  group_by(Aligner, Frag, Graph, Path, Coverage, Error) %>%
  summarize(MeanDiff=mean(MeanDiff), 
            MeanDiffNo0=mean(MeanDiffNo0),
            PropReadsDiff=mean(PropReadsDiff), 
            NumReadsNoMapRef=mean(NumReadsNoMapRef),
            NumReadsNoMapGraph=mean(NumReadsNoMapGraph),
            PropReadsNoMapRef=mean(PropReadsNoMapRef),
            PropReadsNoMapGraph=mean(PropReadsNoMapGraph))


### Breakdown by aligner

graph.stats.av.vgfixed <- graph.stats.av %>%
  filter(!(Aligner=="vg" & Frag=="frag"), !(Aligner=="vg" & Frag=="none"))

diff <- ggplot(graph.stats.av.vgfixed, aes(Graph, MeanDiffNo0)) +
  theme_bw() +
  facet_grid(Error ~ Aligner) +
  geom_hline(yintercept=0) + 
  geom_jitter(aes(color=Path), alpha=0.5, width=0.3, height=0) + 
  scale_color_manual(values=c(palette()[c(3,6)])) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1)) +
  xlab("Graph topology") +
  ylab("Average difference in divergence")

graph.stats.prop <- graph.stats.av.vgfixed %>%
  pivot_longer(cols=c(PropReadsDiff, PropReadsNoMapGraph), 
               names_to="Category", values_to="Value") %>% 
  mutate(Category=case_when(
    Category=="PropReadsDiff" ~ "Divergence difference",
    Category=="PropReadsNoMapGraph" ~ "Unmapped"
  ))

prop <- ggplot(graph.stats.prop, aes(Graph, Value)) +
  theme_bw() +
  facet_grid(Error ~ Aligner) +
  geom_jitter(aes(color=Path, shape=Category), alpha=0.6, width=0.3, height=0) +
  scale_color_manual(values=c(palette()[c(3,6)])) +
  scale_shape_manual(values=c(16,4),
                     name="Read type") +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1)) +
  guides(color=guide_legend(order=1), 
         shape=guide_legend(order=2)) +
  xlab("Graph topology") +
  ylab("Average proportion of reads")

plot_grid(diff, NULL, prop, ncol=1, align=T, axis='lr',
          labels=c("A", "", "B"), rel_heights=c(1,0.05,1))
