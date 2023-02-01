library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

graph.stats.data <- read.csv("data/chapter4/all-graphs-alignment-stats.csv", header=T) 

allele.names <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

graph.stats <- graph.stats.data %>%
  mutate(PropReadsNoMapRef=NumReadsNoMapGraph/NumReads,
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
         Frag=case_when(
           Frag=="fixedfrag" ~ "Fixed",
           Frag=="frag" ~ "Variable",
           Frag=="none" ~ "Not specified"
         ),
         Graph=case_when(
           Graph=="allele-msa-dag" ~ "allele-msa-acyclic",
           Graph=="znf-loop" ~ "zinc-finger-loop",
           Graph=="znf-msa" ~ "zinc-finger-msa",
           Graph=="stacked" ~ "allele-stacked",
           T ~ Graph
         )) %>% 
  mutate(across(matches("Num|Mean|Prop"), ~
                  ifelse(Aligner=="vg" & Graph=="allele-msa", NA, .)))
# alter diff scores for vg allele-msa to NA, since couldn't index
graph.stats$Allele <- factor(graph.stats$Allele, levels=allele.names)
graph.stats$Error <- factor(graph.stats$Error, levels=c("0%", "0.1%", "1%"))
graph.stats$Frag <- factor(graph.stats$Frag, levels=c("Fixed", "Variable", "Not specified"))

graph.stats.av <- graph.stats %>%
  group_by(Aligner, Frag, Graph, Path, Coverage, Error) %>%
  summarize(MeanDiff=mean(MeanDiff), 
            MeanDiffNo0=mean(MeanDiffNo0),
            PropReadsDiff=mean(PropReadsDiff), 
            NumReadsNoMapRef=mean(NumReadsNoMapRef),
            NumReadsNoMapGraph=mean(NumReadsNoMapGraph),
            PropReadsNoMapRef=mean(PropReadsNoMapRef),
            PropReadsNoMapGraph=mean(PropReadsNoMapGraph))

### Improvements
#  vg comparison (paths)
g3 <- graph.stats %>%
  filter(Path=="path", Aligner=="vg")
g3.av <- graph.stats.av %>%
  filter(Path=="path", Aligner=="vg")

vg.box.path <- ggplot(g3, aes(Graph, MeanDiffNo0)) +
  theme_bw() +
  facet_grid(Error ~ Frag, scales="free_y") +
  geom_hline(yintercept=0) +
  geom_boxplot(aes(fill=Graph), outlier.shape=NA, show.legend=T) + 
  geom_jitter(alpha=0.35, aes(shape=20), width=0.2, height=0) + 
  scale_shape_identity() +
  scale_fill_manual(values=c(palette()[c(3:4,6,7)])) + 
  theme(axis.text.x=element_text(angle = 30, hjust = 1, vjust = 1)) +
  xlab("Graph topology") +
  ylab("Average difference in divergence")


g3.av.stats.prop <- g3.av %>%
  pivot_longer(cols=c(PropReadsDiff, PropReadsNoMapGraph), 
               names_to="Category", values_to="Value") %>% 
  mutate(Category=case_when(
    Category=="PropReadsDiff" ~ "Divergence difference",
    Category=="PropReadsNoMapGraph" ~ "Unmapped"
  ))

vg.prop <- ggplot(g3.av.stats.prop, aes(Graph, Value)) +
  theme_bw() +
  facet_grid(Error ~ Frag) +
  geom_jitter(aes(color=Graph, shape=Category), 
              alpha=1, width=0.3, height=0) +
  scale_color_manual(values=c(palette()[c(2:4,6,7)])) +
  scale_shape_manual(values=c(16,4),
                     name="Read type") +
  theme(axis.text.x=element_text(angle = 30, hjust = 1, vjust = 1)) +
  guides(color=F) +
  xlab("Graph topology") +
  ylab("Average proportion of reads")

plot_grid(vg.box.path, NULL, vg.prop, ncol=1, align=T, axis='lr',
          labels=c("A", "", "B"), rel_heights=c(1,0.05,1))
