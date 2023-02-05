library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)

self.stats.data <- read.csv("data/chapter4/all-graphs-self-stats.csv", header=T) 

allele.names <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

self.stats <- self.stats.data %>%
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

self.stats$Allele <- factor(self.stats$Allele, levels=allele.names)
self.stats$Error <- factor(self.stats$Error, levels=c("0%", "0.1%", "1%"))
self.stats$Path <- factor(self.stats$Path, levels=c("yes", "no"))


self.stats.av <- self.stats %>%
  group_by(Aligner, Frag, Graph, Path, Coverage, Error) %>%
  summarize(MeanDiff=mean(MeanDiff), 
            MeanDiffNo0=mean(MeanDiffNo0),
            MeanDiffPos=mean(MeanDiffPos),
            MeanDiffNeg=mean(MeanDiffNeg),
            PropReadsDiff=mean(PropReadsDiff),
            PropPos=mean(PropPos),
            PropNeg=mean(PropNeg),
            NumReadsNoMapRef=mean(NumReadsNoMapRef),
            NumReadsNoMapGraph=mean(NumReadsNoMapGraph),
            PropReadsNoMapRef=mean(PropReadsNoMapRef),
            PropReadsNoMapGraph=mean(PropReadsNoMapGraph))

### Spurious w/ paths
# mg v ga v vg fixedfrag
g1.av <- self.stats.av %>%
  filter(Path=="yes", Error=="0%") %>%
  filter(!(Aligner=="vg" & Frag=="none")) %>%
  filter(!(Aligner=="vg" & Frag=="frag"))

g1.av.prop.diff.sub <- g1.av %>% 
  select(Aligner, Graph, Error, PropPos, PropNeg, PropReadsNoMapGraph) %>% 
  pivot_longer(cols=c(PropPos, PropNeg, PropReadsNoMapGraph), names_to="Divergence", values_to="Value") %>% 
  mutate(Divergence=case_when(
    Divergence=="PropPos" ~ "Improvement",
    Divergence=="PropNeg" ~ "Impairment",
    Divergence=="PropReadsNoMapGraph" ~ "Unmapped"
  ))
g1.av.prop.diff.sub$Divergence <- factor(g1.av.prop.diff.sub$Divergence,
                                         levels=c("Improvement", "Impairment", "Unmapped"))



### Spurious w/o paths
# mg v ga v vg fixedfrag
g2.av <- self.stats.av %>%
  filter(Path=="no", Error=="0%") %>%
  filter(!(Aligner=="vg" & Frag=="none")) %>%
  filter(!(Aligner=="vg" & Frag=="frag"))

g2.av.prop.diff.sub <- g2.av %>% 
  select(Aligner, Graph, Error, PropPos, PropNeg, PropReadsNoMapGraph) %>% 
  pivot_longer(cols=c(PropPos, PropNeg, PropReadsNoMapGraph), names_to="Divergence", values_to="Value") %>% 
  mutate(Divergence=case_when(
    Divergence=="PropPos" ~ "Improvement",
    Divergence=="PropNeg" ~ "Impairment",
    Divergence=="PropReadsNoMapGraph" ~ "Unmapped"
  ))
g2.av.prop.diff.sub$Divergence <- factor(g2.av.prop.diff.sub$Divergence,
                                         levels=c("Improvement", "Impairment",
                                                  "Unmapped"))



### compare path v no path
compare.mean.diff <- g2.av %>% 
  bind_rows(g1.av)

compare.mean.diff.sub <- compare.mean.diff %>% 
  select(Aligner, Graph, Error, MeanDiffPos, MeanDiffNeg) %>% 
  pivot_longer(cols=c(MeanDiffPos, MeanDiffNeg), 
               names_to="Divergence", values_to="Value") %>% 
  mutate(Divergence=case_when(
    Divergence=="MeanDiffPos" ~ "Improvement",
    Divergence=="MeanDiffNeg" ~ "Impairment"
  ))
compare.mean.diff.sub$Divergence  <- factor(compare.mean.diff.sub$Divergence,
                                            levels=c("Improvement", "Impairment"))
compare.mean.diff.sub$Path  <- factor(compare.mean.diff.sub$Path,
                                      levels=c("yes", "no"))


comp.mean <- ggplot(compare.mean.diff.sub, aes(Graph, Value)) +
  theme_bw() +
  facet_grid(Error ~ Aligner) +
  geom_hline(yintercept=0) + 
  geom_jitter(aes(shape=Divergence, fill=Path, color=Path),
              alpha=0.5, width=0.3, height=0) +
  scale_shape_manual(values=c(24,25),
                     name="Divergence difference") +
  scale_fill_manual(values=c(palette()[c(3,6)])) +
  scale_color_manual(values=c(palette()[c(3,6)])) +
  theme(axis.text.x=element_text(angle = 30, hjust = 1, vjust = 1)) +
  guides(shape = guide_legend(order = 2), 
         color = guide_legend(order = 1),
         fill=F) +
  xlab("Graph topology") +
  ylab("Average difference in divergence")


compare.prop.diff <- g2.av.prop.diff.sub %>%
  bind_rows(g1.av.prop.diff.sub)
compare.prop.diff$Path  <- factor(compare.prop.diff$Path,
                                  levels=c("yes", "no"))

comp.diff <- ggplot(compare.prop.diff, aes(Graph, Value)) +
  theme_bw() +
  facet_grid(Error ~ Aligner) +
  geom_jitter(aes(color=Path, shape=Divergence, fill=Path), 
              alpha=0.5, width=0.3, height=0) +
  scale_color_manual(values=c(palette()[c(3,6)])) +
  scale_fill_manual(values=c(palette()[c(3,6)])) +
  scale_shape_manual(values=c(24,25,4),
                     name="Read type") +
  theme(axis.text.x=element_text(angle = 30, hjust = 1, vjust = 1)) +
  guides(shape = guide_legend(order = 2), 
         color = guide_legend(order = 1),
         fill=F) +
  xlab("Graph topology") +
  ylab("Average proportion of reads")

plot_grid(comp.mean, NULL, comp.diff, ncol=1, align=T, axis='lr',
          labels=c("A", "", "B"), rel_heights=c(1,0.05,1))
