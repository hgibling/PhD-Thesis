library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)
library(cowplot)
library(ggtext)


### OHS ###
ohs.gbkc.data <- read.table("data/chapter3/gbkc-pacbio-OHS.tsv", header=F, col.names=c("PacBio.ID", "Method", "k", "gbkcGenotype", "Score"))

ohs.truth.genos.data <- read.table("data/chapter3/truth-genos-OHS.tsv",
                                   header=T)


ohs.truth.genos <- ohs.truth.genos.data %>% 
  rename(PacBio.ID=Sample) %>% 
  mutate(BaseGenotype=gsub("Novel_Similar_|Potential_False_", 
                           "", TruthGenotype)) %>% 
  separate(BaseGenotype, c("TrueFirst", "TrueSecond"), sep="/") %>% 
  mutate(Type=ifelse(TrueFirst==TrueSecond, "Hom", "Het")) 


ohs.calls <- ohs.gbkc.data %>% 
  group_by(PacBio.ID, Method, k) %>%
  mutate(Rank=rank(-Score, ties.method="min")) %>% 
  separate(gbkcGenotype, into=c("gbkcFirst", "gbkcSecond"), sep="/", remove=F) %>% 
  filter(PacBio.ID!="AWA4634-806") %>% 
  full_join(ohs.truth.genos) %>% 
  rowwise() %>% 
  mutate(NumCorrectAlleles=case_when(
    Type=="Hom" & (gbkcFirst==gbkcSecond & gbkcFirst==TrueFirst) ~ 2,
    Type=="Het" & (grepl(gbkcFirst, TruthGenotype) & 
                     grepl(gbkcSecond, TruthGenotype)) ~ 2,
    Type=="Het" & (grepl(gbkcFirst, TruthGenotype) | 
                     grepl(gbkcSecond, TruthGenotype)) ~ 1,
    Type=="Hom" & (gbkcFirst==TrueFirst | gbkcSecond==TrueFirst) ~ 1,
    T ~ 0)) %>% 
  mutate(Sample=as.integer(gsub("AWA4634-", "", PacBio.ID))) %>% 
  arrange(Sample) %>% 
  mutate(Sample=(paste0("AWA4634-", Sample))) %>% 
  mutate(Method=case_when(
    Method=="median" ~ "flank median",
    Method=="mean" ~ "flank mean",
    T ~ Method
  ))

ohs.calls$Sample <- factor(ohs.calls$Sample, levels=c(unique(ohs.calls$Sample)))

ohs.all.correct <- ohs.calls %>% 
  filter(Rank<=1, NumCorrectAlleles==2)

ohs <- ggplot(ohs.all.correct, aes(k, Sample)) +
  theme_bw() +
  facet_wrap(~Method) +
  geom_point( aes(color=Method), alpha=0.5, show.legend=F) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values=c(palette()[c(1,3)])) +
  xlab("*k*") +
  ylab("OHS sample") +
  theme(axis.title.x = ggtext::element_markdown()) +
  xlim(0, 700)



### HPRC ###
hprc.gbkc.data <- read.table("data/chapter3/gbkc-pacbio-HPRC", header=F, col.names=c("Sample", "Method", "k", "gbkcGenotype", "Score"))

hprc.truth.genos.data <- read.table("data/chapter3/truth-genos-HPRC.tsv",
                                    header=T)

hprc.sample.names.data <- read.table("data/chapter3/sample-names-used-HPRC.txt", header=F, col.names=c("Sample"))

hprc.sample.names <- hprc.sample.names.data %>% 
  pull(Sample)


hprc.truth.genos <- hprc.truth.genos.data %>% 
  mutate(BaseGenotype=gsub("Novel_Similar_|Potential_False_", 
                           "", TruthGenotype)) %>% 
  separate(BaseGenotype, c("TrueFirst", "TrueSecond"), sep="/") %>% 
  mutate(Type=ifelse(TrueFirst==TrueSecond, "Hom", "Het")) 


hprc.calls <- hprc.gbkc.data %>% 
  group_by(Sample, Method, k) %>%
  mutate(Rank=rank(-Score, ties.method="min")) %>% 
  separate(gbkcGenotype, into=c("gbkcFirst", "gbkcSecond"), sep="/", remove=F) %>% 
  filter(Sample!="NA19239") %>% 
  full_join(hprc.truth.genos) %>% 
  rowwise() %>% 
  mutate(NumCorrectAlleles=case_when(
    Type=="Hom" & (gbkcFirst==gbkcSecond & gbkcFirst==TrueFirst) ~ 2,
    Type=="Het" & (grepl(gbkcFirst, TruthGenotype) & 
                     grepl(gbkcSecond, TruthGenotype)) ~ 2,
    Type=="Het" & (grepl(gbkcFirst, TruthGenotype) | 
                     grepl(gbkcSecond, TruthGenotype)) ~ 1,
    Type=="Hom" & (gbkcFirst==TrueFirst | gbkcSecond==TrueFirst) ~ 1,
    T ~ 0)) %>% 
  mutate(Method=case_when(
    Method=="median" ~ "flank median",
    Method=="mean" ~ "flank mean",
    T ~ Method
  ))

hprc.calls$Sample <- factor(hprc.calls$Sample, levels=hprc.sample.names)

hprc.all.correct <- hprc.calls %>% 
  filter(Rank<=1, NumCorrectAlleles==2)

hprc <- ggplot(hprc.all.correct, aes(k, Sample)) +
  theme_bw() +
  facet_wrap(~Method) +
  geom_point( aes(color=Method), alpha=0.5, show.legend=F) +
  scale_y_discrete(limits=rev) +
  scale_color_manual(values=c(palette()[c(1,3,4)])) +
  xlab("*k*") +
  ylab("HPRC++ sample") +
  theme(axis.title.x = ggtext::element_markdown()) +
  xlim(0, 700)

length(unique(hprc.all.correct$Sample))

hprc.all.correct %>% 
  ungroup() %>% 
  select(Sample, Method) %>% 
  distinct() %>% 
  group_by(Method) %>% 
  count()


### COMBINE ###

plot_grid(hprc, NULL, ohs, labels=c("A", "", "B"), ncol=3,
          rel_widths=c(6,0.3,5))



