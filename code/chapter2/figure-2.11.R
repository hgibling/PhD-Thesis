library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggtext)

####### HAPLOID #######

### COUNT ###

hap.data <- read.table("data/chapter2/all-average-hap-F1.tsv", header=T)

hap.count <- hap.data %>% 
  filter(Method=="coverage") %>%
  select(Method, Coverage, Error, k, AverageF1) %>%
  mutate(Method="count-coverage",
         Ploidy="haploid")

### DISTANCE ###

hap.distance <- hap.data %>% 
  filter(Method=="max") %>%
  select(Method, Coverage, Error, k, AverageF1) %>%
  mutate(Method="distance-max",
         Ploidy="haploid")


####### DIPLOID #######
### COUNT ###

dip.count.dat <- read.table("data/chapter2/all-count-average-dip-F1.tsv", header=T)

dip.count <- dip.count.dat %>%
  filter(Method=="coverage") %>%
  select(Method, Coverage, Error, k, AverageF1) %>%
  mutate(Method="count-coverage",
         Ploidy="diploid")


### DISTANCE ###
dip.distance.dat <- read.table("data/chapter2/all-distance-average-dip-F1.tsv",
                             header=F, col.names=c("kmer", "Average", "Precision", "Recall",
                                                   "F1Score", "Coverage", "Error", "Iteration",
                                                   "ReadLength", "FragmentLength", "Method"))

dip.distance <- dip.distance.dat %>%
  filter(Method=="max") %>% 
  group_by(kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
  summarize(AverageF1=mean(F1Score)) %>% 
  ungroup() %>%
  select(Method, Coverage, Error, kmer, AverageF1) %>%
  rename(k=kmer) %>%
  mutate(Method="distance-max",
         Ploidy="diploid")


###### COMBINE PLOIDY ########

all.comb <- hap.count %>%
  bind_rows(hap.distance) %>%
  bind_rows(dip.count) %>%
  bind_rows(dip.distance) %>%
  mutate(Coverage=paste0(Coverage, "X"),
         Error=case_when(
           Error==0 ~ "0%",
           Error==0.001 ~ "0.1%",
           Error==0.01 ~ "1%"
         ))

all.comb$Coverage <- factor(all.comb$Coverage, 
                            levels=c("20X", "40X", "60X", "80X", "100X"))
all.comb$Error <- factor(all.comb$Error,
                         levels=c("0%", "0.1%", "1%"))
all.comb$Ploidy <- factor(all.comb$Ploidy, levels=c("haploid", "diploid"))

ggplot(all.comb, aes(k, AverageF1)) +
  facet_grid(Coverage~Error) +
  theme_bw() +
  geom_point(aes(color=Method, shape=Ploidy, alpha=Ploidy)) +
  geom_line(aes(color=Method, linetype=Ploidy, alpha=Ploidy)) +
  scale_color_manual(values=c(palette()[c(1:2)])) +
  scale_alpha_manual(values=c(0.5,0.5)) +
  scale_shape_manual(values=c(19,1)) +
  xlab("*k*") +
  ylab("Average F1 score") +
  theme(axis.title.x = ggtext::element_markdown())
