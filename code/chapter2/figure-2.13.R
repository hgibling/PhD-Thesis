library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggtext)

## COUNT
new.count.dat <- read.table("data/chapter2/all-count-average-dip-F1.tsv", header=T)

new.count <- new.count.dat %>%
  mutate(Coverage=paste0(Coverage, "X"),
         Error=case_when(
           Error==0 ~ "0%",
           Error==0.001 ~ "0.1%",
           Error==0.01 ~ "1%"
         ))


new.count$Coverage <- factor(new.count$Coverage, 
                             levels=c("20X", "40X", "60X", "80X", "100X"))
new.count$Error <- factor(new.count$Error, 
                          levels=c("0%", "0.1%", "1%"))




## DISTANCE

dip.distance.dat <- read.csv("data/chapter2/all-distance-average-dip-F1.tsv",
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
         Ploidy="diploid") %>%
  mutate(Coverage=paste0(Coverage, "X"),
         Error=case_when(
           Error==0 ~ "0%",
           Error==0.001 ~ "0.1%",
           Error==0.01 ~ "1%"
         ))

dip.distance$Coverage <- factor(dip.distance$Coverage, 
                                levels=c("20X", "40X", "60X", "80X", "100X"))
dip.distance$Error <- factor(dip.distance$Error, 
                             levels=c("0%", "0.1%", "1%"))

new.count.dist <- new.count %>% 
  filter(Method=="coverage") %>%
  mutate(Method=paste0("count-", Method)) %>% 
  bind_rows(dip.distance %>% 
              filter(Method=="max") %>%
              mutate(Method=paste0("distance-", Method)))






dat.com <- read.csv("data/chapter2/av-combine-diploid-flank10k.csv",
                    header=F, col.names=c("kmer", "Average", "Precision", "Recall",
                                          "F1Score", "Coverage", "Error", "Iteration", "ReadLength", 
                                          "FragmentLength", "Method"))

dat.com.av <- dat.com %>%
  filter(Method=="max") %>% 
  group_by(kmer, ReadLength, FragmentLength, Coverage, Error, Method) %>%
  summarize(AveragePrecision=mean(Precision), 
            AverageRecall=mean(Recall), 
            AverageF1=mean(F1Score)) %>%
  rename(k=kmer) %>%
  ungroup() %>%
  select(-AveragePrecision, -AverageRecall) %>%
  mutate(Coverage=paste0(Coverage, "X"),
         Error=case_when(
           Error==0 ~ "0%",
           Error==0.001 ~ "0.1%",
           Error==0.01 ~ "1%"
         ))



### NICE COMBINED PLOT

dat.count.dist.com <- new.count %>%
  ungroup() %>%
  filter(Method=="coverage") %>%
  select(Method, k, Coverage, Error, AverageF1) %>%
  mutate(Method="count-coverage") %>%
  bind_rows(dip.distance %>% 
              ungroup() %>%
              mutate(Method="distance-max") %>% 
              select(Method, k, Coverage, Error, AverageF1)) %>% 
  bind_rows(dat.com.av %>% 
              ungroup() %>%
              mutate(Method="count-coverage\n& distance-max") %>% 
              select(Method, k, Coverage, Error, AverageF1))

dat.count.dist.com$Coverage <- factor(dat.count.dist.com$Coverage,
                                      levels=c("20X", "40X", "60X", "80X", "100X"))
dat.count.dist.com$Error <- factor(dat.count.dist.com$Error,
                                      levels=c("0%", "0.1%", "1%"))


ggplot(dat.count.dist.com, aes(k, AverageF1)) +
  facet_grid(Coverage~Error) +
  theme_bw() +
  geom_point(aes(color=Method), alpha=0.5) +
  geom_line(aes(color=Method), alpha=0.5) +
  scale_color_manual(values=c(palette()[c(1,4,2)])) +
  xlab("*k*") +
  ylab("Average F1 score") +
  theme(axis.title.x = ggtext::element_markdown()) +
  coord_cartesian(ylim=c(0.3,0.8))
