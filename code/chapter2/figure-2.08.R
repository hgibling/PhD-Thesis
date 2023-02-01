library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggtext)

new.count.dat <- read.table("data/chapter2/all-count-average-dip-F1.tsv", header=T)
new.count.geno.dat <- read.table("data/chapter2/all-count-genotype-dip-F1.tsv", header=T)

new.count <- new.count.dat %>%
  mutate(Coverage=paste0(Coverage, "X"),
         Error=case_when(
           Error==0 ~ "0%",
           Error==0.001 ~ "0.1%",
           Error==0.01 ~ "1%"
         ))

new.count.geno <- new.count.geno.dat %>%
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

new.count.geno$Coverage <- factor(new.count.geno$Coverage, 
                                  levels=c("20X", "40X", "60X", "80X", "100X"))
new.count.geno$Error <- factor(new.count.geno$Error, 
                               levels=c("0%", "0.1%", "1%"))


pa <- ggplot(new.count %>% filter(Method=="coverage"), aes(k, AverageF1)) +
  theme_bw() +
  facet_grid(Error ~.) +
  geom_point(aes(color=Coverage), alpha=0.5) +
  geom_line(aes(color=Coverage), alpha=0.5) +
  scale_color_manual(values=c(palette()[c(1, 2,3, 4, 6)])) +
  xlab("*k*") +
  ylab("Average F1 score") +
  theme(axis.title.x = ggtext::element_markdown())

pb <- ggplot(new.count.geno %>% filter(k==51, Coverage=="100X", Method=="coverage"), 
             aes(Error, GenotypeF1)) +
  theme_bw() +
  facet_grid(~Coverage) +
  geom_jitter(aes(color=Error),height=0, alpha=0.25, show.legend=F) +
  scale_color_manual(values=c(palette()[c(2:4)])) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
  xlab("Sequencing error rate") + 
  ylab("Average F1 score")


bottom <- plot_grid(pb,NULL,
                    nrow=1, rel_widths=c(6.3,1))


plot_grid(pa,NULL,bottom, ncol=1, align="lr", axis="hv",
          rel_heights=c(3,0.1,2),
          labels=c("A","","B"))
