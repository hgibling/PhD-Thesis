library(dplyr)
library(tidyr)
library(ggplot2)

get.lambda <- function(length, kmer, coverage, error) {
  lam <- (length - kmer + 1) * (coverage / length) * ((1 - error)^kmer)
  return(round(lam, digits=2))
}

get.error <- function(length, kmer, coverage, lambda) {
  err <- 1 - (lambda / ((length - kmer + 1) * (coverage / length)))^(1/kmer)
  return(err)
}

flank.lambdas <- read.csv("data/chapter2/simulated-manual-lambda-flank-mean.csv", 
                          stringsAsFactors=F, header=F,
                          col.names=c("Sample", "Error", "Iteration", "k", "ExpectedLam", "FlankMean")) %>%
  mutate(ExpectedLam=round(ExpectedLam*2, digits=2), FlankMean=FlankMean*2,
         Sample=gsub("-", "/", Sample)) %>% 
  filter(Iteration==1)

all.dat <- read.csv("data/chapter2/simulated-manual-lambda-all-values.csv", 
                stringsAsFactors=F, header=F,
                col.names=c("Sample", "Error", "Iteration", "lambda", "k", "Genotype", "Score")) 
                
all.iter1 <- all.dat %>%
  filter(Iteration==1) %>% 
  group_by(Sample, Error, k, lambda) %>%
  mutate(Rank=rank(-Score, ties.method = "first"), 
         Sample=gsub("-", "/", Sample)) %>%
  mutate(Called=case_when(
    Rank == 1 & Sample == Genotype ~ "True",
    Rank > 1 ~ "False"),
    ExpectedLam=get.lambda(100, k, 100, Error)) %>%
  arrange(Rank) %>% 
  full_join(flank.lambdas) %>%
  mutate(Error=case_when(
    Error==0 ~ "0%",
    Error==0.001 ~ "0.1%",
    Error==0.01 ~"1%"
  ),
  k=paste0(k, "-mers"))

all.iter1$Called <- factor(all.iter1$Called, levels=c("True", "False"))
all.iter1$Error <- factor(all.iter1$Error, levels=c("0%", "0.1%", "1%"))

all.AA <- all.iter1 %>%
  filter(Sample=="A/A", Genotype=="A/A")

all.AL37 <- all.iter1 %>%
  filter(Sample=="A/L37", Genotype=="A/L37")

aa <- ggplot(all.AA, aes(lambda, Rank)) +
  theme_bw() +
  facet_grid(k ~ Error) +
  geom_point(alpha=0.5, aes(color=Called)) +
  geom_vline(aes(xintercept=ExpectedLam, linetype="coverage"), color="black") +
  geom_vline(aes(xintercept=FlankMean, linetype="flank mean"), color=palette()[c(7)]) +
  scale_color_manual(values=c(palette()[c(4,2)]), labels=c("True", "False")) +
  scale_linetype_manual(name="Estimated λ model", values=c(1, 2),
                        guide=guide_legend(override.aes=list(color=c(palette()[c(1,7)])))) +
  guides(color=guide_legend(title="A/A called correctly", order=1)) +
  xlab("Manually specified λ") +
  ylab("Rank of correct genotype")

al <- ggplot(all.AL37, aes(lambda, Rank)) +
  theme_bw() +
  facet_grid(k ~ Error) +
  geom_point(alpha=0.5, aes(color=Called)) +
  geom_vline(aes(xintercept=ExpectedLam, linetype="coverage"), color="black") +
  geom_vline(aes(xintercept=FlankMean, linetype="flank mean"), color=palette()[c(7)]) +
  scale_color_manual(values=c(palette()[c(4,2)]), labels=c("True", "False")) +
  scale_linetype_manual(name="Estimated λ model", values=c(1, 2),
                        guide=guide_legend(override.aes=list(color=c(palette()[c(1,7)])))) +
  guides(color=guide_legend(title="A/L37 called correctly", order=1)) +
  xlab("Manually specified λ") +
  ylab("Rank of correct genotype")

plot_grid(aa, NULL, al, labels=c("A","","B"), ncol=1, align=T, rel_heights=c(1,0.1,1))