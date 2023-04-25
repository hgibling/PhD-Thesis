library(dplyr)
library(tidyr)
library(tidytext)
library(ggplot2)
library(cowplot)

get.lambda <- function(length, kmer, coverage, error) {
  lam <- (length - kmer + 1) * (coverage / length) * ((1 - error)^kmer)
  return(round(lam, digits=2))
}

rcoverage <- read.table("data/chapter2/manual-lambda/sample-coverage.tsv",
                       stringsAsFactors=F, header=F,
                       col.names=c("Sample", "ReadLength", "Coverage"))


rflex <- read.csv("data/chapter2/manual-lambda/all-flexi-lam-results.csv",
                 stringsAsFactors=F, header=F,
                 col.names=c("Sample", "ReadLength", "lambda", "k", 
                             "Tested", "Score"))

rflank.lambdas <- read.table("data/chapter2/manual-lambda/all-counts-mean-lambda.tsv",
                            stringsAsFactors=F,
                            col.names=c("Sample", "ReadLength", "k", "Flanklambda"))

rflank <- read.csv("data/chapter2/manual-lambda/all-mean-lam-results.csv", 
                  stringsAsFactors=F, header=F,
                  col.names=c("Sample", "ReadLength", "lambda", "k", 
                              "Tested", "Score")) %>%
  select(-lambda) %>%
  full_join(rflank.lambdas) %>%
  mutate(Flanklambda=round(Flanklambda, 2))



### Lambda range plots

rflex.filtered <- rflex %>%
  full_join(rcoverage) %>%
  mutate(Coverage=round(Coverage, 1)) %>%
  mutate(LowerLam=get.lambda(ReadLength, k, Coverage, 0.01),
         HigherLam=get.lambda(ReadLength, k, Coverage, 0)) %>%
  filter(lambda > LowerLam & lambda < HigherLam)

rflex.max <- rflex.filtered %>%
  select(-c(LowerLam, HigherLam)) %>%
  group_by(Sample, ReadLength, lambda, k) %>%
  mutate(Rank=rank(-Score, ties.method="first"), Max=max(Score)) %>%
  mutate(Called=case_when(
    Score == Max ~ T,
    Score != Max ~ F),
    ExpectedLam_0=get.lambda(ReadLength, k, Coverage, 0),
    ExpectedLam_0.001=get.lambda(ReadLength, k, Coverage, 0.001),
    ExpectedLam_0.01=get.lambda(ReadLength, k, Coverage, 0.01)) %>%
  full_join(rflank %>% select(Sample, ReadLength, k, Tested, Flanklambda)) %>%
  arrange(Rank)

### scores for correct genotype for all lambdas

rflex.predicted <- rflex.max %>%
  filter((Sample=="HG002" & Tested=="A/A") | 
           (Sample=="HG003" & Tested=="A/A") |
           (Sample=="HG004" & Tested=="A/L37")) %>%
  arrange(Sample, ReadLength, Rank) %>%
  mutate(k=paste0(k, "-mers"),
         k=factor(k, levels=c("51-mers", "71-mers", "91-mers", "111-mers", "131-mers")))



# Ranks
raw <- ggplot(rflex.predicted %>% filter(ReadLength==150), aes(lambda, Rank)) +
  theme_bw() +
  facet_grid(k ~ Sample) +
  geom_point(alpha=0.5, aes(color=Called)) +
  geom_vline(aes(xintercept=ExpectedLam_0, linetype="coverage 0.0% error"), color="black") +
  geom_vline(aes(xintercept=ExpectedLam_0.001, linetype="coverage 0.1% error"), color="black") +
  geom_vline(aes(xintercept=ExpectedLam_0.01, linetype="coverage 1.0% error"), color="black") +
  geom_vline(aes(xintercept=Flanklambda, linetype="flank mean"), color=palette()[c(7)]) +
  scale_color_manual(values=c(palette()[c(2)]), labels=c("False")) +
  scale_linetype_manual(name="Estimated λ model", values=c(1, 2, 3, 2),
                        guide=guide_legend(override.aes=list(color=c(rep("black", 3), palette()[c(7)])))) +
  guides(color=guide_legend(title="Correctly called", order=1)) +
  xlab("Manually specified λ") +
  ylab("Rank of truth genotype") +
  ggtitle("Raw reads")

rflex.predicted %>% filter(ReadLength==150) %>% ungroup() %>% 
  group_by(Sample, ReadLength) %>% summarize(min(Rank))



coverage <- read.table("data/chapter2/manual-lambda/sample-coverage.tsv",
                       stringsAsFactors=F, header=F,
                       col.names=c("Sample", "ReadLength", "Coverage"))


flex <- read.csv("data/chapter2/manual-lambda/all-flexi-lam-results-aligncorrect.csv",
                 stringsAsFactors=F, header=F,
                 col.names=c("Sample", "ReadLength", "EditDistance", "lambda", "k", 
                             "Tested", "Score"))

flank.lambdas <- read.table("data/chapter2/manual-lambda/all-counts-mean-lambda-aligncorrect.tsv",
                            stringsAsFactors=F,
                            col.names=c("Sample", "ReadLength", "EditDistance", "k", "Flanklambda"))

flank <- read.csv("data/chapter2/manual-lambda/all-mean-lam-results-aligncorrect.csv", 
                  stringsAsFactors=F, header=F,
                  col.names=c("Sample", "ReadLength", "EditDistance", "lambda", "k", 
                              "Tested", "Score")) %>%
  select(-lambda) %>%
  full_join(flank.lambdas) %>%
  mutate(Flanklambda=round(Flanklambda, 2))



### Lambda range plots

flex.filtered <- flex %>%
  full_join(coverage) %>%
  mutate(Coverage=round(Coverage, 1)) %>%
  mutate(LowerLam=get.lambda(ReadLength, k, Coverage, 0.01),
         HigherLam=get.lambda(ReadLength, k, Coverage, 0)) %>%
  filter(lambda > LowerLam & lambda < HigherLam)

flex.max <- flex.filtered %>%
  select(-c(LowerLam, HigherLam)) %>%
  group_by(Sample, ReadLength, EditDistance, lambda, k) %>%
  mutate(Rank=rank(-Score, ties.method="first"), Max=max(Score)) %>%
  mutate(Called=case_when(
    Score == Max ~ T,
    Score != Max ~ F),
    ExpectedLam_0=get.lambda(ReadLength, k, Coverage, 0),
    ExpectedLam_0.001=get.lambda(ReadLength, k, Coverage, 0.001),
    ExpectedLam_0.01=get.lambda(ReadLength, k, Coverage, 0.01)) %>%
  full_join(flank %>% select(Sample, ReadLength, EditDistance, k, Tested, Flanklambda)) %>%
  arrange(Rank)


### scores for correct genotype for all lambdas

flex.predicted.cor <- flex.max %>%
  filter((Sample=="HG002" & Tested=="A/A") | 
           (Sample=="HG003" & Tested=="A/A") |
           (Sample=="HG004" & Tested=="A/L37")) %>%
  arrange(Sample, ReadLength, Rank) %>%
  mutate(EditDistance=as.factor(EditDistance),
         k=paste0(k, "-mers"),
         k=factor(k, levels=c("51-mers", "71-mers", "91-mers", "111-mers", "131-mers")))

flex.predicted.cor$Called <- factor(flex.predicted.cor$Called, levels=c(T,F))

# Ranks

cor <- ggplot(flex.predicted.cor %>% filter(ReadLength==150, EditDistance==40), aes(lambda, Rank)) +
  theme_bw() +
  facet_grid(k ~ Sample) +
  geom_point(alpha=0.5, aes(color=Called)) +
  geom_vline(aes(xintercept=ExpectedLam_0, linetype="coverage 0.0% error"), color="black") +
  geom_vline(aes(xintercept=ExpectedLam_0.001, linetype="coverage 0.1% error"), color="black") +
  geom_vline(aes(xintercept=ExpectedLam_0.01, linetype="coverage 1.0% error"), color="black") +
  geom_vline(aes(xintercept=Flanklambda, linetype="flank mean"), color=palette()[c(7)]) +
  scale_color_manual(values=c(palette()[c(4,2)]), labels=c("True", "False")) +
  scale_linetype_manual(name="Estimated λ model", values=c(1, 2, 3, 2),
                        guide=guide_legend(override.aes=list(color=c(rep("black", 3), palette()[c(7)])))) +
  guides(color=guide_legend(title="Correctly called", order=1)) +
  xlab("Manually specified λ") +
  ylab("Rank of truth genotype") +
  ggtitle("AlignCorrected reads")

plot_grid(raw, cor, ncol=1, labels=c("A", "B"))


