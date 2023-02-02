### analyze pseudo-F1 scores for graph k-mer precision and recall

# load libraries
library(dplyr)
library(tidyr)

# load data file for allele kmers
allele <- read.table("data/chapter4/graph-kmers/all-alleles-unique-k.tsv", col.names=c("k", "kmer"))

# get number of distinct kmers for each k
allele.count <- allele %>% 
  arrange(k, kmer) %>%
  group_by(k) %>%
  count(name="allele.num")

# load data files for graph kmers
allele.msa <- read.table("allele-msa-no-all-k-unique.tsv", 
                         col.names=c("graph", "k", "kmer"))
allele.msa.dag <- read.table("allele-msa-dag-no-all-k-unique.tsv", 
                             col.names=c("graph", "k", "kmer"))
stacked <- read.table("stacked-no-all-k-unique.tsv", 
                      col.names=c("graph", "k", "kmer"))
znf.loop <- read.table("znf-loop-no-all-k-unique.tsv", 
                       col.names=c("graph", "k", "kmer"))
znf.msa <- read.table("znf-msa-no-all-k-unique.tsv", 
                      col.names=c("graph", "k", "kmer"))
grc <- read.table("GRCh38-all-k-unique.tsv", 
                  col.names=c("graph", "k", "kmer"))

# function to get list of kmers in both allele list and graph list and calculate pseudo-F1 score
get.matches <- function(graph.data) {
  graph.count <- graph.data %>%
    group_by(graph, k) %>%
    count(name="graph.num")
  out <- allele %>% 
    inner_join(graph.data) %>% 
    group_by(graph, k) %>% 
    count(name="match.num") %>% 
    full_join(allele.count) %>% 
    full_join(graph.count) %>% 
    rowwise() %>% 
    mutate(precision=match.num/graph.num,
         recall=match.num/allele.num) %>% 
    mutate(F1=2*((precision*recall)/(precision+recall)))
  return(out)
}

# get F1 scores for each graph and GRCh38
allele.msa.matches <- get.matches(allele.msa)
allele.msa.dag.matches <- get.matches(allele.msa.dag)
stacked.matches <- get.matches(stacked)
znf.loop.matches <- get.matches(znf.loop)
znf.msa.matches <- get.matches(znf.msa)
grc.matches <- get.matches(grc)

# combine into dataframe
all.graphs <- bind_rows(allele.msa.matches, allele.msa.dag.matches, stacked.matches,
                        znf.loop.matches, znf.msa.matches, grc.matches)

# save as table
write.table(all.graphs, "data/chapter4/all-graphs-kmer-F1.tsv", quote=F, row.names=F, sep="\t")