library(dplyr)
library(tidyr)
#library(ggplot2)

args <- commandArgs(trailingOnly=T)

#setwd("prdm9/graphs/alignments/jq")
setwd("/.mounts/labs/simpsonlab/users/hgibling/alignments/compare")

graph.dir <- "/.mounts/labs/simpsonlab/users/hgibling/alignments/"

graph.names <- c("stacked", "allele-msa", "allele-msa-dag", "znf-loop", "znf-msa")
allele.names <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38)) 

get.scores <- function(allele) {
    file <- paste0(graph.dir, allele, "-graph-alignments.csv")
    data <- read.csv(file, header=F, stringsAsFactors=F,
                     col.names=c("Allele", "Name", "Length", "Score", "Graph", "Coverage", "Error", "Iteration")) %>%
        mutate(Score=as.integer(Score)) %>%
        distinct()
    #head(data)
    ref <- data %>%
        filter(Graph=="hg19") %>%
        mutate(NormRefScore=Score/Length) %>%
        select(-Graph, -Length, -Score)
    graphs <- data %>% 
        filter(Graph!="hg19") %>%
        mutate(NormScore=Score/Length) %>%
        select(-Score) %>%
        spread(Graph, NormScore) %>%
        full_join(ref) %>%
        gather("Graph", "NormScore", graph.names) %>%
        mutate(Difference=NormScore-NormRefScore) %>%
        group_by(Allele, Coverage, Error, Graph)
    graphs.mean <- graphs %>%
        summarize(MeanDiff=mean(Difference))
    graphs.reads <- graphs %>%
        summarize(NumReads=n())
    graphs.mean.no.0 <- graphs %>%
        filter(Difference!=0) %>%
        summarize(MeanDiffNo0=mean(Difference))
    graphs.mean.no.0.reads <- graphs %>%
        filter(Difference!=0) %>%
        summarize(NumReadsNo0=n())
    graphs.both <- graphs.mean %>%
        full_join(graphs.mean.no.0) %>%
        full_join(graphs.reads) %>%
        full_join(graphs.mean.no.0.reads) %>%
        mutate(PropReadsDiff=NumReadsNo0/NumReads)
    return(graphs.both)
}

save.name <- paste0(args[1], ".csv")

output <- get.scores(args[1])

write.csv(output, save.name,  quote=F, row.names=F)
