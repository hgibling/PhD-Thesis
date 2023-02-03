library(dplyr)
library(tidyr)
#library(ggplot2)

args <- commandArgs(trailingOnly=T)

alignment.dir <- "/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/graphs/alignments/"
setwd(alignment.dir)

graph.names <- c("allele-msa", "allele-msa-dag", "znf-msa", "znf-loop", "stacked")
allele.names <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

get.scores <- function(allele) {
    graph.file <- paste0(alignment.dir, allele, "-graph-alignments-znf.csv")
    allele.ref.file <- paste0(alignment.dir, allele, "-self-alignments-znf.csv")
    ref <- read.csv(allele.ref.file, header=F, stringsAsFactors=F,
                     col.names=c("Aligner", "Frag", "Allele", "Graph", "Path", "Coverage", "Error", "Iteration", "Name", "Deviation")) %>%
        # if a read has multiple alignments, keep the one with the lowest deviation
        arrange(Deviation) %>%
        distinct(Aligner, Frag,  Allele, Graph, Path, Coverage, Error, Iteration, Name, .keep_all=T) %>%
        select(-Graph) %>%
        rename(SelfDeviation=Deviation)
    graphs <- read.csv(graph.file, header=F, stringsAsFactors=F,
                       col.names=c("Aligner", "Frag", "Allele", "Graph", "Path", "Coverage", "Error", "Iteration", "Name", "Deviation")) %>%
        filter(Graph!="hg19") %>%
        filter(!(Aligner=="vg" & Graph=="allele-msa")) %>%
        # if a read has multiple alignments, keep the one with the lowest deviation
        arrange(Deviation) %>%
        distinct(Aligner, Frag,  Allele, Graph, Path, Coverage, Error, Iteration, Name, .keep_all=T) %>%
        spread(Graph, Deviation) %>%
        full_join(ref) %>%
        gather("Graph", "Deviation", graph.names) %>%
        mutate(Difference=SelfDeviation-Deviation) %>%
        filter(Difference>=0.2 | Difference<=-0.2)

    return(graphs)
}

save.name <- paste0("indiv-reads/", args[1], "-self-znf.csv")

output <- get.scores(args[1])

write.csv(output, save.name,  quote=F, row.names=F)

