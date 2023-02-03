library(dplyr)
library(tidyr)
#library(ggplot2)

args <- commandArgs(trailingOnly=T)

alignment.dir <- "/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/graphs/alignments/"
setwd(alignment.dir)

graph.names <- c("stacked", "allele-msa", "allele-msa-dag", "znf-loop", "znf-msa")
allele.names <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

get.scores <- function(allele) {
    file <- paste0(alignment.dir, allele, "-graph-alignments-znf.csv")
    data <- read.csv(file, header=F, stringsAsFactors=F,
                     col.names=c("Aligner", "Frag", "Allele", "Graph", "Path", "Coverage", "Error", "Iteration", "Name", "Deviation")) %>%
        # if a read has multiple alignments, keep the one with the lowest deviation
        arrange(Deviation) %>%
        distinct(Aligner, Frag,  Allele, Graph, Path, Coverage, Error, Iteration, Name, .keep_all=T)
    ref <- data %>%
        filter(Graph=="hg19") %>%
        mutate(RefDeviation=Deviation) %>%
        select(-Graph, -Deviation)
    graphs <- data %>%
        filter(Graph!="hg19") %>%
        filter(!(Aligner=="vg" & Graph=="allele-msa")) %>%
        spread(Graph, Deviation) %>%
        full_join(ref) %>%
        gather("Graph", "Deviation", graph.names) %>%
        mutate(Difference=RefDeviation-Deviation) %>%
        group_by(Aligner, Frag, Allele, Graph, Path, Coverage, Error)
    graphs.num.reads <- graphs %>% 
        summarize(NumReadsTotal=n())
    no.ref <- graphs %>%
        filter(is.na(RefDeviation)) %>%
        summarize(NumReadsNoMapRef=n())
    no.graph <- graphs %>%
        filter(is.na(Deviation)) %>%
        summarize(NumReadsNoMapGraph=n())

## remove reads not aligned in graphs or ref and separate
    graphs.map <- graphs %>%
        filter(!is.na(RefDeviation), !is.na(Deviation))
    graphs.mean <- graphs.map %>%
        summarize(MeanDiff=mean(Difference))
    graphs.reads <- graphs.map %>%
        summarize(NumReads=n())
    graphs.mean.no.0 <- graphs.map %>%
        filter(Difference!=0) %>%
        summarize(MeanDiffNo0=mean(Difference))
    graphs.mean.no.0.reads <- graphs.map %>%
        filter(Difference!=0) %>%
        summarize(NumReadsNo0=n())
    graphs.all <- graphs.mean %>%
        full_join(graphs.mean.no.0) %>%
        full_join(graphs.reads) %>%
        full_join(graphs.mean.no.0.reads) %>%
        full_join(no.ref) %>%
        full_join(no.graph) %>%
        full_join(graphs.num.reads) %>%
        mutate(PropReadsDiff=NumReadsNo0/NumReads, 
           PropReadsNoMapRef=NumReadsNoMapRef/NumReadsTotal,
           PropReadsNoMapGraph=NumReadsNoMapGraph/NumReadsTotal,
           PropReadsDiffTotal=NumReadsNo0/NumReadsTotal) %>%
        mutate_all(~replace(., is.na(.), 0))
    return(graphs.all)
}

save.name <- paste0("differences/", args[1], "-znf.csv")

output <- get.scores(args[1])

write.csv(output, save.name,  quote=F, row.names=F)
