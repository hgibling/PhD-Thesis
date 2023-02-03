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
        group_by(Aligner, Frag, Allele, Graph, Path, Coverage, Error)
    graphs.num.reads <- graphs %>% 
        summarize(NumReadsTotal=n())
    no.ref <- graphs %>%
        filter(is.na(SelfDeviation)) %>%
        summarize(NumReadsNoMapRef=n())
    no.graph <- graphs %>%
        filter(is.na(Deviation)) %>%
        summarize(NumReadsNoMapGraph=n())

    ## remove reads not aligned in graphs or ref and separate
    graphs.map <- graphs %>%
        filter(!is.na(SelfDeviation), !is.na(Deviation))
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
    graphs.mean.pos <- graphs.map %>%
        filter(Difference>0) %>%
        summarize(MeanDiffPos=mean(Difference))
    graphs.mean.pos.reads <- graphs.map %>%
        filter(Difference>0) %>%
        summarize(NumReadsPos=n())
    graphs.mean.neg <- graphs.map %>%
        filter(Difference<0) %>%
        summarize(MeanDiffNeg=mean(Difference))
    graphs.mean.neg.reads <- graphs.map %>%
        filter(Difference<0) %>%
        summarize(NumReadsNeg=n())
    graphs.all <- graphs.mean %>%
        full_join(graphs.mean.no.0) %>%
        full_join(graphs.reads) %>%
        full_join(graphs.mean.no.0.reads) %>%
        full_join(graphs.mean.pos) %>%
        full_join(graphs.mean.neg) %>%
        full_join(graphs.mean.pos.reads) %>%
        full_join(graphs.mean.neg.reads) %>%
        full_join(no.ref) %>%
        full_join(no.graph) %>%
        full_join(graphs.num.reads) %>%
        mutate(PropReadsDiff=NumReadsNo0/NumReads, 
           PropPos=NumReadsPos/NumReads, 
           PropNeg=NumReadsNeg/NumReads,
           PropReadsNoMapRef=NumReadsNoMapRef/NumReadsTotal,
           PropReadsNoMapGraph=NumReadsNoMapGraph/NumReadsTotal,
           PropReadsDiffTotal=NumReadsNo0/NumReadsTotal,
           PropPosTotal=NumReadsPos/NumReadsTotal, 
           PropNegTotal=NumReadsNeg/NumReadsTotal) %>%
        mutate_all(~replace(., is.na(.), 0))
    return(graphs.all)
}

save.name <- paste0("differences/", args[1], "-self-znf.csv")

output <- get.scores(args[1])

write.csv(output, save.name,  quote=F, row.names=F)
