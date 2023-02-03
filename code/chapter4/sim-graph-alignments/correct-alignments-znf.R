library(dplyr)
library(tidyr)
#library(ggplot2)

args <- commandArgs(trailingOnly=T)

alignment.dir <- "/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/graphs/alignments/"
setwd(alignment.dir)

graph.names <- c("allele-msa", "allele-msa-dag", "znf-msa", "znf-loop", "stacked")
graph.ref.names <- c("allele-msa", "allele-msa-dag", "znf-msa", "znf-loop", "stacked", "hg19")
allele.names <- c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

sims.nums <- read.table("/.mounts/labs/awadallalab/private/hgibling/PRDM9-Project/THESIS/PRDM9-36/graphs/alignments/all-sims-nums-znf.tsv", header=F, col.names=c("Allele", "TotalZnfReads"))

get.scores <- function(allele) {
    graph.file <- paste0(alignment.dir, allele, "-graph-alignments-znf.csv")
    allele.ref.file <- paste0(alignment.dir, allele, "-self-alignments-znf.csv")

    self.data <- read.csv(allele.ref.file, header=F, stringsAsFactors=F,
                     col.names=c("Aligner", "Frag", "Allele", "Graph", "Path", "Coverage", "Error", "Iteration", "Name", "Deviation")) 
                     
    self.num.spurious <- self.data %>%
        filter(Error==0, Deviation!=0) %>%
        arrange(Deviation) %>%
        distinct() %>%
        group_by(Aligner, Frag, Allele, Path, Coverage, Error) %>%
        count(name="NumReadsSpurious")

    self <- self.data %>%
        # keep only error-free reads with no deviation
        filter(Error==0, Deviation==0) %>%
        distinct() %>%
        select(-Graph, -Deviation)

    data <- read.csv(graph.file, header=F, stringsAsFactors=F,
                       col.names=c("Aligner", "Frag", "Allele", "Graph", "Path", "Coverage", "Error", "Iteration", "Name", "Deviation")) %>%
        filter(Error==0) %>%
        filter(!(Aligner=="vg" & Graph=="allele-msa")) %>%
        # if a read has multiple alignments, keep the one with the lowest deviation
        arrange(Deviation) %>%
        distinct(Aligner, Frag,  Allele, Graph, Path, Coverage, Error, Iteration, Name, .keep_all=T) %>%
        spread(Graph, Deviation) %>%
        # remove reads that had spurious alignment score differences
        right_join(self) %>%
        gather("Graph", "Deviation", all_of(graph.ref.names))
        
    ref <- data %>%
        filter(Graph=="hg19") %>%
        rename(RefDeviation=Deviation)

    graphs <- data %>%
        filter(Graph!="hg19") %>%
        filter(!(Aligner=="vg" & Graph=="allele-msa")) %>%
        spread(Graph, Deviation) %>%
        full_join(ref) %>%
        gather("Graph", "Deviation", all_of(graph.names)) %>%
        mutate(Difference=RefDeviation-Deviation) %>%
        group_by(Aligner, Frag, Allele, Graph, Path, Coverage, Error)

    yes.ref <- graphs %>%
        filter(!is.na(RefDeviation)) %>%
        summarize(NumReadsYesMapRef=n())
    yes.graph <- graphs %>%
        filter(!is.na(Deviation)) %>%
        summarize(NumReadsYesMapGraph=n())
    
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
        full_join(yes.ref) %>%
        full_join(yes.graph) %>%
        full_join(self.num.spurious) %>%
        # join read number info
        left_join(sims.nums) %>%
        mutate_all(~replace(., is.na(.), 0)) %>%
        mutate(NumReadsNoMapGraph=TotalZnfReads- 
            (NumReadsYesMapGraph +
            NumReadsSpurious),
            NumReadsNoMapRef=TotalZnfReads- 
            (NumReadsYesMapRef +
            NumReadsSpurious),
            PropReadsDiff=NumReadsNo0/NumReadsYesMapGraph, 
            PropReadsDiffTotal=NumReadsNo0/TotalZnfReads,
            PropReadsNoMapRef=NumReadsNoMapRef/TotalZnfReads,
            PropReadsNoMapGraph=NumReadsNoMapGraph/TotalZnfReads,
            PropReadsYesMapRef=NumReadsYesMapRef/TotalZnfReads,
            PropReadsYesMapGraph=NumReadsYesMapGraph/TotalZnfReads,
            PropReadsSpuriousTotal=NumReadsSpurious/TotalZnfReads)

    return(graphs.all)
}

save.name <- paste0("differences/", args[1], "-corrected-znf.csv")

output <- get.scores(args[1])

write.csv(output, save.name,  quote=F, row.names=F)