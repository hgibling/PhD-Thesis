library(dplyr)
library(purrr)

setwd("~/awadalla/PRDM9-Project/HMM/sim-reads/pe-lam/100")

read.scores <- partial(read.csv, header=F, stringsAsFactors=F,
                       col.names=c("Simulated", "Tested", "Score", "Iteration", "Type"))

A.av <- read.scores("A-all-i-87-k-scores-av.csv")
A.sum <- read.scores("A-all-i-87-k-scores-sum.csv")

L2.av <- read.scores("L2-all-i-87-k-scores-av.csv")
L2.sum <- read.scores("L2-all-i-87-k-scores-sum.csv")

L33.av <- read.scores("L33-all-i-87-k-scores-av.csv")
L33.sum <- read.scores("L33-all-i-87-k-scores-sum.csv")


all <- bind_rows(A.av, A.sum, L2.av, L2.sum, L33.av, L33.sum) %>%
    group_by(Type, Simulated, Iteration) %>%
    filter(Score==max(Score)) %>%
    ungroup() %>%
    group_by(Type, Simulated) %>%
    count(Tested)

all.x <- bind_rows(A.av, A.sum, L2.av, L2.sum, L33.av, L33.sum) %>%
    group_by(Type, Simulated, Iteration) %>%
    filter(Score==max(Score)) %>%
    filter(Simulated!=Tested)


### Scores

read.indiv.scores <- partial(read.csv, header=F, stringsAsFactors=F,
                             col.names=c("Simulated", "Tested", "ID", "First",
                                         "Second", "Score", "Distances"))

A.av.ind <- read.indiv.scores("compare-av-sum/A-flank1000-100-bp-100-x-e-0-i-3-F-pe-87-k-scores-verbose-av.csv") %>%
    rename(ScoreAv=Score, DistancesAv=Distances) %>%
    filter(!is.na(ScoreAv)) 

A.av.ind.sp <- A.av.ind %>%
    spread(Tested, ScoreAv)

A.sum.ind <- read.indiv.scores("compare-av-sum/A-flank1000-100-bp-100-x-e-0-i-3-F-pe-87-k-scores-verbose-sum.csv") %>%
    rename(ScoreSum=Score, DistancesSum=Distances) %>%
    filter(!is.na(ScoreSum)) 

A.sum.ind.sp <- A.sum.ind %>%
    spread(Tested, ScoreSum)

A.ind <- full_join(A.av.ind, A.sum.ind, by=c("Simulated", "Tested", "ID", "First", "Second"))
