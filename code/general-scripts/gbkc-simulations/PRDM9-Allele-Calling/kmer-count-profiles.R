library(dplyr)
library(tidyr)
library(ggplot2)

setwd("~/Documents/PRDM9-Allele-Calling")

znf <- read.csv("~/Documents/PRDM9-Project/allele-znf-content-aligned.csv",
                stringsAsFactors=F, header=F, col.names=c("Tested", "ZnF"))

znf$ZnF <- gsub('-', '', znf$ZnF)

znf <- znf %>%
    group_by(Tested) %>%
    mutate(Length=nchar(ZnF), Repeat=length(which(duplicated(unlist(strsplit(ZnF, ''))))==T)) %>%
    arrange(Length)
znf$Length <- as.factor(znf$Length)
znf$Repeat <- as.factor(znf$Repeat)



order <-c(LETTERS[1:5], paste0("L", 1:24), paste0("L", 32:38))

scores.N <- read.csv("all-alleles-flank1000-lamx-e1-e001-cov40-k51count-score.csv",
                     header=F, stringsAsFactors=F,
                     col.names=c("Simulated", "Tested", "Score"))

scores <- scores.N %>%
    filter(Simulated!="N", Tested!="N") %>%
    spread(Tested, Score) %>%
    arrange(factor(Simulated, levels=rev(order))) %>%
    select(order)
rownames(scores) <- rev(colnames(scores))
    
levelplot(t(as.matrix(scores)), scale=list(x=list(rot=45)),
          col.regions=colorRampPalette(c("#004949","#009292","#ff6db6"))(1e4),
          aspect="fill", 
          xlab="Tested Allele \n",
          ylab="Simulated Allele",
          colorkey=list(space='bottom'),
          main="Scores of tested alleles (log probability)"
)

scores.binary <- t(apply(scores, 1, function(z){
    1 * (z==max(z))
}))

levelplot(t(as.matrix(scores.binary)), scale=list(x=list(rot=45)),
          col.regions=colorRampPalette(c("white", "black"))(2),
          aspect="fill", 
          xlab="Tested Allele \n",
          ylab="Simulated Allele",
          colorkey=list(space='bottom'),
          main="Tested alleles with highest score"
)

all.scores.N <- group_by(scores.N, Simulated) %>%
    arrange(desc(Score)) %>%
    arrange(Simulated) %>%
    mutate(Correct=ifelse(Simulated==Tested, T, F)) %>%
    group_by(Simulated) %>%
    mutate(Unique=ifelse(Score %in% Score[duplicated(Score)], F, T))

all.scores <- filter(all.scores.N, Simulated!="N", Tested!="N")  %>%
    full_join(znf, by='Tested') %>%
    arrange(desc(Score))

all.scores.max <- filter(all.scores, Score==max(Score))
View(all.scores.max)

norm <- all.scores %>%
    mutate(NormScore=Score*as.numeric(Length))


# test

ggplot(filter(all.scores, Simulated=="L6"), 
       aes(reorder(Tested, -Score), Score)) +
    geom_point(size=5, aes(color=Correct)) +
    theme_bw() +
    xlab("Allele Tested") +
    ylab("Log Score")

# plot all
for (i in unique(all.scores$Simulated)) {
    ggplot(filter(all.scores, Simulated==i), 
           aes(reorder(Tested, -Score), Score)) +
        geom_point(size=5, aes(color=Correct)) +
        scale_shape_manual(values=c(15, 16)) +
        theme_bw() +
        xlab("Allele Tested") +
        ylab("Log Score")
    ggsave(paste0("plots/kmer/lam9e1/", i, "-score.pdf"), width=15, height=5)
}


ggplot(filter(norm, Simulated=="A"), 
       aes(reorder(Tested, -NormScore), NormScore)) +
    geom_point(size=5, aes(color=Length, shape=Correct)) +
    theme_bw() +
    xlab("Allele Tested (Norm)") +
    ylab("Log Score")


scores.norm <- norm %>%
    spread(Tested, NormScore) %>%
    arrange(factor(Simulated, levels=rev(order))) %>%
    ungroup() %>%
    select(order)
rownames(scores.norm) <- rev(colnames(scores.norm))

levelplot(t(as.matrix(scores.norm)), scale=list(x=list(rot=45)),
          col.regions=colorRampPalette(c("#004949","#009292","#ff6db6"))(1e4),
          aspect="fill", 
          xlab="Tested Allele \n",
          ylab="Simulated Allele",
          colorkey=list(space='bottom'))


# kmers

allele.count <- read.csv("all-alleles.k51.count",
         header=F, stringsAsFactors=F,
         col.names=c("Allele", "kmer", "AlleleCount"))

read.count <- read.csv("all.flank1000.perfect.100bp.cov40.F.reads.kmers.count",
                         header=F, stringsAsFactors=F,
                         col.names=c("Allele", "kmer", "ReadCount"))

ar.count <- full_join(allele.count, read.count, by=c("Allele", "kmer"))
ar.count[is.na(ar.count)] <- 0



universal <- filter(ar.count, AlleleCount==1) %>%
    count(kmer) %>%
    arrange(desc(n))

universal.all <- filter(universal, n==36)

universal.all.ar <- filter(ar.count, kmer %in% universal.all$kmer) %>%
    group_by(Allele) %>%
    summarize(MedianRead=median(ReadCount)) %>%
    arrange(desc(MedianRead))


ggplot(filter(ar.count, Allele=="N"), 
       aes(reorder(kmer, -Difference), Difference)) +
    geom_point(size=2) +
    theme_bw() +
    geom_hline(yintercept=0, color='red') +
    xlab("kmer") +
    ylab("Difference ((Allele count *20) - Read count)") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


count.sum <- allele.count %>%
    count(Allele) %>%
    mutate(Tested=Allele)

norm.kmer <- all.scores %>%
    full_join(count.sum, by="Tested") %>%
    mutate(NormScore=Score/(n))
norm.kmer$n <- as.factor(norm.kmer$n)

ggplot(filter(all.scores, Simulated=="A"), 
       aes(reorder(Tested, -Score), Score)) +
    geom_point(size=5, aes(color=Length, shape=Correct)) +
    theme_bw() +
    xlab("Allele Tested") +
    ylab("Log Score")

ggplot(filter(norm.kmer, Simulated=="A"), 
       aes(reorder(Tested, -NormScore), NormScore)) +
    geom_point(size=5, aes(color=n, shape=Correct)) +
    theme_bw() +
    xlab("Allele Tested (Norm)") +
    ylab("Log Score")



# check sim reads

sim <- read.csv("all.alleles.pos.mysim.flank.perfect.100bp.F.reads", header=F,
                stringsAsFactors=F, col.names=c('Allele', 'Pos', 'Read'))

sim <- mutate(sim, n=1:n())

ggplot(filter(sim, Allele=="E"), aes(Pos, reorder(n, -Pos))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(breaks=sim$n, labels=sim$Reads) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

ggplot(filter(sim, Allele=="L22"), aes(Pos, reorder(n, -Pos))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(breaks=sim$n, labels=sim$Reads) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())



sim.vg <- read.csv("test-E-sim.csv", header=F,
                stringsAsFactors=F, col.names=c('Read', 'Pos'))

sim.vg <- mutate(sim.vg, n=1:n())
sim.vg$Pos <- as.numeric(sim.vg$Pos)

ggplot(sim.vg, aes(Pos, reorder(n, -Pos))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(breaks=sim.vg$n, labels=sim.vg$Reads) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())


sim.vg.L22 <- read.csv("test-L22-sim.csv", header=F,
                   stringsAsFactors=F, col.names=c('Read', 'Pos', 'Node'))
sim.vg.L22$Pos <- as.numeric(sim.vg.L22$Pos)

sim.vg.L22 <- mutate(sim.vg.L22, n=1:n(), 
                     TruePos=ifelse(Node==2, Pos+688, Pos))

ggplot(sim.vg.L22, aes(TruePos, reorder(n, -TruePos))) +
    geom_point() +
    theme_bw() +
    scale_y_discrete(breaks=sim.vg.L22$n, labels=sim.vg.L22$Reads) +
    theme(axis.text.y=element_blank(),
          axis.ticks.y=element_blank())



# kmers for flank1000

flank.count <- read.csv("both-flank1000.k51.count",
                         header=F, stringsAsFactors=F,
                         col.names=c("Allele", "kmer", "AlleleCount"))

read.count <- read.csv("all.flank1000.perfect.100bp.cov40.F.reads.kmers.count",
                       header=F, stringsAsFactors=F,
                       col.names=c("Allele", "kmer", "ReadCount"))

universal <- flank.count$kmer

universal.fr <- filter(read.count, kmer %in% universal) %>%
    group_by(Allele) %>%
    summarize(MedianRead=median(ReadCount)) %>%
    arrange(desc(MedianRead))
