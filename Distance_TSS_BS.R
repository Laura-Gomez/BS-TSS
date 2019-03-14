
library(ggplot2)

args = commandArgs(trailingOnly=TRUE)

dir.data = args[7]
bs.file = paste(dir.data, args[1], sep="/")
tf.name = args[2]
effect.file <- paste(dir.data, args[3], sep="/")
target.file <- paste(dir.data, args[4], sep="/")
tss.file <- paste(dir.data, args[5], sep="/")
ht.tss <- paste(dir.data, args[6], sep="/")

#bs.file <- c("BindingSiteSet.txt")
#tf.name <- c("OmpR")
#effect.file <- c("OmpR-RNASeq-PeakTargets.csv")
#target.file <- c("OmpR-All.csv")
#tss.file <- c("pm_w_first_transc_g__w_ids_noHT.txt")
#ht.tss <- c("StorzG_TSS_Table_M63_0.4.txt")
#### READ BINDING SITE DATA FOR TF

bs.rdb <- read.table(bs.file, header=F, sep="\t", stringsAsFactors = F)
names(bs.rdb) <- c("TF.ID", "TF", "TFBS.ID", "TF.LEFT", "TF.RIGTH", "TFBS.SATRND", "TF.GENE.ID", "TU", "EFFECT", "PROMOTER", "DIST.TSS", "SEQ", "EVIDENCE", "EVIDENCE.LEVEL")

bs.ompr <- subset(bs.rdb, TF == tf.name)


#### BINDING SITES SUBSETS

bs.ompr.strong <- subset(bs.ompr, EVIDENCE.LEVEL == "Strong")
bs.ompr.activator <- subset(bs.ompr, EFFECT == "+")
bs.ompr.repressor <- subset(bs.ompr, EFFECT == "-")


#### PLOTS FROM DISTANCE TO TSS FOR THIS TF

pdf(file = paste(tf.name, "pdf", sep="."), onefile = T)

{hist(bs.ompr$DIST.TSS, main = "All Binding Sites",
      xlab = "", ylab = "", col = "gold", cex.axis=2, cex.lab = 3, xaxt="n")}

{hist(bs.ompr$DIST.TSS, breaks = seq(-450, 150, by = 100), main = "REGION -450:150",
      xlab = "", ylab = "", col = "gold", cex.axis=2, cex.lab = 3, xaxt="n")
  axis(1, at=seq(-400, 100, by = 100), labels=seq(-400, 100, by = 100), cex.axis=2, las = 2)}


{hist(bs.ompr.strong$DIST.TSS, breaks = seq(-450, 150, by = 100), main = "Strong Evidence",
      xlab = "Distance", ylab = "Number of TFBS", col = "Maroon", cex.axis = 0.8, las = 2, xaxt="n")
  axis(1, at=seq(-400, 100, by = 100), labels=seq(-400, 100, by = 100), cex=0.8, las = 2)}

{hist(bs.ompr.activator$DIST.TSS, breaks = seq(-450, 150, by = 100), main = "Activator-All Evidence",
      xlab = "Distance", ylab = "Number of TFBS", col = "ForestGreen", cex.axis = 0.8, las = 2, xaxt="n")
  axis(1, at=seq(-400, 100, by = 100), labels=seq(-400, 100, by = 100), cex=0.8, las = 2)}

{hist(bs.ompr.repressor$DIST.TSS, breaks = seq(-450, 150, by = 100), main = "Repressor-All  Evidence",
      xlab = "Distance", ylab = "Number of TFBS", col = "Maroon", cex.axis = 0.8, las = 2, xaxt="n")
  axis(1, at=seq(-400, 100, by = 100), labels=seq(-400, 100, by = 100), cex=0.8, las = 2)}

dev.off()

#### MINIMUM AN MAXIMUM BINDING SITES DISTANCES 
activator.min <- quantile(bs.ompr.activator$DIST.TSS, probs = seq(0, 1, 0.1))[2]
activator.max <- quantile(bs.ompr.activator$DIST.TSS, probs = seq(0, 1, 0.1))[10]

represor.min <- quantile(bs.ompr.repressor$DIST.TSS, probs = seq(0, 1, 0.1), na.rm = TRUE)[2]
represor.max <- quantile(bs.ompr.repressor$DIST.TSS, probs = seq(0, 1, 0.1), na.rm = TRUE)[10]


#### YOU HAVE:
#### - A SET OF BINDING SITES ASSOCIATED WITH GENES 
#### - A SET OF EFFECTS UPON GENES
#### - A SET OF TSS??S ASSOCITAED WITH GENES
#### YOU WANT TO IDENTIFY:
#### THE MOST PROBABLE RELATIONSHIP BS-TSS

### Associate an effect with each target gene 

peak <- read.csv(effect.file, header=F, stringsAsFactors = F)
names(peak) <- c("Run", "Sample", "TF", "Target", "LogFoldFPKM",  "LogFoldTPM",	"FPKM",	"Counts","WildTypeFPKM","WildTypeTPM","WTCounts")

peak$effect[peak$LogFoldFPKM > 0] <- "+"
peak$effect[peak$LogFoldFPKM < 0] <- "-"


### Associate each binding site with its respective target gene 

all <- read.csv(target.file, header=F, stringsAsFactors = F)
names(all) <- c("Exp", "sample", "TF", "type", "ID1", "Start", "Stop", "PeakPos", "Height", "No1", "No2", "Shift", "Bnumber", "Gene","Dist")

##### Read TSS information (RegulonDB). Associate each TSS with a gene

tss <- read.table(tss.file, header=F, stringsAsFactors = F, sep = "\t")
names(tss) <- c("ID", "Promoter.Name", "TSS", "Sigma", "Strand", "GI", "Gene", "PosLeft", "PosRigth", "Evidence", "Bnumber", "PromSeq", "DistPromGene" )
tss <- subset(tss, !is.na(TSS) )

##### Calculate distance TSS-BS. For every gene (separate activated or repressed genes):

##   Calculate number of TSSs per gene

gene.repressed <- peak$Target[which(peak$effect == "-")]
gene.activated <- peak$Target[which(peak$effect == "+")]

no.tss.repressed <- sapply(gene.repressed, function(x,tss){ nrow(subset(tss, Gene == x))}, tss = tss, simplify = T)
no.tss.activated <- sapply(gene.activated, function(x,tss){ nrow(subset(tss, Gene == x))}, tss = tss, simplify = T)

# For all genes with at least one associated TSS:  
#  * Look all BS ChIP-Seq associated with that gene
# * Calculate distance between each BS and each TSS

calculate_distance <- function(x, tss, bs){
  bs.gene <- subset(bs, Gene == x)
  tss.gene <- subset(tss, Gene == x)
  if (nrow(tss.gene) > 0){
    distance.all <- sapply(bs.gene$PeakPos, function(x,tss){
      if (tss.gene$Strand[1] == "forward"){
        distance <- tss.gene$TSS - x  
      }else{
        distance <- x - tss.gene$TSS 
      }
      distance
    } ,simplify = T, tss = tss)
  }else{
    return (NA)
  }
}

TSS.activated <- sapply(gene.activated[no.tss.activated > 0], calculate_distance, tss = tss, bs = all, simplify =T)
TSS.repressed <- sapply(gene.repressed[no.tss.repressed > 0], calculate_distance, tss = tss, bs = all, simplify =T)


## DISTANCE BETWEEN EACH BS AND EACH TSS  

## - Una columna para cada BS  
## - Una fila para cada TSS 

#### PRINT ALL POSSIBLE INTERACTIONS
dir.create("TSS")
dir.create("TSS/Activated")
dir.create("TSS/Repressed")


print_list <- function(x, names, list, type, dir){
  write.table(list[[x]], file = paste(dir, paste(names[x], type, sep="_"), sep="/"), col.names = F, row.names = F, quote = F)
}

TSS.repressed.names <- names(TSS.repressed)
lapply(seq(1, length(TSS.repressed)), print_list, names=TSS.repressed.names, list=TSS.repressed, type="repressed.txt", dir="TSS/Repressed")

TSS.activated.names <- names(TSS.activated)
lapply(seq(1, length(TSS.activated)), print_list, names=TSS.activated.names, list=TSS.activated, type="activated.txt", dir="TSS/Activated")


##### Filter BS-TSS interactions

TSS.repressed.no <- sapply(TSS.repressed, function(x){length(as.vector(x))}, simplify=T)
write.table(TSS.repressed.no, "TSS.repressed.all", quote = F, col.names = F, row.names = T, sep="\t")

pass.repressed <- lapply(TSS.repressed, function(x, min, max){ x[x>min & x<max]}, min = -250 , max = 50 )

no.repressed <- sapply(pass.repressed, length, simplify=T)
write.table(no.repressed, "TSS.repressed.pass", quote = F, col.names = F, row.names = T, sep="\t")


TSS.activated.no <- sapply(TSS.activated, function(x){length(as.vector(x))}, simplify=T)
write.table(TSS.activated.no, "TSS.activated.all", quote = F, col.names = F, row.names = T, sep="\t")

pass.activated <- lapply(TSS.activated, function(x, min, max){ x[x>min & x<max]}, min = -250 , max = -50 )

no.activated <- sapply(pass.activated, length, simplify=T)
write.table(no.activated, "TSS.activated.pass", quote = F, col.names = F, row.names = T, sep="\t")


#### PRINT FILTER BY DISTANCE POSSIBLE INTERACTIONS

dir.create("TSS-PASS")
dir.create("TSS-PASS/Activated")
dir.create("TSS-PASS/Repressed")


pass.repressed.names <- names(pass.repressed)
lapply(seq(1, length(pass.repressed)), print_list, names=pass.repressed.names, list=pass.repressed, type="repressed.txt", dir="TSS-PASS/Repressed")

pass.activated.names <- names(pass.activated)
lapply(seq(1, length(pass.activated)), print_list, names=pass.activated.names, list=pass.activated, type="activated.txt", dir="TSS-PASS/Activated")


#### YOU HAVE:
#### - A SET OF BINDING SITES ASSOCIATED WITH GENES 
#### - A SET OF EFFECTS UPON GENES
#### - NO TSS??S ASSOCITAED WITH SOME GENES
#### YOU WANT TO IDENTIFY:
#### ASSOCIATE HT TSS??S WITH THOS GENES
#### THE MOST PROBABLE RELATIONSHIP BS:HT-TSS

 

#### Identify most probable sets BS-TSSs for genes with unknown TSS

##### Identify genes with no associated TSSs  

activated.noTSS <- names(no.tss.activated[no.tss.activated == 0])
repressed.noTSS <- names(no.tss.repressed[no.tss.repressed == 0])


##### Identify TF Binding Sites associated to those genes  

activated.noTSS.bs <- subset(all, Gene %in% activated.noTSS)
repressed.noTSS.bs <- subset(all, Gene %in% repressed.noTSS)

##### Look for TSSs (only sense TSSs) in HT data associated to each of the genes with no associated TSS

M63 <- read.table(ht.tss, stringsAsFactors = F, header=F, sep='\t', quote="")
names(M63) <- c("TSSPosition", "RPKM", "Promoter", "Strand", "RelPos", "Gene", "Bnumber", "LeftGene", "RigthGene", "Orientation", "TSSClass", "Enrichment", "evidence")

## These lines are specific to Gisella Storz data, keep only sense TSS
type.sense <- c("intragenic/sense", "upstream/sense")
sense <- subset(M63, Orientation %in% type.sense)


##### Count the number of possible associated TSSs in Storz data

tss.no.activated <- sapply(activated.noTSS, function(x,sense){ nrow(subset(sense, Gene == x))}, sense = sense, simplify = T)
tss.no.repressed <- sapply(repressed.noTSS, function(x,sense){ nrow(subset(sense, Gene == x))}, sense = sense, simplify = T)

##### For all genes with at least one associated TSS (from Storz data):  
## * Look all BS ChIP-Seq associated with that gene
## * Calculate distance between each BS and each TSS


TSS_BS_distance<- function(x, tss, bs){
  bs.gene <- subset(bs, Gene == x)
  tss.gene <- subset(tss, Gene == x)
  if (nrow(tss.gene) > 0){
    distance.all <- sapply(bs.gene$PeakPos, function(x,tss.gene){
    if (tss.gene$Strand[1] == "-"){
      distance <- tss.gene$TSSPosition - x  
    }else{
      distance <- x - tss.gene$TSSPosition
    }
    distance
  }, tss.gene = tss.gene, simplify=T)
  return(distance.all)
  }else{
    return (NA)
  }
}

### DISTANCE BETWEEN EACH BS AND EACH TSS  
# - Una columna para cada BS  
# - Una fila para cada TSS  

TSS_BS.dist.activated <- sapply(activated.noTSS[tss.no.activated > 0], TSS_BS_distance, tss=sense, bs=all, simplify = T)
TSS_BS.dist.repressed <- sapply(repressed.noTSS[tss.no.repressed > 0], TSS_BS_distance, tss=sense, bs=all, simplify = T)

#### PRINT ALL POSIBLE INTERACTIONS

dir.create("HT-TSS")
dir.create("HT-TSS/Activated")
dir.create("HT-TSS/Repressed")


TSS_BS.dist.activated.names <- names(TSS_BS.dist.activated)
lapply(seq(1, length(TSS_BS.dist.activated)), print_list, names=TSS_BS.dist.activated.names, list=TSS_BS.dist.activated, type="repressed.txt", dir="HT-TSS/Repressed")

TSS_BS.dist.repressed.names <- names(TSS_BS.dist.repressed)
lapply(seq(1, length(TSS_BS.dist.repressed)), print_list, names=TSS_BS.dist.repressed.names, list=TSS_BS.dist.repressed, type="activated.txt", dir="HT-TSS/Activated")


##### Filter BS-TSS interactions
repressor.all <- sapply(TSS_BS.dist.repressed, function(x){length(as.vector(x))}, simplify=T)
write.table(repressor.all, "HT-TSS.repressed.all", quote = F, col.names = F, row.names = T, sep="\t")

pass.repressor <- lapply(TSS_BS.dist.repressed, function(x, min, max){ x[x>min & x<max]}, min = -250 , max = 50 )
no.repressor <- sapply(pass.repressor, length, simplify=T)
write.table(no.repressor, "HT-TSS.repressed.pass", quote = F, col.names = F, row.names = T, sep="\t")


activator.all <- sapply(TSS_BS.dist.activated, function(x){length(as.vector(x))}, simplify=T)
write.table(activator.all, "HT-TSS.activated.all", quote = F, col.names = F, row.names = T, sep="\t")

pass.activated <- lapply(TSS_BS.dist.activated, function(x, min, max){ x[x>min & x<max]}, min = -250 , max = -50 )
no.activated <- sapply(pass.activated, length, simplify=T)
write.table(no.activated, "HT-TSS.activated.pass", quote = F, col.names = F, row.names = T, sep="\t")


#### PRINT FILTER BY DISTANCE INTERACTIONS

dir.create("HT-TSS-PASS")
dir.create("HT-TSS-PASS/Activated")
dir.create("HT-TSS-PASS/Repressed")

pass.activated.names <- names(pass.activated)
lapply(seq(1, length(pass.activated)), print_list, names=pass.activated.names, list=pass.activated, type="repressed.txt", dir="HT-TSS-PASS/Repressed")

pass.repressor.names <- names(pass.repressor)
lapply(seq(1, length(pass.repressor)), print_list, names=pass.repressor.names, list=pass.repressor, type="activated.txt", dir="HT-TSS-PASS/Activated")

