k        <- 3
num_bins <- 5
topicPrffile   <- '3_topicPrf.csv'         # input
dysreg_file    <- '4_dysreg_genes.csv'     # output
genelist1_file <- '5_versus2v1.csv'        # output
genelist2_file <- '5_versus3v1.csv'        # output




################################   Examine Grouped Normal vs Tumor Profile  ############################################################
topicprf <- read.table(file=topicPrffile, sep=',', header=T, row.names=1)
exceed <- which(abs(topicprf$prf_diff) > (num_bins*0.25))   # threshold is 1/4 of bins (e.g. threshold = 5 if bins = 20)
(genelistno <- length(exceed))
genelist <- topicprf[exceed,]
genelist <- cbind(genelist, genesymbol=rep(NA, genelistno))

library('biomaRt')
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
conversym <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"), values= rownames(genelist), mart= mart)
write.table(conversym, file=dysreg_file, sep=',', row.names=F, col.names=F)








################################   Examine individual Topic Profile  ############################################################
#Examine difference between topics
#E.g. tumor sample has highest prob in Topic_2. Thus examine difference of T2 vs with base reference T1
versus2v1 <- topicprf[, 2] - topicprf[, 1]
mostpath  <- which(abs(versus2v1) > (num_bins*0.25))        # threshold is 1/4 of bins (e.g. threshold = 5 if bins = 20)
(mostlistno <- length(mostpath))

mostpath <- data.frame(ENSG=rownames(topicprf)[mostpath], topicprf[mostpath, c(2,1)], topic_diff = versus2v1[mostpath])
conversym <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"), values=rownames(mostpath), mart=mart)
mostpath <- merge(mostpath, conversym, by.x = 'ENSG', by.y='ensembl_gene_id')
mostpath$overlap <- NA

write.table(mostpath, file=genelist1_file, sep=',', row.names=F, col.names=T)
latestage <- as.character(mostpath[,1])



#E.g. normal sample has highest prob in Topic_3. Thus examine difference of T3 vs with base reference T1
versus3v1 <- topicprf[,3]-topicprf[,1]
rarepath  <- which(abs(versus3v1) > (num_bins*0.25))        # threshold is 1/4 of bins (e.g. threshold = 5 if bins = 20)
(rarelistno <- length(rarepath))

rarepath <- data.frame(ENSG=rownames(topicprf)[rarepath], topicprf[rarepath, c(3,1)], topic_diff = versus3v1[rarepath])
conversym <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id", "hgnc_symbol"), values= rownames(rarepath), mart=mart)
rarepath <- merge(rarepath, conversym, by.x = 'ENSG', by.y='ensembl_gene_id')
rarepath$overlap <- NA

write.table(rarepath, file=genelist2_file, sep=',', row.names=F, col.names=T)
prestage <- as.character(rarepath[,1])



################################  Examine genelist of versus2v1 & versus3v1
overlap <- intersect(prestage,latestage)
(intersectno <- length(overlap))          #17 genes are overlap

for (i in 1:mostlistno){
  if (mostpath[i,1] %in% overlap)    mostpath$overlap[i] <- 1    else    mostpath$overlap[i] <- 0
}
for (i in 1:rarelistno){
  if (rarepath[i,1] %in% overlap)    rarepath$overlap[i] <- 1    else    rarepath$overlap[i] <- 0
}


## generate Venn Diagram for comparison of two genelists
library(VennDiagram)
draw.pairwise.venn(area1=length(prestage), area2=length(latestage), cross.area=intersectno, cat.pos=c(0,0), cat.dist=c(-0.02,0.02), scaled=F, fill = c("blue", "yellow"))
#breast cancer = 1168, 593, 559
#lung cancer = 482, 106, 42




