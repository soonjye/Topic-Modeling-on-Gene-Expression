inputfile   <- 'TCGA_Cancer/TCGA_LungCancer.csv'    # datafile
num_normal  <- 49    # Number of normal samples: breast = 112, lung = 49
num_tumor   <- 49    # Number of tumor samples:  breast = 117, lung = 49
num_bins    <- 5     # number of bins
outputfile1 <- "1_exp2rank.csv"
outputfile2 <- "2_rank2word.csv"



################## Read data files
fpkmdata <- read.table(file=inputfile, sep=',', row.names=1)
## Change inputfile to BreastCancer or LungCancer

(totalsample <- dim(fpkmdata)[2])
(totalgene <- dim(fpkmdata)[1])
(portion <- floor(totalgene/num_bins))



################## Divide the genes into different bins (j) based on expression level
for (i in 1:totalsample) {
  if (i == 1) {
    cat('Converting Expression level to Bins: sample ', i)
  } else {
    cat(' ', i)
  }
  
  eachprf <- fpkmdata[,i]
  names(eachprf) <- rownames(fpkmdata)
  eachprf <- sort(eachprf)
  rankprf <- eachprf
  k <- 1
  
  #### divides the genes to 20 bins
  for (j in 1:num_bins) {
    rankprf[(portion*(j-1)+k):(portion*(j-1)+portion)] <- j
    k<-1
    
    #### the genes may not be evenly divide into j bins, thus, some bins might have more (range = [0, j]) genes
    while (eachprf[(portion*(j-1)+portion+k)] == eachprf[(portion*(j-1)+portion)]){
      rankprf[(portion*(j-1)+portion+k)] <- j
      k <- k+1
    }
  }
  rankprf[(totalgene-(totalgene%%5)):totalgene] <- num_bins
  
  
  if (i==1) {
    allrank <- rankprf
  } else {
    allrank <- merge(allrank,rankprf,by='row.names',all=T)
    rownames(allrank) <- allrank[,1]
    allrank <- allrank[,-1]
    colnames(allrank) <- 1:i
  }
}

# output: converted expression level into bins
write.table(allrank, file=outputfile1, sep=',', row.names=T, col.names=F)

#heatmap(as.matrix(allrank), ColSideColors=c(rep('blue',num_normal), rep('red',num_tumor)))





################## Convert into word "gene-bin"
word <- matrix(NA, ncol=totalgene, nrow=totalsample)
label <- c(rep('normal',num_normal), rep('tumor', num_tumor))

for (i in 1:totalsample) {
  if (i == 1) {
    cat('Converting Gene & Bin to Word: sample', i)
  } else {
    cat(' ', i)
  }
  
  for (j in 1:totalgene) {
    word[i,j] <- paste(rownames(allrank)[j] , '-', allrank[j,i] , sep='')
  }
}

write.table(word, file=outputfile2, sep=',', col.names=F, row.names=F)
