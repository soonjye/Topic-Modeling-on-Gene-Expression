num_normal <- 49  # Number of normal samples: breast = 112, lung = 49
num_tumor  <- 49  # Number of tumor samples:  breast = 117, lung = 49
num_bins   <- 5   # number of bins
k          <- 3   # number of topic
termProbfile  <- paste('2_k', k, 'termProb.txt', sep='')    # input
topicProbfile <- paste('2_k', k, 'topicProb.txt', sep='')   # input
topicPrffile  <- '3_topicPrf.csv'                         # output


######## impute bin of each gene in each topic (AKA topic profile)
termProb     <- as.matrix(read.table(file = termProbfile, sep=',', header=T))
dim(termProb)   #fill=TRUE for filling blank value; number of column depends on first row

genes <- unique(substring(termProb[,1], 1, 15))

topicprf <- data.frame()
for (t in 1:k) {
  tokens <- termProb[, t*2-1]
  tokeng <- substring(tokens, 1, 15)
  tokenb <- as.numeric(gsub("^.*-", "", tokens))
  tokenp <- as.numeric(termProb[, t*2])
  
  word_topic <- matrix(NA, nrow=length(genes), ncol=num_bins)
  rownames(word_topic) <- genes
  for (i in 1:length(tokens)) {
    word_topic[tokeng[i], tokenb[i]] <- tokenp[i]
  }

  rowm <- apply(word_topic, 1, function(x) sum(x*1:num_bins, na.rm=T) / sum(x, na.rm=T))
  # formula: sum(prob*rank) / sum(prob)
  if (nrow(topicprf) == 0) {
    topicprf <- data.frame(t1=rowm)
  } else {
    topicprf <- cbind(topicprf, data.frame(rowm))
  }
}

colnames(topicprf) <- paste("topic", 1:k, sep='_')




######## impute difference of bin between grouped normal and tumor samples
doctopic <- as.matrix(read.table(file=topicProbfile, sep=',', header=T))

(normalmean <- colMeans(doctopic[1:num_normal,]))
(cancermean <- colMeans(doctopic[(num_normal+1):(num_normal+num_tumor),]))

topicprf$normal_prf <- apply(topicprf, 1, function(x) sum(normalmean * x[1:k]))
topicprf$tumor_prf  <- apply(topicprf, 1, function(x) sum(cancermean * x[1:k]))
topicprf$prf_diff   <- apply(topicprf, 1, function(x) x['normal_prf'] - x['tumor_prf'])
write.table(topicprf, file=topicPrffile, sep=',', row.names=T, col.names=T)



