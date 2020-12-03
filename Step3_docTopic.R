num_normal <- 49  # Number of normal samples: breast = 112, lung = 49
num_tumor  <- 49  # Number of tumor samples:  breast = 117, lung = 49
num_bins   <- 5   # number of bins
k          <- 3   # number of topic
topicProbfile <- paste('2_k', k, 'topicProb.txt', sep='')    # input


###### read input file
doctopicf <- as.matrix(read.table(file = topicProbfile, sep=',', header=T, fill=TRUE))
dim(doctopicf)   #fill=TRUE for filling blank value; number of column depends on first row


normaltopic <- doctopicf[1:num_normal,]
cancertopic <- doctopicf[(num_normal+1):(num_normal+num_tumor),]
(normalmean <- colMeans(normaltopic))
(cancermean <- colMeans(cancertopic))


###### generate bar plot of grouped samples (normal vs tumor)
par(mfrow=c(1,2))
barplot(normalmean, ylim=c(0,1), main='Normal Sample', xlab="Topic", ylab='Probability', names.arg=paste('Topic', 1:k))
barplot(cancermean, ylim=c(0,1), main='Cancer Sample', xlab="Topic", ylab='Probability', names.arg=paste('Topic', 1:k))


###### generate heatmap
heatmap(doctopicf, col=heat.colors(75, alpha=1), RowSideColors=c(rep('blue', num_normal), rep('red', num_tumor)), labCol=paste('Topic', 1:k))

 
###### generate bar plot of invididual samples
cols <- c('blue',"red")
posstat <- rep(T, num_normal+num_tumor)
posstat[1:num_normal] <- F
op <- par(mfrow = c(k,1), oma = c(5,4,0,0) + 0.1, mar = c(0,0,2,0) + 0.1)
for (i in 1:k){
  barplot(doctopicf[,i], ylim=c(0,1), xlab=paste('Topic',i,sep=' '), space=0, col=cols[posstat+1], border=cols[posstat+1])
}

