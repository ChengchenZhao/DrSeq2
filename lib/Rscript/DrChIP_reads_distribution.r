a <- commandArgs(T)
infile <- a[1]
outname <- a[2]

# library(mclust)


reads_number_per_cell <- read.table(infile)
cell_barcode_id <- as.numeric(reads_number_per_cell[,2])
reads_number <- as.numeric(reads_number_per_cell[,1])

# x <- log10(reads_number)
# m <- densityMclust(x,model='V')
# att <- m$parameters
# pdf(paste(outname,'_Figure7_readsDistribution.pdf',sep=""),width=6,height=6)
# #hist(x, freq = FALSE);
# plot(density(x),col='red',main='Distribution of reads number per cell',xlab='Reads Number',xaxt='n',xlim=c(0,6),lwd=2)
# axis(1,at=c(0:5),labels=c('0','10','10^2','10^3','10^4','10^5'))
# cutoff <- c()
# colors <- c('lightblue','yellow')
# for(i in 1:att$variance$G)
# {
# curve(att$pro[i] * dnorm(x, att$mean[i], sqrt(att$variance$sigmasq[i])),n = 300, add = TRUE,col=colors[i],lwd=2,lty=2);
# cutoff <- c(cutoff,qnorm(0.99, att$mean[i], sqrt(att$variance$sigmasq[i])))
# # abline(v=qnorm(0.99, att$mean[i], sqrt(att$variance$sigmasq[i])))
# }
# legend("topright",c("Background","individual cell"),col=colors[1:2],lwd=3,bty="n")
# # print(sum(((reads_number <= 10**cutoff[2]) == TRUE) & ((reads_number >= 10**cutoff[1]) == TRUE)))
# # cell_barcode_id[which(((reads_number <= 10**cutoff[2]) == TRUE) & ((reads_number >= 10**cutoff[1]) == TRUE))]
# dev.off()

x <- log10(reads_number)
pdf(paste(outname,'_Figure7_readsDistribution.pdf',sep=""),width=6,height=6)
h <- hist(x, plot=FALSE,breaks=100)
h$counts=h$counts/sum(h$counts)
# plot(h,col="grey", main='Distribution of reads number per cell', xlab='Reads Number',xlim=c(0,6),ylim=c(0,1))
# lines(density(x), col="green")
plot(density(x), col="red", main='Distribution of reads number per cell', xlab='Reads Number',xlim=c(0,6),ylim=c(0,1),lwd=2)

# hist(x, freq = TRUE, breaks=200);
# plot(density(x),col='red',main='Distribution of reads number per cell',xlab='Reads Number',xaxt='n',xlim=c(0,6),lwd=2)
# # points(density(x),col='red',main='Distribution of reads number per cell',xlab='Reads Number',xaxt='n',xlim=c(0,6),lwd=2)
# axis(1,at=c(0:5),labels=c('0','10','10^2','10^3','10^4','10^5'))
dev.off()