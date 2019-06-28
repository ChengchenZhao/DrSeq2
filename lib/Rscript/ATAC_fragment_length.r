a<-commandArgs(T)

infile <- a[1]
outname <- a[2]

fragment_length <- read.table(infile,header=T)
fl_x <- as.numeric(fragment_length[,1])
fl_y <- as.numeric(fragment_length[,2])


pdf(paste(outname,"_Figure5_fragment_length_distribution.pdf",sep=""),height=6,width=6)
plot(fl_x[order(fl_x)],fl_y[order(fl_x)]/sum(fl_y),type="l",xlim=c(0,1200),lwd=2,col="#CE0013",ylim=c(0,max(fl_y/sum(fl_y),fl_y/sum(fl_y))),xlab="fragment_length",ylab="Count number/total number")
dev.off()