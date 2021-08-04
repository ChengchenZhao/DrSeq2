sys_argv <- commandArgs(T)
#sys_argv<-c("G:/A-Dr.Seq-Project/readsnumber/test_ATAC.readsnumber_sort.txt","test")
cccol <- c("#CE0013","#16557A","#C7A609","#87C232","#008792","#A14C94","#15A08C","#8B7E75","#1E7CAF","#EA425F","#46489A","#E50033","#0F231F","#1187CD")
cccol50 <- c("#CE001350","#16557A50","#C7A60950","#87C23250","#00879250","#A14C9450","#15A08C50","#8B7E7550","#1E7CAF50","#EA425F50","#46489A50","#E5003350","#0F231F50","#1187CD50")

data<-read.table(sys_argv[1],header = F ,sep = "\t")
out_name<-sys_argv[2]
colnames(data)<-c("barcode","number")
data<-data[which(data$number>10),]
barcode<-data[,1]
number<-data[,2]

id<-seq(length(barcode))
reads_num<- number/max(number)
data2<-as.data.frame(cbind(id,reads_num))
k<-log10(reads_num[2:length(reads_num)])-log10(reads_num[1:length(reads_num)-1])
k[1:100]<-max(k)
sk<-smooth.spline(k)$y
n<-which(sk==min(sk))

png(paste(out_name,"_CellReadsDistribution.png",sep=""))
plot(seq(length(reads_num)),(reads_num),type="l",col=cccol[2],lwd=2,xlab="Cell barcode id",ylab="Reads number(log10)",main=out_name);box(lwd=2)
abline(v=n,lty=2,lwd=2)
text(n,max((reads_num))*0.9,n,col=cccol[1])
dev.off()
data[1:n,1]
write.table(data[1:n,1],paste(out_name,"_cells_info_barcodes.txt",sep=""),row.names = F,col.names = FALSE, quote = FALSE)
