a<-commandArgs(T)
# a <- c("A01.peak,A02.peak,A03.peak,A04.peak,A05.peak,A06.peak,A07.peak,A08.peak,A09.peak,A10.peak,A11.peak,A12.peak,A13.peak,A14.peak,A15.peak,A16.peak,A17.peak,A18.peak,A19.peak,A20.peak,A21.peak,A22.peak,A23.peak,A24.peak,A25.peak,A26.peak,A27.peak,A28.peak,A29.peak,A30.peak,A31.peak,A32.peak,A33.peak,A34.peak,A35.peak,A36.peak,A37.peak,A38.peak,A39.peak,A40.peak,A41.peak,A42.peak,A43.peak,A44.peak,A45.peak,A46.peak,A47.peak,A48.peak,A49.peak,A50.peak,A51.peak,A52.peak,A53.peak,A54.peak,A55.peak,A56.peak,A57.peak,A58.peak,A59.peak,A60.peak,A61.peak,A62.peak,A63.peak,A64.peak,A65.peak,A66.peak,A67.peak,A68.peak,A69.peak,A70.peak,A71.peak,A72.peak,A73.peak,A74.peak,A75.peak,A76.peak,A77.peak,A78.peak,A79.peak,A80.peak,B01.peak,B02.peak,B03.peak,B04.peak,B05.peak,B06.peak,B07.peak,B08.peak,B09.peak,B10.peak,B11.peak,B12.peak,B13.peak,B14.peak,B15.peak,B16.peak,B17.peak,B18.peak,B19.peak,B20.peak,B21.peak,B22.peak,B23.peak,B24.peak,B25.peak,B26.peak,B27.peak,B28.peak,B29.peak,B30.peak,B31.peak,B32.peak,B33.peak,B34.peak,B35.peak,B36.peak,B37.peak,B38.peak,B39.peak,B40.peak,B41.peak,B42.peak,B43.peak,B44.peak,B45.peak,B46.peak,B47.peak,B48.peak,B49.peak,B50.peak,B51.peak,B52.peak,B53.peak,B54.peak,B55.peak,B56.peak,B57.peak,B58.peak,B59.peak,B60.peak,B61.peak,B62.peak,B63.peak,B64.peak,B65.peak,B66.peak,B67.peak,B68.peak,B69.peak,B70.peak,B71.peak,B72.peak,B73.peak,B74.peak,B75.peak,B76.peak,B77.peak,B78.peak,B79.peak,B80.peak,C01.peak,C02.peak,C03.peak,C04.peak,C05.peak,C06.peak,C07.peak,C08.peak,C09.peak,C10.peak,C11.peak,C12.peak,C13.peak,C14.peak,C15.peak,C16.peak,C17.peak,C18.peak,C19.peak,C20.peak,C21.peak,C22.peak,C23.peak,C24.peak,C25.peak,C26.peak,C27.peak,C28.peak,C29.peak,C30.peak,C31.peak,C32.peak,C33.peak,C34.peak,C35.peak,C36.peak,C37.peak,C38.peak,C39.peak,C40.peak,C41.peak,C42.peak,C43.peak,C44.peak,C45.peak,C46.peak,C47.peak,C48.peak,C49.peak,C50.peak,C51.peak,C52.peak,C53.peak,C54.peak,C55.peak,C56.peak,C57.peak,C58.peak,C59.peak,C60.peak,C61.peak,C62.peak,C63.peak,C64.peak,C65.peak,C66.peak,C67.peak,C68.peak,C69.peak,C70.peak,C71.peak,C72.peak,C73.peak,C74.peak,C75.peak,C76.peak,C77.peak,C78.peak,C79.peak,C80.peak","scATAC_test","scATAC_test_total_peaks.location","0","5","5","3")
infiles <- a[1]
outname <- a[2]
total_peak_file <- a[3]
cut_height <- as.numeric(a[4])
cell_cutoff <- as.numeric(a[5])
peak_cutoff <- as.numeric(a[6])
given_cluster_number <- as.numeric(a[7])

cells <- strsplit(infiles,",")[[1]]
peaks <- as.vector(read.table(total_peak_file)[,1])

out_matrix <- matrix(0, nrow = length(peaks), ncol = length(cells))
rownames(out_matrix) <- peaks
out_col <- c()
for (each_file in cells){
  tmp <- strsplit(each_file,"/")[[1]]
  cell_name <- strsplit(tmp[length(tmp)],"\\.")[[1]][1]
  out_col <- c(out_col,cell_name)
}
colnames(out_matrix) <- out_col

for (each_file in cells){
	tmp <- strsplit(each_file,"/")[[1]]
	cell_name <- strsplit(tmp[length(tmp)],"\\.")[[1]][1]
	if (file.info(each_file)$size > 0){
		tmp_info <- as.vector(read.table(each_file)[,1])
		out_matrix[tmp_info,cell_name] <- 1
	}
}

write.table(out_matrix,file=paste(outname,"peakMatrix.txt",sep=""),col.names=T,row.names=T,quote=F)

pdf(paste(outname,"_Figure4_peak_number_distribution.pdf",sep=""),width=6,height=6)
hist(apply(out_matrix,2,sum),border="#CE0013",col="#CE0013",breaks=100,main="",xlab="peak number")
dev.off()

RemoveUnwantedLine <- function(signal_per_cell,peak_cutoff){
	return(sum(signal_per_cell) >= peak_cutoff)
}
NormalizeSignal <- function(x){
	if (sum(x)!= 0){
		x <- x/sum(x)
	}
	else{
		x <- x
	}
}

pc_signal <- out_matrix
cp_signal <- t(pc_signal) # cell*peak

signal <- cp_signal[names(which(apply(cp_signal,1,RemoveUnwantedLine,peak_cutoff) == TRUE)),]
signal <- signal[,names(which(apply(signal,2,RemoveUnwantedLine,cell_cutoff) == TRUE))]
signal <- signal[apply(signal,1,sum)!=0,apply(signal,2,sum)!=0]
signal <- apply(signal,2,NormalizeSignal)

if (length(signal) == 0){
	print ("NONE of cells pass your cutoff,please loose your cutoff.")
}


cell_cor <- dist(cor(t(signal)))
# cell_cov <- dist(cov(t(signal)))
hc <- hclust(cell_cor)
# library(amap)
# hc <- hcluster(signal, method = "correlation")

pdf(paste(outname,"_Figure6_cell_clusting.pdf",sep=""),width=10,height=4)
plot(hc, labels = FALSE, main = "hclust based on peak", xlab = "cells",)
dev.off()

if (cut_height == 0){
    i <- 0.8
    while (i <= 1){
      cut_height <- quantile(hc$height,i)
      mycl <- cutree(hc, h=cut_height)
      cluster_number <- length(unique(mycl))
      if (cluster_number <= given_cluster_number){
        break
      }
      else{
        i <- i + 0.001
      }
    }
}else{
    mycl <- cutree(hc, h=cut_height)
}


silhouette_evaluate <- function(clusmat,cluster){
    library(cluster)
    C <- cluster
    D <- dist(clusmat)
    return(silhouette(C,D))
}
silScore <- silhouette_evaluate(signal,mycl)
use_cluster_mat <- cbind(mycl)
if (is.na(silScore)){
    use_cluster_mat_sil <- cbind(use_cluster_mat,rep('NA',nrow(use_cluster_mat)))
}else{
    use_cluster_mat_sil <- cbind(use_cluster_mat,silScore[,'sil_width'])
    colnames(use_cluster_mat_sil) <- c("cluster","sil_width")
    pdf(file=paste(outname,"_Figure7_silhouetteScore.pdf",sep=""))
    plot(silScore,main="")
    dev.off()
}
write.table(use_cluster_mat_sil,file = paste(outname,"_cluster_with_silhouette_score.txt",sep=""),quote=FALSE)

# k <- 10
# set.seed(1)
# km <- Kmeans(t(signal),k,method = "correlation")

clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
png(file = paste(outname,"_Figure8_heatmap.png",sep=""), width = 600, height = 600);
par(oma=c(0.5,0.5,0.5,0.5),mar=c(2,2,2,2))
layout(matrix(c(rep(1,2),rep(2,1),rep(3,8),rep(4,1)),ncol=1,nrow=12))
plot(hc, labels = FALSE, main = "hclust based on peak", xlab = "cells")
image(seq(length(mycl)),1,matrix(data=mycl[hc$order], nrow=length(mycl),ncol=1),col=rainbow(length(unique(mycl))), xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n")
allPlotMatrix <- t(signal)[,hc$order]
# allPlotMatrix <- c()
# for (each in seq(k)){
# 	plotMatrix <- t(signal)[names(which(km$cluster==each)),hc$order]
# 	allPlotMatrix <- rbind(allPlotMatrix,plotMatrix)
# }
all_exp <- c(as.matrix(allPlotMatrix))
zmax <- max(all_exp)
zmin <- min(all_exp)
ColorRamp <- colorRampPalette(c("#F2FAB0","red"), bias=1)(10000)   #color list
ColorLevels <- seq(to=zmax,from=zmin, length=10000)   #number sequencemodGenes <- names(which(km$cluster==each))
image(1:ncol(allPlotMatrix), 1:nrow(allPlotMatrix), t(allPlotMatrix), xaxt="n", yaxt="n", col=ColorRamp, xlab="", ylab="")
image(ColorLevels,1,matrix(data=ColorLevels, nrow=length(ColorLevels),ncol=1),col=ColorRamp, xlab="",ylab="",cex.axis=2,xaxt="n",yaxt="n",useRaster=T)
axis(side=1,c(zmin,round((zmax-zmin)/2,1),zmax),labels=c(round(zmin,2),round((zmax-zmin)/2,1),round(zmax,1)))
dev.off()

for (i in seq(length(unique(mycl)))){
  peak_cover <- apply(signal[names(mycl[mycl==i]),],2,sum)
  other_cell_peak_cover <- apply(signal[names(mycl[mycl!=i]),],2,sum)
  tmp_specific_peak <- c()
  for (j in seq(length(peak_cover))){
    if (peak_cover[j] > 0 & other_cell_peak_cover[j] == 0) {
      tmp_specific_peak <- c(tmp_specific_peak,chartr("_", "\t",names(peak_cover)[j]))
    }
  }
  write.table(cbind(names(mycl[mycl==i])),file=paste(outname,"_cluster",i,"_cells.txt",sep=""),quote=F,col.names=F,row.names=F)
	write.table(cbind(tmp_specific_peak),file=paste(outname,"_cluster",i,"_specific.peak.bed",sep=""),quote=F,col.names=F,row.names=F)
}
