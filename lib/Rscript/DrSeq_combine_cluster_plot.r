library(cluster)
a<-commandArgs(T)

cluster_mat <- read.table(a[1],row.names=1,header=T)
outname <- a[2]
use_cluster_mat <- cluster_mat[which(cluster_mat[,3]>0),]
unuse_cluster_mat <- cluster_mat[which(cluster_mat[,3]==0),]
silhouette_evaluate <- function(clusmat){
    C <- clusmat[,3]
    D <- dist(clusmat[,1:2])
    return(silhouette(C,D))
}
silScore <- silhouette_evaluate(use_cluster_mat)

if (is.na(silScore)){
    use_cluster_mat_sil <- cbind(use_cluster_mat,rep('NA',nrow(use_cluster_mat)))
}else{
    use_cluster_mat_sil <- cbind(use_cluster_mat,silScore[,'sil_width'])
    pdf(file=paste(outname,"_Figure12_silhouetteScore.pdf",sep=""))
    plot(silScore,main="")
    dev.off()
}
colnames(use_cluster_mat_sil)[4] <- 'silhouette' 
if(nrow(unuse_cluster_mat) == 0){
    cluster_mat_sil  <- use_cluster_mat_sil[rownames(cluster_mat),]
}else{
    unuse_cluster_mat_sil <- cbind(unuse_cluster_mat,rep('NA',nrow(unuse_cluster_mat)))
    colnames(unuse_cluster_mat_sil)[4] <- 'silhouette' 
    cluster_mat_sil <- rbind(use_cluster_mat_sil,unuse_cluster_mat_sil)[rownames(cluster_mat),]
}