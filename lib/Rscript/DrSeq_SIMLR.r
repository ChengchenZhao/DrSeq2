a<-commandArgs(T)

inmatrix <- a[1]
outname <- a[2]
hvZ <- as.numeric(a[3])
RDnumber <- as.numeric(a[4])
maxKnum <- as.numeric(a[5])
cortableY <- as.numeric(a[6])
clustering_method <- as.numeric(a[7])
custom_k <- as.numeric(a[8])
custom_d <- as.numeric(a[9])
R_dir <- a[10]
  
# required external packages for SIMLR
library(Matrix)
library(parallel)
# load the igraph package to compute the NMI
library(igraph)
# load the palettes for the plots
library(grDevices)
# load the SIMLR R package
source(paste(R_dir,"SIMLR.R",sep=""))
source(paste(R_dir,"compute.multiple.kernel.R",sep=""))
source(paste(R_dir,"network.diffusion.R",sep=""))
source(paste(R_dir,"utils.simlr.R",sep=""))
source(paste(R_dir,"tsne.R",sep=""))
dyn.load(paste(R_dir,"projsplx_R.so",sep=""))
set.seed(11111)

### internal function
stable_gap_decideK <- function(indata,SD,maxK,Gapplot){

    set.seed(SD)
    b <- clusGap(indata,kmeans,maxK+20)
    GapK <- b$Tab[,3]

    Knum = 0
    GapCUTOFF <- mean(GapK[maxK:(maxK+20)]) - 5*sd(GapK[maxK:(maxK+20)])
    for(i in seq(maxK-1)){
        if(GapK[i] > GapCUTOFF & GapK[i] > GapK[i+1]){
            Knum = i
            break}
    }
    pdf(file=Gapplot)
    plot(b)
    abline(v=Knum,lwd=2,col="blue")
    abline(h=GapCUTOFF,lwd=2,col="darkgreen")
    legend("bottomright",legend=c("selected K value","cutoff for stable gap score"),col=c("blue","darkgreen"),bty="n",lwd=2)
    dev.off()

    return(Knum)
}
maxSE_decideK <- function(indata,SD,maxK,Gapplot){
    set.seed(SD)
    b <- clusGap(indata,kmeans,maxK)
    Knum <- maxSE(b$Tab[,3],b$Tab[,4],method="Tibs2001SEmax",SE.factor=0.01)
    pdf(file=Gapplot)
    plot(b)
    abline(v=Knum,lwd=2,col="blue")
    legend("topleft",legend=c("selected K value"),col=c("blue"),bty="n",lwd=2)
    dev.off()
    return(Knum)
}

getTPM <- function(indata){
	return(log10(indata*1e4/sum(indata) + 1))
}
getcoverGnum <- function(indata){
	return(length(which(indata > 0)))
}
selct_high_var_gene <- function(inputdata,highvarZ){
	Gave <- apply(inputdata,1,mean)
	Gvar <- apply(inputdata,1,var)
	Gdp <- Gvar/Gave
	Gname <- rownames(inputdata)
	each <- floor(length(Gname)/20)
	g1 <- Gname[order(Gave)][(1+each*0):(each*1)]
	g2 <- Gname[order(Gave)[(1+each*1):(each*2)]]
	g3 <- Gname[order(Gave)[(1+each*2):(each*3)]]
	g4 <- Gname[order(Gave)[(1+each*3):(each*4)]]
	g5 <- Gname[order(Gave)[(1+each*4):(each*5)]]
	g6 <- Gname[order(Gave)[(1+each*5):(each*6)]]
	g7 <- Gname[order(Gave)[(1+each*6):(each*7)]]
	g8 <- Gname[order(Gave)[(1+each*7):(each*8)]]
	g9 <- Gname[order(Gave)[(1+each*8):(each*9)]]
	g10 <- Gname[order(Gave)[(1+each*9):(each*10)]]
	g11 <- Gname[order(Gave)[(1+each*10):(each*11)]]
	g12 <- Gname[order(Gave)[(1+each*11):(each*12)]]
	g13 <- Gname[order(Gave)[(1+each*12):(each*13)]]
	g14 <- Gname[order(Gave)[(1+each*13):(each*14)]]
	g15 <- Gname[order(Gave)[(1+each*14):(each*15)]]
	g16 <- Gname[order(Gave)[(1+each*15):(each*16)]]
	g17 <- Gname[order(Gave)[(1+each*16):(each*17)]]
	g18 <- Gname[order(Gave)[(1+each*17):(each*18)]]
	g19 <- Gname[order(Gave)[(1+each*18):(each*19)]]
	g20 <- Gname[order(Gave)[(1+each*19):(length(Gname))]]
	highvargene <- c(g2[which((scale(Gdp[g2]))>highvarZ)],g3[which((scale(Gdp[g3]))>highvarZ)],g4[which((scale(Gdp[g4]))>highvarZ)],g5[which((scale(Gdp[g5]))>highvarZ)],g6[which((scale(Gdp[g6]))>highvarZ)],g7[which((scale(Gdp[g7]))>highvarZ)],g8[which((scale(Gdp[g8]))>highvarZ)],g9[which((scale(Gdp[g9]))>highvarZ)],g10[which((scale(Gdp[g10]))>highvarZ)],g11[which((scale(Gdp[g11]))>highvarZ)],g12[which((scale(Gdp[g12]))>highvarZ)],g13[which((scale(Gdp[g13]))>highvarZ)],g14[which((scale(Gdp[g14]))>highvarZ)],g15[which((scale(Gdp[g15]))>highvarZ)],g16[which((scale(Gdp[g16]))>highvarZ)],g17[which((scale(Gdp[g17]))>highvarZ)],g18[which((scale(Gdp[g18]))>highvarZ)],g19[which((scale(Gdp[g19]))>highvarZ)],g20[which((scale(Gdp[g20]))>highvarZ)])
	return(highvargene)
}

# ref functions
# givenK_kmeans <- function(indata,SD,Knum,clusterplot){
#     set.seed(SD)
#     pdf(file=clusterplot)
#     x_lim <- c(quantile(indata[,1],0.01),quantile(indata[,1],0.99))
#     y_lim <- c(quantile(indata[,2],0.01),quantile(indata[,2],0.99))
#     tmp_indata <- indata[indata[,1]<x_lim[2]&indata[,1]>x_lim[1]&indata[,2]<y_lim[2]&indata[,2]>y_lim[1],]
#     km <- kmeans(tmp_indata,Knum)
#     plot(tmp_indata,pch=16,xlab="t-SNE 1",ylab="t-SNE 2",main=paste("k-means (k=", Knum, ")",sep=""),xlim=x_lim,ylim=y_lim)
#     rain <- rainbow(length(km$size))
#     for( i in 1:length(km$size)){
#         points(tmp_indata[which(km$cluster==i),],col=rain[i],pch=16)    
#     }
#     text(km$centers,labels=seq(nrow(km$centers)))
#     dev.off()
#     cluster_result <- cbind(tmp_indata,km$cluster)
#     rownames(cluster_result) <- row.names(tmp_indata)
#     return(cluster_result)
# }
givenK_kmeans <- function(indata,SD,Knum,clusterplot){
    set.seed(SD)
    pdf(file=clusterplot)
    km <- kmeans(indata,Knum)
    plot(indata,pch=16,xlab="dimension 1",ylab="dimension 2",main=paste("k-means (k=", Knum, ")",sep=""))
    rain <- rainbow(length(km$size))
    for( i in 1:length(km$size)){
        points(indata[which(km$cluster==i),],col=rain[i],pch=16)    
    }
    text(km$centers,labels=seq(nrow(km$centers)))
    dev.off()
    cluster_result <- cbind(indata,km$cluster)
    return(cluster_result)
}
givenE_dbscan <- function(indata,EPS,clusterplot){
    pdf(file=clusterplot)
    tmp_indata <- indata
    ds <- ref_dbscan(tmp_indata,eps=EPS)
    cluster_result <- cbind(tmp_indata,ds$cluster)
    kmsize <- sort(cluster_result[,3],decreasing=T)[1]
    plot(tmp_indata,pch=16,xlab="dimension 1",ylab="dimension 2",main=paste("DBSCAN (eps=", EPS, ")",sep=""))
    rain <- rainbow(kmsize)
    for(i in 1:kmsize){
        points(cluster_result[which(cluster_result[, 3]==i),1:2],col=rain[i],pch=16)
    }
    dev.off()
    rownames(cluster_result) <- row.names(tmp_indata)
    return(cluster_result)
}
ref_dbscan <- function(data, eps, MinPts = 5, scale = FALSE, method = c("hybrid", 
    "raw", "dist"), seeds = TRUE, showplot = FALSE, countmode = NULL) 
{
    distcomb <- function(x, data) {
        data <- t(data)
        temp <- apply(x, 1, function(x) {
            sqrt(colSums((data - x)^2))
        })
        if (is.null(dim(temp))) 
            matrix(temp, nrow(x), ncol(data))
        else t(temp)
    }
    method <- match.arg(method)
    data <- as.matrix(data)
    n <- nrow(data)
    if (scale) 
        data <- scale(data)
    classn <- cv <- integer(n)
    isseed <- logical(n)
    cn <- integer(1)
    for (i in 1:n) {
        if (i %in% countmode) 
            cat("Processing point ", i, " of ", n, ".\n")
        unclass <- (1:n)[cv < 1]
        if (cv[i] == 0) {
            if (method == "dist") {
                reachables <- unclass[data[i, unclass] <= eps]
            }
            else {
                reachables <- unclass[as.vector(distcomb(data[i, 
                  , drop = FALSE], data[unclass, , drop = FALSE])) <= 
                  eps]
            }
            if (length(reachables) + classn[i] < MinPts) 
                cv[i] <- (-1)
            else {
                cn <- cn + 1
                cv[i] <- cn
                isseed[i] <- TRUE
                reachables <- setdiff(reachables, i)
                unclass <- setdiff(unclass, i)
                classn[reachables] <- classn[reachables] + 1
                while (length(reachables)) {
                  if (showplot) 
                    plot(data, col = 1 + cv, pch = 1 + isseed)
                  cv[reachables] <- cn
                  ap <- reachables
                  reachables <- integer()
                  if (method == "hybrid") {
                    tempdist <- distcomb(data[ap, , drop = FALSE], 
                      data[unclass, , drop = FALSE])
                    frozen.unclass <- unclass
                  }
                  for (i2 in seq(along = ap)) {
                    j <- ap[i2]
                    if (showplot > 1) 
                      plot(data, col = 1 + cv, pch = 1 + isseed)
                    if (method == "dist") {
                      jreachables <- unclass[data[j, unclass] <= 
                        eps]
                    }
                    else if (method == "hybrid") {
                      jreachables <- unclass[tempdist[i2, match(unclass, 
                        frozen.unclass)] <= eps]
                    }
                    else {
                      jreachables <- unclass[as.vector(distcomb(data[j, 
                        , drop = FALSE], data[unclass, , drop = FALSE])) <= 
                        eps]
                    }
                    if (length(jreachables) + classn[j] >= MinPts) {
                      isseed[j] <- TRUE
                      cv[jreachables[cv[jreachables] < 0]] <- cn
                      reachables <- union(reachables, jreachables[cv[jreachables] == 
                        0])
                    }
                    classn[jreachables] <- classn[jreachables] + 
                      1
                    unclass <- setdiff(unclass, j)
                  }
                }
            }
        }
        if (!length(unclass)) 
            break
    }
    rm(classn)
    if (any(cv == (-1))) {
        cv[cv == (-1)] <- 0
    }
    if (showplot) 
        plot(data, col = 1 + cv, pch = 1 + isseed)
    out <- list(cluster = cv, eps = eps, MinPts = MinPts)
    if (seeds && cn > 0) {
        out$isseed <- isseed
    }
    class(out) <- "dbscan"
    out
}

# main
Rdata <- read.table(inmatrix,row.names=1,header=T)
### transform to TPM, transcript per million reads(log10 scale)
Ndata <- apply(Rdata,2,getTPM)
covered_gene_number <- apply(Rdata,2,getcoverGnum)
highvargene <- selct_high_var_gene(Ndata,hvZ)

maxKnum <- min(maxKnum, ncol(Rdata) - 21)
tmp_data <- Rdata[highvargene,]
SIMLR_res = SIMLR(X=tmp_data,c=custom_k)
# Dimension reduction ressult : SIMLR_res$ydata

if(cortableY == 1){
    cortable <- cor(Ndata)
    write.table(cortable,file=paste(outname,'_correlation_table.txt',sep=""),row.names=T,col.names=T,sep="\t",quote=F)
}
cluster_plot <- paste(outname,'_Figure11_cluster.pdf',sep="")

if (clustering_method == 4){
    ### if user choose dbscan
    final_result <- givenE_dbscan(SIMLR_res$ydata,custom_d,cluster_plot)
}else{	
    KNUM <- custom_k
    final_result <- givenK_kmeans(SIMLR_res$ydata,RDnumber,KNUM,cluster_plot)
}

SpecificGeneScore <- function(in_group,out_group){
    in_group <- as.numeric(in_group)
    out_group <- as.numeric(out_group)
    in_mean <- mean(in_group)
    in_var <- var(in_group)
    out_mean <- mean(out_group)
    out_var <- var(out_group)
    if (out_mean == 0){
        if(in_mean == 0){
            return(0)
        }
        else{
            if (in_var == 0){
                in_var <- var(c(in_group[1]+1,in_group[2:length(in_group)]))
            }
            return((in_mean^2)/in_var)
        }
    }
    else{
        if (out_var == 0){
            out_var <- var(c(out_group[1]+1,out_group[2:length(out_group)]))
        }
        if (in_var == 0){
            in_var <- var(c(in_group[1]+1,in_group[2:length(in_group)]))
        }
        return((in_mean^2)/(in_var*out_mean^2*out_var))
    }
}
ClusterSpecificGene <- function(Rdata,highvargene,final_result,outname){
    for (each_cluster in unique(final_result[,3])){
        if (length(which(final_result[,3]!=each_cluster)) >= 3 & length(which(final_result[,3]==each_cluster)) >= 3){
            in_group <- Rdata[highvargene,which(final_result[,3]==each_cluster)]
            out_group <- Rdata[highvargene,which(final_result[,3]!=each_cluster)]
            all_genes <- highvargene
            all_score <- c()
            for (each_gene in all_genes){
                tmp_score <- SpecificGeneScore(in_group[each_gene,],out_group[each_gene,])
                all_score <- c(all_score,tmp_score)
            }
            write.table(cbind(all_genes[order(all_score,decreasing=T)][1:min(300,length(highvargene))]),file=paste(outname,"_specific_genes_of_cell_cluster",each_cluster,sep=""),quote=F,row.names=F,col.names=F)
        }
    }
}
ClusterSpecificGene(Rdata,highvargene,final_result,outname)

row.names(final_result) <- colnames(tmp_data)
if (clustering_method == 4){
    colnames(final_result) <- c("dimension_d1","dimension_d2","dbscan_cluster")
}else{
    colnames(final_result) <- c("dimension_d1","dimension_d2","kmeans_cluster")
}

write.table(final_result,file=paste(outname,'_cluster.txt',sep=""),row.names=T,col.names=T,sep="\t",quote=F)
