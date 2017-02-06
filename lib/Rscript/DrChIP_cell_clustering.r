a <-commandArgs(T)
# a <- c("H3K4me3_20130423_6_signal.txt","H3K4me3_20130423_6","0","20","20","/mnt/Storage/home/zhaocc/Work/0.DrSeq2_paper/DrChIP_Final/H3K4me3_20130423_6/mapping/H3K4me3_20130423_6.readsnumber.txt",2)
input_file <- a[1]
outname <- a[2]
cut_height <- as.numeric(a[3])
cell_cutoff <- as.numeric(a[4])
peak_cutoff <- as.numeric(a[5])
reads_number_file <- a[6]
given_cluster_number <- as.numeric(a[7])

ifEmptyCell <- function(signal_per_cell){
	return(sum(signal_per_cell) >= 3)
}

ref_tsne <-
function(X, initial_config = NULL, k=2, initial_dims=30, perplexity=30, max_iter = 1000, min_cost=0, epoch_callback=NULL,whiten=TRUE, epoch=100 ){
	#### original code from http://lvdmaaten.github.io/tsne/
	#### L.J.P. van der Maaten and G.E. Hinton. Visualizing High-Dimensional Data Using t-SNE. 2008
	if (class(X) == 'dist') { 
		n = attr(X,'Size')
		}
	else 	{
		X = as.matrix(X)
		X = X - min(X)
		X = X/max(X)
		initial_dims = min(initial_dims,ncol(X))
		if (whiten) X<-.whiten(as.matrix(X),n.comp=initial_dims)
		n = nrow(X)
	}

	momentum = .5
	final_momentum = .8
	mom_switch_iter = 250

	epsilon = 500
	min_gain = .01
	initial_P_gain = 4

	eps = 2^(-52) # typical machine precision

	if (!is.null(initial_config) && is.matrix(initial_config)) { 		
		if (nrow(initial_config) != n | ncol(initial_config) != k){
			stop('initial_config argument does not match necessary configuration for X')
		}
		ydata = initial_config
		initial_P_gain = 1
		
	} else {
		ydata = matrix(rnorm(k * n),n)
	}
	
	P = .x2p(X,perplexity, 1e-5)$P
	# P[is.nan(P)]<-eps
	P = .5 * (P + t(P))

	P[P < eps]<-eps
	P = P/sum(P)
	

	
	P = P * initial_P_gain
	grads =  matrix(0,nrow(ydata),ncol(ydata))
	incs =  matrix(0,nrow(ydata),ncol(ydata))
	gains = matrix(1,nrow(ydata),ncol(ydata))

	
	for (iter in 1:max_iter){
		if (iter %% epoch == 0) { # epoch
			cost =  sum(apply(P * log((P+eps)/(Q+eps)),1,sum))
			message("Epoch: Iteration #",iter," error is: ",cost)
			if (cost < min_cost) break
			if (!is.null(epoch_callback)) epoch_callback(ydata)

		}


		sum_ydata = apply(ydata^2, 1, sum)
		num =  1/(1 + sum_ydata +    sweep(-2 * ydata %*% t(ydata),2, -t(sum_ydata))) 
		diag(num)=0
		Q = num / sum(num)
		if (any(is.nan(num))) message ('NaN in grad. descent')
		Q[Q < eps] = eps
		stiffnesses = 4 * (P-Q) * num
		for (i in 1:n){
			grads[i,] = apply(sweep(-ydata, 2, -ydata[i,]) * stiffnesses[,i],2,sum)
		}
		
		gains = (gains + .2) * abs(sign(grads) != sign(incs)) + gains * .8 * abs(sign(grads) == sign(incs))		
		gains[gains < min_gain] = min_gain

		incs = momentum * incs - epsilon * (gains * grads)
		ydata = ydata + incs
		ydata = sweep(ydata,2,apply(ydata,2,mean))
		if (iter == mom_switch_iter) momentum = final_momentum
		
		if (iter == 100 && is.null(initial_config)) P = P/4
		

	
		
	}
	ydata
}

.Hbeta <-
function(D, beta){
	P = exp(-D * beta)
	sumP = sum(P)
	if (sumP == 0){
		H = 0
		P = D * 0
	} else {
		H = log(sumP) + beta * sum(D %*% P) /sumP
		P = P/sumP
	}
	r = {}
	r$H = H
	r$P = P
	r
}

.x2p <-
function(X,perplexity = 15,tol = 1e-5){
	if (class(X) == 'dist') {
		D = X
		n = attr(D,'Size')
	} else{
		D = dist(X)
		n = attr(D,'Size')
	}

	D = as.matrix(D)
	P = matrix(0, n, n )		
	beta = rep(1, n)
	logU = log(perplexity)
	
	for (i in 1:n){
		betamin = -Inf
		betamax = Inf
		Di = D[i, -i]
		hbeta = .Hbeta(Di, beta[i])
		H = hbeta$H; 
		thisP = hbeta$P
		Hdiff = H - logU;
		tries = 0;

		while(abs(Hdiff) > tol && tries < 50){
			if (Hdiff > 0){
				betamin = beta[i]
				if (is.infinite(betamax)) beta[i] = beta[i] * 2
				else beta[i] = (beta[i] + betamax)/2
			} else{
				betamax = beta[i]
				if (is.infinite(betamin))  beta[i] = beta[i]/ 2
				else beta[i] = ( beta[i] + betamin) / 2
			}
			
			hbeta = .Hbeta(Di, beta[i])
			H = hbeta$H
			thisP = hbeta$P
			Hdiff = H - logU
			tries = tries + 1
		}	
			P[i,-i]  = thisP	
	}	
	
	r = {}
	r$P = P
	r$beta = beta
	sigma = sqrt(1/beta)
	
	message('sigma summary: ', paste(names(summary(sigma)),':',summary(sigma),'|',collapse=''))

	r 
}

.whiten <-
function(X, row.norm=FALSE, verbose=FALSE, n.comp=ncol(X))
{  
	n.comp; # forces an eval/save of n.comp
	if (verbose) message("Centering")
   n = nrow(X)
	p = ncol(X)
	X <- scale(X, scale = FALSE)
   X <- if (row.norm) 
       t(scale(X, scale = row.norm))
   else t(X)

   if (verbose) message("Whitening")
   V <- X %*% t(X)/n
   s <- La.svd(V)
   D <- diag(c(1/sqrt(s$d)))
   K <- D %*% t(s$u)
   K <- matrix(K[1:n.comp, ], n.comp, p)
   X = t(K %*% X)
	X
}

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

reads_number_per_cell <- read.table(reads_number_file)
cell_barcode_id <- as.numeric(reads_number_per_cell[,2])
cell_labels <- c()
for (each in cell_barcode_id){
	cell_labels <- c(cell_labels,paste("C",each,sep=""))
}
reads_number <- as.numeric(reads_number_per_cell[,1])
names(reads_number) <- cell_labels
multiple_cells_cut <- quantile(reads_number,0.9)
no_cell_cut <- quantile(reads_number,0.1)
reads_cells <- cell_labels[which(reads_number<multiple_cells_cut & reads_number > no_cell_cut)]

pc_signal <- read.table(input_file,row.names=1,header=T) # peak*cell
cp_signal <- t(pc_signal) # cell*peak
selected_cells <- intersect(rownames(cp_signal),reads_cells)
if (length(selected_cells) == 0){
	print ("NONE of cells pass your cutoff,please loose your cutoff.")
}
cp_signal <- cp_signal[selected_cells,]

signal <- cp_signal[names(which(apply(cp_signal,1,RemoveUnwantedLine,peak_cutoff) == TRUE)),]
signal <- signal[,names(which(apply(signal,2,RemoveUnwantedLine,cell_cutoff) == TRUE))]
signal <- signal[apply(signal,1,sum)!=0,apply(signal,2,sum)!=0]
signal <- apply(signal,1,NormalizeSignal)
if (length(signal) == 0){
	print ("NONE of cells pass your cutoff,please loose your cutoff.")
}

signal <- t(signal)
cell_cor <- dist(cor(t(signal)))
# cell_cov <- dist(cov(t(signal)))
hc <- hclust(cell_cor)

pdf(paste(outname,"_Figure8_cell_clusting.pdf",sep=""),width=10,height=4)
# library("cluster")
# da <- diana(cov(t(signal)), diss = T)
# plot(da)
plot(hc, labels = FALSE, main = "hclust based on peak", xlab = "cells")
dev.off()

# k <- 10
# set.seed(1)
# km <- kmeans(t(signal),k)

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
    pdf(file=paste(outname,"_Figure9_silhouetteScore.pdf",sep=""))
    plot(silScore,main="")
    dev.off()
}
write.table(use_cluster_mat_sil,file = paste(outname,"_cluster_with_silhouette_score.txt",sep=""),quote=FALSE)

names(mycl) <- rownames(signal)
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
png(file = paste(outname,"_Figure10_heatmap.png",sep=""), width = 600, height = 600);
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
  # peak_cover <- apply(signal[names(mycl[mycl==i]),],2,sum)
  # other_cell_peak_cover <- apply(signal[names(mycl[mycl!=i]),],2,sum)
  # tmp_specific_peak <- c()
  # for (j in seq(length(peak_cover))){
  #   if (peak_cover[j] > 0 & other_cell_peak_cover[j] == 0) {
  #     tmp_specific_peak <- c(tmp_specific_peak,chartr("_", "\t",names(peak_cover)[j]))
  #   }
  # }
  # write.table(cbind(tmp_specific_peak),file=paste(outname,"_cluster",i,"_specific.peak.bed",sep=""),quote=F,col.names=F,row.names=F)
  write.table(cbind(names(mycl[mycl==i])),file=paste(outname,"_cluster",i,"_cells.txt",sep=""),quote=F,col.names=F,row.names=F)
}