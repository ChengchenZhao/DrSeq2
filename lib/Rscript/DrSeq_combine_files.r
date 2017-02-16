a <- commandArgs(T)
files <- strsplit(a[1],",")[[1]]
cell <- a[2]
outname <- a[3]

all_gene <- c()
all_cell <- c()
all_exp <- c()
exp1 <- read.table(files[1],header=T,row.names=1)

if (cell == "row"){
	exp1 <- t(exp1)
}

gene1 <- rownames(exp1)
cell1 <- colnames(exp1)
	
all_gene <- gene1
all_cell <- cell1
all_exp <- exp1

if (length(files) >= 2) {
	for (each_file in files[2:length(files)]){
		tmp_exp <- read.table(each_file,header=T,row.names=1)
		if (cell == "row"){
			tmp_exp <- t(tmp_exp)
		}
		tmp_gene <- rownames(tmp_exp)
		tmp_cell <- colnames(tmp_exp)
		comm_gene <- intersect(tmp_gene,all_gene)
		all_exp <- cbind(all_exp[comm_gene,],tmp_exp[comm_gene,])
		all_gene <- comm_gene
		all_cell <- c(all_cell,tmp_cell)
	}
}

colnames(all_exp) <- all_cell
rownames(all_exp) <- all_gene

write.table(all_exp,file=paste(outname,"_expmat.txt",sep=''),quote=F)