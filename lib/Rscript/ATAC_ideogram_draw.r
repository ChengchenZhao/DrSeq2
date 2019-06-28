a<-commandArgs(T)
cytoBandfile <- a[1]
outname <- a[2]
k <- as.numeric(a[3])
sample_files <- strsplit(a[4],",")[[1]]
org <- a[5]

if (org == "hs"){
    organism <- "hsa"
}
if (org == "mm"){
    organism <- "mmu"
}
transCoordinate <- function (chrompos = NULL,  organism, sexChromosomes = FALSE, units = "hg19") 
{
    if (!missing(organism)) {
        organism <- match.arg(organism, c("hsa", "mmu", "rno"))
        chrom.n <- switch(organism, hsa = 22, mmu = 19, rno = 20)
        if (is.null(chrompos)) {
            chroms = c(1:chrom.n, if (sexChromosomes) c(98, 99) else NULL)
            chrnames = characterCHR(chroms, "chr")
            mapinfo = rep(0, length(chroms))
            chrompos = cbind(CHR = chroms, MapInfo = mapinfo)
            rownames(chrompos) <- chrnames
        }
        chrs2 <- factor(numericCHR(chrompos[, "CHR"]), levels = c(1:chrom.n, 
            if (sexChromosomes) c(98, 99) else NULL))
        if (organism %in% c("hsa", "mmu")) 
            lens <- lengthChromosome(levels(chrs2), units = units)
        else lens <- sapply(split(chrompos[, "MapInfo"], chrs2), 
            function(x) max(c(0, x)))
        names(lens) <- characterCHR(names(lens))

        dwidth <- NULL
        for (i in 1:(chrom.n%/%2)) dwidth[i] <- lens[i] + lens[chrom.n + 
            1 - i]
        if (chrom.n%%2 == 1) 
            dwidth <- c(dwidth, lens[chrom.n%/%2 + 1])
        if (sexChromosomes) 
            dwidth <- c(dwidth, lens["X"] + lens["Y"])
        maxdwidth <- max(dwidth) * 1.05
        leftrow <- c(if (sexChromosomes) "X" else NULL, ((chrom.n + 
            1)%/%2):1)
        rightrow <- c(if (sexChromosomes) "Y" else NULL, if (chrom.n%%2 == 
            1) "" else NULL, ((chrom.n + 1)%/%2 + 1):chrom.n)

        dchrompos <- matrix(0, nrow = nrow(chrompos), ncol = 2, 
            dimnames = list(rownames(chrompos), c("CHR", "MapInfo")))
        for (i in 1:length(rightrow)) if (rightrow[i] != "") {
            probes <- characterCHR(chrompos[, "CHR"]) == rightrow[i]
            dchrompos[probes, 2] <- chrompos[probes, "MapInfo"] + 
                maxdwidth - lens[rightrow[i]]
            dchrompos[probes, 1] <- i
        }
        for (i in 1:length(leftrow)) {
            probes <- characterCHR(chrompos[, "CHR"]) == leftrow[i]
            dchrompos[probes, 2] <- chrompos[probes, "MapInfo"]
            dchrompos[probes, 1] <- i
        }
    }
    dchrompos
}
prepareGenomePlot <- function (chrompos = NULL, cols = "grey50", paintCytobands = FALSE, 
    bleach = 0, topspace = 1, organism, sexChromosomes = FALSE, 
    units = "hg19", ...) 
{
    cytobandWidth <- 0.075
    par(mar = c(1, 4, 2, 3) + 0.1)
    if (!missing(organism)) {
        organism <- match.arg(organism, c("hsa", "mmu", "rno"))
        chrom.n <- switch(organism, hsa = 22, mmu = 19, rno = 20)
        if (is.null(chrompos)) {
            chroms = c(1:chrom.n, if (sexChromosomes) c(98, 99) else NULL)
            chrnames = characterCHR(chroms, "chr")
            mapinfo = rep(0, length(chroms))
            chrompos = cbind(CHR = chroms, MapInfo = mapinfo)
            rownames(chrompos) <- chrnames
        }
        chrs2 <- factor(numericCHR(chrompos[, "CHR"]), levels = c(1:chrom.n, 
            if (sexChromosomes) c(98, 99) else NULL))
        if (organism %in% c("hsa", "mmu")) 
            lens <- lengthChromosome(levels(chrs2), units = units)
        else lens <- sapply(split(chrompos[, "MapInfo"], chrs2), 
            function(x) max(c(0, x)))
        names(lens) <- characterCHR(names(lens))
        cols <- rep(cols, length.out = length(lens))
        names(cols) <- names(lens)
        dwidth <- NULL
        for (i in 1:(chrom.n%/%2)) dwidth[i] <- lens[i] + lens[chrom.n + 
            1 - i]
        if (chrom.n%%2 == 1) 
            dwidth <- c(dwidth, lens[chrom.n%/%2 + 1])
        if (sexChromosomes) 
            dwidth <- c(dwidth, lens["X"] + lens["Y"])
        maxdwidth <- max(dwidth) * 1.05
        leftrow <- c(if (sexChromosomes) "X" else NULL, ((chrom.n + 
            1)%/%2):1)
        rightrow <- c(if (sexChromosomes) "Y" else NULL, if (chrom.n%%2 == 
            1) "" else NULL, ((chrom.n + 1)%/%2 + 1):chrom.n)
        plot(c(0, maxdwidth), c(0.5, 0.5 + length(dwidth) + topspace), 
            type = "n", ylab = "Chromosome", xlab = "", axes = FALSE, 
            las = 2, ...)
        axis(2, c(1:length(dwidth)), characterCHR(leftrow), las = 2)
        axis(4, c(1:length(dwidth)), characterCHR(rightrow), 
            las = 2)
        if (paintCytobands && organism %in% c("hsa", "mmu")) {
            for (i in 1:length(dwidth)) {
                if (lens[leftrow[i]] > 0) 
                  paintCytobands(leftrow[i], c(0, i + cytobandWidth/2), 
                    units = units, width = cytobandWidth, length.out = lens[leftrow[i]], 
                    legend = FALSE, bleach = bleach)
                if (rightrow[i] != "" && lens[rightrow[i]] > 
                  0) 
                  paintCytobands(rightrow[i], c(maxdwidth - lens[rightrow[i]], 
                    i + cytobandWidth/2), units = units, width = cytobandWidth, 
                    length.out = lens[rightrow[i]], legend = FALSE, 
                    bleach = bleach)
            }
        }
        else {
            for (i in 1:length(dwidth)) {
                lines(c(0, lens[leftrow[i]]), c(i, i), col = cols[leftrow[i]], 
                  lwd = 2)
                if (rightrow[i] != "") 
                  lines(c(maxdwidth - lens[rightrow[i]], maxdwidth), 
                    c(i, i), col = cols[rightrow[i]], lwd = 2)
            }
        }
        dchrompos <- matrix(0, nrow = nrow(chrompos), ncol = 2, 
            dimnames = list(rownames(chrompos), c("CHR", "MapInfo")))
        for (i in 1:length(rightrow)) if (rightrow[i] != "") {
            probes <- characterCHR(chrompos[, "CHR"]) == rightrow[i]
            dchrompos[probes, 2] <- chrompos[probes, "MapInfo"] + 
                maxdwidth - lens[rightrow[i]]
            dchrompos[probes, 1] <- i
        }
        for (i in 1:length(leftrow)) {
            probes <- characterCHR(chrompos[, "CHR"]) == leftrow[i]
            dchrompos[probes, 2] <- chrompos[probes, "MapInfo"]
            dchrompos[probes, 1] <- i
        }
    }
    else {
        chrs2 <- factor(numericCHR(chrompos[, "CHR"]))
        lens <- sapply(split(chrompos[, "MapInfo"], chrs2), max)
        m <- length(lens)
        cols <- rep(cols, length.out = m)
        maxdwidth <- max(lens)
        plot(c(0, maxdwidth), c(0.5, m + 0.5 + topspace), type = "n", 
            ylab = "Chromosome", xlab = "", axes = FALSE, las = 2, 
            ...)
        axis(2, c(m:1), characterCHR(names(lens)), las = 2)
        for (i in 1:m) lines(c(0, lens[i]), c(m + 1 - i, m + 
            1 - i), col = cols[as.numeric(names(lens))], lwd = 2)
        dchrompos <- chrompos
        dchrompos[, 1] <- m + 1 - as.numeric(chrs2)
    }
    dchrompos
}
characterCHR <- function (CHR, prefix = "") 
{
    CHR <- as.character(CHR)
    CHR[CHR == "98"] <- "X"
    CHR[CHR == "99"] <- "Y"
    CHR[CHR == "100"] <- "XY"
    CHR[CHR == "101"] <- "MT"
    CHR[CHR == "102"] <- "-"
    paste(prefix, CHR, sep = "")
}
numericCHR <- function (CHR, prefix = "chr") 
{
    CHR <- sub(paste0("^", prefix), "", as.character(CHR))
    CHR[grep("_", CHR)] <- "102"
    CHR[CHR == "X"] <- "98"
    CHR[CHR == "Y"] <- "99"
    CHR[CHR == "XY"] <- "100"
    CHR[CHR %in% c("M", "MT")] <- "101"
    as.numeric(CHR)
}
lengthChromosome <- function (chrom, units = "hg19") 
{
    if (class(units) == "data.frame") {
        chrombands <- .convertUCSCcytoband(units)
    }
    else {
        chrombands <- .getChrombands(units)
    }
    chrom <- characterCHR(chrom)
    chromdata <- subset(chrombands, chrombands$chr %in% chrom)
    if (nrow(chromdata) == 0) {
        warning("chromosomes not found in chromosomal data")
        res = rep(NA, length(chrom))
    }
    else {
        chromlengths <- aggregate(chromdata$segend, list(chromdata$chr), 
            function(x) x[length(x)])
        res <- chromlengths[match(chrom, chromlengths[, 1]), 
            2]
    }
    names(res) <- chrom
    return(res)
}
.convertUCSCcytoband <- function (ucscdata) 
{
    data.frame(chr = characterCHR(numericCHR(ucscdata[, 1])), 
        segstart = ucscdata[, 2], segend = ucscdata[, 3], stain = ucscdata[, 
            5], band = substring(ucscdata[, 4], 2), arm = substring(ucscdata[, 
            4], 1, 1), stringsAsFactors = FALSE)
}
paintCytobands <- function (chrom, pos = c(0, 0), units = "hg19", width = 0.4, 
    length.out, bands = "major", orientation = c("h", "v"), legend = TRUE, 
    cex.leg = 0.7, bleach = 0, ...) 
{
    bleacher <- function(x) {
        (x * (1 - bleach)) + bleach
    }
    if (class(units) == "data.frame") {
        chrombands <- .convertUCSCcytoband(units)
    }
    else {
        chrombands <- .getChrombands(units)
    }
    chrom <- switch(as.character(chrom), `98` = "X", `99` = "Y", 
        as.character(chrom))
    orientation <- match.arg(orientation)
    if (length(pos) == 1) 
        pos <- c(0, pos)
    chromdata <- subset(chrombands, chrombands$chr == chrom)
    if (nrow(chromdata) > 0) {
        lc <- nchar(chromdata$band)
        sel <- !(substr(chromdata$band, lc, lc) %in% letters)
        if (bands != "major") 
            sel <- !sel
        chromdata <- chromdata[sel, ]
        rm(lc, sel)
        type.b <- match(chromdata$stain, c("acen", "gneg", "gpos", 
            "gvar", "stalk", "gpos25", "gpos50", "gpos75", "gpos100", 
            "gpos33", "gpos66"))
        bandpos <- chromdata[, c("segstart", "segend")]
        bandcol <- gray(bleacher(c(0.5, 1, 0.2, 0.6, 0.75, 0.7, 
            0.5, 0.3, 0.1, 0.6, 0.4)))[type.b]
        banddens <- c(30, -1, -1, -1, 10, -1, -1, -1, -1, -1, 
            -1)[type.b]
        bandbord <- gray(bleacher(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 
            0, 0)))[type.b]
        if (!missing(length.out)) {
            bandpos <- (bandpos/max(bandpos)) * length.out
        }
        n <- nrow(chromdata)
        centromere <- which(chromdata$arm[-n] != chromdata$arm[-1])
        if (length(centromere == 1)) {
            idx <- c(2:(centromere - 1), (centromere + 2):(n - 
                1))
        }
        else {
            idx <- c(2:(n - 1))
        }
        if (orientation == "h") {
            rect(pos[1] + bandpos[idx, 1], pos[2], pos[1] + bandpos[idx, 
                2], pos[2] - width, col = bandcol[idx], density = banddens[idx], 
                border = bandbord[idx])
            qs.semicircle(pos[1] + bandpos[1, 2], pos[2] - width, 
                width, bandpos[1, 2] - bandpos[1, 1], 2, col = bandcol[1], 
                density = banddens[1], border = bandbord[1], 
                ...)
            qs.semicircle(pos[1] + bandpos[n, 1], pos[2] - width, 
                width, bandpos[n, 2] - bandpos[n, 1], 4, col = bandcol[n], 
                density = banddens[n], border = bandbord[n], 
                ...)
            if (length(centromere == 1)) {
                qs.semicircle(pos[1] + bandpos[centromere, 1], 
                  pos[2] - width, width, bandpos[centromere, 
                    2] - bandpos[centromere, 1], 4, col = bandcol[centromere], 
                  density = banddens[centromere], border = bandbord[centromere], 
                  ...)
                qs.semicircle(pos[1] + bandpos[centromere + 1, 
                  2], pos[2] - width, width, bandpos[centromere + 
                  1, 2] - bandpos[centromere + 1, 1], 2, col = bandcol[centromere + 
                  1], density = banddens[centromere + 1], border = bandbord[centromere + 
                  1], ...)
                centromere.size = 0.6 * 0.5 * width/yinch(1)
                symbols(pos[1] + bandpos[centromere, 2], pos[2] - 
                  0.5 * width, circles = 1, inches = centromere.size, 
                  add = TRUE, fg = gray(bleacher(0)), bg = "white", 
                  ...)
            }
            if (legend) 
                text(pos[1] + (bandpos[, 1] + bandpos[, 2])/2, 
                  pos[2] + 0.5 * width, paste(chromdata[, "arm"], 
                    chromdata[, "band"], sep = ""), adj = c(0, 
                    0.5), srt = 90, cex = cex.leg, ...)
        }
        else {
            rect(pos[1], pos[2] - bandpos[idx, 1], pos[1] - width, 
                pos[2] - bandpos[idx, 2], col = bandcol[idx], 
                density = banddens[idx], border = bandbord[idx], 
                ...)
            qs.semicircle(pos[1] - width, pos[2] - bandpos[1, 
                2], width, bandpos[1, 2] - bandpos[1, 1], 3, 
                col = bandcol[1], density = banddens[1], border = bandbord[1], 
                ...)
            qs.semicircle(pos[1] - width, pos[2] - bandpos[n, 
                1], width, bandpos[n, 2] - bandpos[n, 1], 1, 
                col = bandcol[n], density = banddens[n], border = bandbord[n], 
                ...)
            if (length(centromere == 1)) {
                qs.semicircle(pos[1] - width, pos[2] - bandpos[centromere, 
                  1], width, bandpos[centromere, 2] - bandpos[centromere, 
                  1], 1, col = bandcol[centromere], density = banddens[centromere], 
                  border = bandbord[centromere], ...)
                qs.semicircle(pos[1] - width, pos[2] - bandpos[centromere + 
                  1, 2], width, bandpos[centromere + 1, 2] - 
                  bandpos[centromere + 1, 1], 3, col = bandcol[centromere + 
                  1], density = banddens[centromere + 1], border = bandbord[centromere + 
                  1], ...)
                centromere.size = 0.6 * 0.5 * width/xinch(1)
                symbols(pos[1] - 0.5 * width, pos[2] - bandpos[centromere, 
                  2], circles = 1, inches = centromere.size, 
                  add = TRUE, fg = gray(bleacher(0)), bg = "white", 
                  ...)
            }
            if (legend) 
                text(pos[1] + 0.5 * width, pos[2] - (bandpos[, 
                  1] + bandpos[, 2])/2, paste(chromdata[, "arm"], 
                  chromdata[, "band"], sep = ""), adj = c(0, 
                  0.5), srt = 0, cex = cex.leg, ...)
        }
    }
    else {
        warning(paste("Chromosome", chrom, "is not plotted because cytoband data is not available"))
    }
}
qs.semicircle <- function (base.x, base.y, base.length, height = base.length, 
    side = 1, orientation = NULL, plottype = "poly", ...) 
{
    radius <- base.length/2
    x <- radius * seq(-1, 1, length = 40)
    y <- height/radius * sqrt(radius^2 - x^2)
    if (is.null(orientation)) {
        co <- as.integer(cos(pi * (3 - side)/2))
        so <- as.integer(sin(pi * (3 - side)/2))
    }
    else {
        co <- cos(orientation)
        so <- sin(orientation)
    }
    tx <- co * x - so * y
    ty <- so * x + co * y
    if (is.null(orientation)) {
        if (side == 1 || side == 3) {
            base.x <- base.x + radius
        }
        else if (side == 2 || side == 4) {
            base.y <- base.y + radius
        }
    }
    x <- base.x + tx
    y <- base.y + ty
    switch(plottype, poly = polygon(x, y, ...), line = lines(x, 
        y, ...))
}
# temp <- tempfile(fileext = ".txt.gz")
# download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz",temp)
cytobands <- read.table(cytoBandfile,sep="\t")

png(paste(outname,"_Figure9_ideogram.png",sep=""),width=8,height=8)
prepareGenomePlot(paintCytobands=TRUE,organism=organism,sexChromosomes=F,units=cytobands)
sample <- seq(k)
cols <- rainbow(k,alpha = 0.5)
names(cols) <- sample
for(each in sample){
    d <- read.table(sample_files[each], header=F, sep='\t')
    d <- d[which(d[,1]!='chrX' & d[,1]!='chrY'),]
    CHR <- as.numeric(gsub('chr', '', d[,1]))
    MapInfo <- (d[,2] + d[,3])/2
    tmp <- transCoordinate(data.frame(CHR,MapInfo),organism=organism,sexChromosomes=F,units=cytobands)
    points(tmp[,2],tmp[,1]+0.5,pch="|",col=cols[each])
}
# legend('bottomright',legend=sample,col=cols,lwd=2,bty='n',cex=0.8)
title('specific region on genomic ideogram')
dev.off()
