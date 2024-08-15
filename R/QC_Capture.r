## Load data
library(RColorBrewer)
library(tidyverse)
library(viridis)
library(miscTools)

arguments = commandArgs(trailingOnly = TRUE)
args <- commandArgs(TRUE)
QCFILE <- as.character(arguments[1])

qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

QCdata <- read.table(paste0(QCFILE,".QC.txt"), header=T)
attach(QCdata)

for (i in SAMPLE) {
    Temp <- read.table(paste(".Temp_Frag/", i, ".txt", sep=""), header=T)
    assign(i, Temp)
}

nsamples <- length(SAMPLE)
samplecolours <- viridis(length(SAMPLE))

fragdata <- data.frame(matrix(ncol = 70, nrow = 0))

### frag calculations for text files
k=1
for (i in SAMPLE) {
    Temp <- get(i)
    Temp600 <- Temp %>% filter(Size < 600)
    Temp2 <- i
    Temp2 <- append(Temp2, nrow(Temp))
    Temp2 <- append(Temp2, mean(Temp600$Size))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 0 & Size < 90))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 0 & Size < 150))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 90 & Size < 150))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 140 & Size < 190))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 150 & Size < 600))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 0 & Size < 600))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 0 & Size < 10))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 10 & Size < 20))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 20 & Size < 30))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 30 & Size < 40))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 40 & Size < 50))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 50 & Size < 60))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 60 & Size < 70))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 70 & Size < 80))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 80 & Size < 90))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 90 & Size < 100))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 100 & Size < 110))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 110 & Size < 120))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 120 & Size < 130))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 130 & Size < 140))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 140 & Size < 150))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 150 & Size < 160))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 160 & Size < 170))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 170 & Size < 180))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 180 & Size < 190))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 190 & Size < 200))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 200 & Size < 210))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 210 & Size < 220))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 220 & Size < 230))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 230 & Size < 240))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 240 & Size < 250))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 250 & Size < 260))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 260 & Size < 270))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 270 & Size < 280))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 280 & Size < 290))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 290 & Size < 300))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 300 & Size < 310))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 310 & Size < 320))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 320 & Size < 330))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 330 & Size < 340))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 340 & Size < 350))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 350 & Size < 360))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 360 & Size < 370))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 370 & Size < 380))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 380 & Size < 390))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 390 & Size < 400))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 400 & Size < 410))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 410 & Size < 420))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 420 & Size < 430))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 430 & Size < 440))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 440 & Size < 450))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 450 & Size < 460))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 460 & Size < 470))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 470 & Size < 480))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 480 & Size < 490))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 490 & Size < 500))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 500 & Size < 510))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 510 & Size < 520))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 520 & Size < 530))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 530 & Size < 540))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 540 & Size < 550))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 550 & Size < 560))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 560 & Size < 570))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 570 & Size < 580))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 580 & Size < 590))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 590 & Size < 600))/nrow(Temp))
    Temp2 <- append(Temp2, nrow(Temp %>% filter(Size > 600))/nrow(Temp))
    if ( k == 1 ) {
        fragdata <- as.data.frame(t(Temp2))
    } else {
        fragdata <- rbind(fragdata, as.data.frame(t(Temp2)))
    }
    k=k+1
}
colnames(fragdata) <- c("SAMPLE", "TOTALPAIREDREADS", "AVGSIZE", "0-90", "0-150", "90-150", "140-190", "150-600", "0-600", "0-10", "10-20", "20-30", "30-40", "40-50", "50-60", "60-70", "70-80", "80-90", "90-100", "100-110", "110-120", "120-130", "130-140", "140-150", "150-160", "160-170", "170-180", "180-190", "190-200", "200-210", "210-220", "220-230", "230-240", "240-250", "250-260", "260-270", "270-280", "280-290", "290-300", "300-310", "310-320", "320-330", "330-340", "340-350", "350-360", "360-370", "370-380", "380-390", "390-400", "400-410", "410-420", "420-430", "430-440", "440-450", "450-460", "460-470", "470-480", "480-490", "490-500", "500-510", "510-520", "520-530", "530-540", "540-550", "550-560", "560-570", "570-580", "580-590", "590-600", "600+") 

write.table(fragdata, file=paste0(QCFILE, ".Frag.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

## Fragmentation plots per sample
limitchromosomes=1:24

for (i in SAMPLE) {
    Temp <- get(i)
    Temp <- Temp %>% filter(Size < 600)
    Temp$Chr <- gsub("chrX","23",Temp$Chr)
    Temp$Chr <- gsub("chrY","24",Temp$Chr)
    Temp$Chr <- gsub("chr","",Temp$Chr)
    Temp=Temp[Temp$Chr %in% limitchromosomes, ]
    Temp$Chr <- as.numeric(Temp$Chr)
    Temp <- Temp[order(Temp$Chr, Temp$BP),]
    NCHR=length(unique(Temp$Chr))

    BinFrag <- data.frame(matrix(ncol = 4, nrow = 0))

    for (k in unique(Temp$Chr)) {
        ChrTemp <- Temp %>% filter(Chr == k)
        for (j in seq(from=1, to=max(ChrTemp$BP), by=250000)) {
            l=j+999999
            BinFrag <- rbind(BinFrag, as.data.frame(t(c(k, j, mean(ChrTemp[ChrTemp$BP > j & ChrTemp$BP < l, ]$Size), nrow(ChrTemp[ChrTemp$BP > j & ChrTemp$BP < l, ])))))
        }
    }
    colnames(BinFrag) <- c("Chr", "BP", "Size", "Mass")
    BinFrag$MassNorm <- (BinFrag$Mass-min(BinFrag$Mass))/(max(BinFrag$Mass)-min(BinFrag$Mass))
    BinFrag <- BinFrag[complete.cases(BinFrag), ]
    write.table(BinFrag, file=paste0(i, ".Frag.txt"), sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

    if (NCHR==1) {
        Temp$BP2=Temp$BP
        png(file=paste0(i, ".Frag.png"), width=2000, height=1000)
	    par(mar=c(5,6,2,2))
		plot(Temp$BP2, Temp$Size, pch=21, cex=0.1, ylim=c(0,600), col="#5b5b5b08", bg="#5b5b5b08", xlab=paste("Chromosome",unique(Temp$Chr),"position"), ylab="Fragmentation", xaxt="n", cex.lab=2)
		axis(1, at = seq(head(Temp$BP2, n=1), tail(Temp$BP2, n=1), 500000)) 
        plot(BinFrag$BP2, BinFrag$Size, pch=21, cex=BinFrag$MassNorm*3, ylim=c(0,600), col="#5b5b5b08", bg="#5b5b5b08", xlab=paste("Chromosome",unique(BinFrag$Chr),"position"), ylab="Fragmentation", xaxt="n", cex.lab=2)
		axis(1, at = seq(head(BinFrag$BP2, n=1), tail(BinFrag$BP2, n=1), 500000)) 
        dev.off()
    } else {
        lastbase=0
	    Temp$BP2=NA
	    ticks=NULL
	    for (k in unique(Temp$Chr)) {
		    if (k==1) {
			    Temp[Temp$Chr==k, ]$BP2=Temp[Temp$Chr==k, ]$BP 
		    } else { 
			    lastbase=lastbase+tail(subset(Temp,Chr==k-1)$BP, 1)
			    Temp[Temp$Chr==k, ]$BP2=Temp[Temp$Chr==k, ]$BP+lastbase 
		    }
	        ticks=c(ticks, Temp[Temp$Chr==k, ]$BP2[floor(length(Temp[Temp$Chr==k, ]$BP2)/2)+1])
        }
        lastbase=0
	    BinFrag$BP2=NA
        for (k in unique(BinFrag$Chr)) {
            if (k==1) {
			    BinFrag[BinFrag$Chr==k, ]$BP2=BinFrag[BinFrag$Chr==k, ]$BP 
		    } else { 
			    lastbase=lastbase+tail(subset(BinFrag,Chr==k-1)$BP, 1)
			    BinFrag[BinFrag$Chr==k, ]$BP2=BinFrag[BinFrag$Chr==k, ]$BP+lastbase 
		    }
        }
        png(file=paste0(i, ".Frag.png"), width=2000, height=1000)
	    par(mar=c(5,6,2,2))
	    with(Temp, plot(BP2, Size, ylim=c(0,600), ylab="Fragmentation", xlab="Chromosome", cex.lab=2, cex.axis=2, xaxt="n", type="n"))
	    axis(1, at=ticks, lab=unique(Temp$Chr))
		for (k in seq(head(unique(Temp$Chr),1),tail(unique(Temp$Chr),1),2)) {
			with(Temp[Temp$Chr==k, ], rect(head(BP2,1),0,tail(BP2,1),600, border="grey90", col="grey90"))
		}
	    icol=1
		for (k in unique(Temp$Chr)) {
		    with(Temp[Temp$Chr==k, ],points(BP2, Size, pch=21, col="#5b5b5b08", bg="#5b5b5b08", cex=0.1))
		    icol=icol+1
        }
  	    icol=1
		for (k in unique(BinFrag$Chr)) {
		    with(BinFrag[BinFrag$Chr==k, ],points(BP2, Size, pch=21, col="black", bg="red", cex=BinFrag$MassNorm*3))
		    icol=icol+1
        }
        garbage <- dev.off()
    }
}
## Global QC plot
png(file=paste0(QCFILE, ".QC.png"), width=500+(nsamples*200), height=1000)
par(mfrow = c(3,3), oma = c(2,2,1,1) + 0.1, mar = c(1,1,1,1) + 0.1)
barplot(TOTAL_READS, main="Total Reads", cex.main=1.4, ylim=c(0,max(c(10000000, max(TOTAL_READS)+1000000))), col=samplecolours)
abline(5000000,0)
barplot(MEAN_TARGET_COVERAGE, main="Coverage", cex.main=1.4, ylim=c(0,max(c(1000, max(MEAN_TARGET_COVERAGE)))), col=samplecolours)
abline(500,0)
barplot(FOLD_80_BASE_PENALTY, ylim=c(0,max(c(8, max(FOLD_80_BASE_PENALTY)))), main="Fold80", cex.main=1.4, col=(samplecolours))
abline(3.5,0)
barplot(rbind(PCT_TARGET_BASES_100X,PCT_PF_UQ_READS_ALIGNED,STRAND_BALANCE), beside=T, ylim=c(0,1), main="Incl Percts", cex.main=1.4, col=brewer.pal(3, "Greys"))
abline(0.5,0)
abline(0.9,0)
legend("bottomright", c("% bases 100x", "% Aligned", "Strand"), inset=.02, horiz=F, cex=1.4, fill=brewer.pal(3, "Greys"), bg="#ffffffb3")
barplot(rbind(PERCENT_DUPLICATION,PCT_OFF_BAIT,PCT_EXC_OVERLAP,PCT_EXC_MAPQ,PCT_EXC_BASEQ), beside=T, ylim=c(0,1), main="Excl Percts", cex.main=1.4, col=brewer.pal(5, "Greys"))
abline(0.2,0)
abline(0.05,0)
legend("topright", c("% Dups", "% OffBait", "% Overlap", "% LowMap", "% LowQ"), inset=.02, horiz=F, cex=1.4, fill=brewer.pal(5, "Greys"), bg="#ffffffb3")
barplot(rbind(PCT_PF_READS_IMPROPER_PAIRS,contamination,PCT_ADAPTER,PF_MISMATCH_RATE,PF_INDEL_RATE,PCT_CHIMERAS), beside=T, ylim=c(0,max(0.03,PCT_PF_READS_IMPROPER_PAIRS,contamination,PCT_ADAPTER,PF_MISMATCH_RATE,PF_INDEL_RATE,PCT_CHIMERAS)), main="Misc", cex.main=1.4, col=brewer.pal(6, "Greys"))
abline(0.01,0)
legend("topright", c("% Unpaired", "Contamination", "% Adapter", "% Mismatches", "% Indels", "% Chimeras"), inset=.02, horiz=F, cex=1.4, fill=brewer.pal(6, "Greys"), bg="#ffffffb3")
barplot(MEAN_INSERT_SIZE, main="Insert Size", cex.main=1.4, ylim=c(0,600), col=samplecolours)
abline(300,0)
k=1
upperlimit=0
for (i in SAMPLE) {
    Temp <- get(i)
    Temp2 <- Temp %>% filter(Size < 600)
    Temp3 <- density(Temp2$Size)
    upperlimittemp <- max(Temp3$y)
    if ( upperlimittemp > upperlimit) {
        upperlimit <- upperlimittemp
    }
}
for (i in SAMPLE) {
    Temp <- get(i)
    Temp2 <- Temp %>% filter(Size < 600)
    Temp3 <- density(Temp2$Size)
    if ( k == 1 ) {
        plot(Temp3, main="Fragmentation", cex.main=1.4, ylim=c(0,max(upperlimit)+0.0004), col=samplecolours[k], lwd=2)
    } else {
        lines(Temp3, col=samplecolours[k], lwd=2) 
    }
    k=k+1
}
plot.new()
legend("center", as.vector(SAMPLE), inset=.02, horiz=F, cex=1.4, fill=samplecolours, bg="#ffffffb3", ncol=2)
garbage <- dev.off()
