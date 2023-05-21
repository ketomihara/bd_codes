#load libraries
library(dplyr)
library(ggplot2)
library(stringr)
library(reshape2)

###read and reformat vcf files###
list <- c("bd_bd", "bd_plus", "p50T")
chr_num <- paste(rep("Bomo_Chr", 28), c(1:28), sep="")

for(i in list){
	x <- read.table(paste("E:/workspace/bd_RNASeq/", i, "_processed/", i, "_filtered.vcf", sep=""), comment.char="#", 
		sep="\t", header=F, quote="", stringsAsFactors = F)
	colnames(x) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "DESCRIP")
	y <- x[grep("0/1|1/2", x$DESCRIP),]
	z <- table(y$CHROM[y$CHROM %in% chr_num])
	variable_name1 <- str_sub(i)
	variable_name2 <- str_sub(paste(i, "_hetero", sep=""))
	variable_name3 <- str_sub(paste(i, "_hetero_count", sep=""))
	assign(variable_name1, x)
	assign(variable_name2, y)
	assign(variable_name3, z)
}

###count heterozygous SNPs and visualize them###
hetero_count <- rbind(bd_bd_hetero_count, bd_plus_hetero_count, p50T_hetero_count)
colnames(hetero_count) <- str_sub(colnames(hetero_count), start=9, end=-1)
head(hetero_count)

chr_length <- read.table("E:/workspace/reference_genome/Bmori_genome_length.txt", 
		sep="\t", header=T, quote="", stringsAsFactors = F)
chr_length2 <- chr_length[c(1:28),]
chr_length2$seqid <- str_sub(chr_length2$seqid, start=9, end=-1)

hetero_count_melt <- melt(hetero_count)
hetero_count_melt2 <- transform(hetero_count_melt, Var2= factor(Var2, 
	levels = c(1:28)))
hetero_count_chrlength <- merge(hetero_count_melt2, chr_length2, by.x="Var2", by.y="seqid", all.x=T)
hetero_count_chrlength$density <- (hetero_count_chrlength$value/hetero_count_chrlength$length)*1000000

g_whole_chr_het_density <- ggplot(hetero_count_chrlength, aes(x = Var2, y= density, fill = Var1))+
	geom_bar(stat="identity", position="dodge")+
	labs(x="Chromosome", y="Density of heterozygous SNPs (/Mb)")+
	scale_fill_manual(name = "Genotype", labels=c(bd_bd_hetero_count="bd/bd",bd_plus_hetero_count="bd/+",
		p50T_hetero_count="WT (Dazao)"), values=c(bd_bd_hetero_count="#56B4E9",bd_plus_hetero_count="#E69F00",
		p50T_hetero_count="#000000"))+
	scale_y_continuous(expand = c(0,0))+
	theme_classic()

pdf("whole_chr_het_density_revise.pdf", height=3.2, width=6)
g_whole_chr_het_density
dev.off()

svg("whole_chr_het_density_revise.svg", height=3.2, width=6)
g_whole_chr_het_density
dev.off()


###visualize hetero SNPs across chr9###
list2 <- list(bd_bd_hetero, bd_plus_hetero, p50T_hetero)
data_chr9 <- data.frame(matrix(rep(NA, 2), nrow=1))[numeric(0), ]
	colnames(data_chr9) <- c("sample", "position")
for(j in 1:3){
	x <- list2[[j]]
	y <- x[x$CHROM %in% "Bomo_Chr9", ]
	data_chr9 <- rbind(data_chr9, cbind(sample=rep(as.character(list[j]), nrow(y)),
		position=as.numeric(y$POS)), stringsAsFactors = FALSE)
	data_chr9$position <- as.numeric(data_chr9$position)
}

g9 <- ggplot(data_chr9, aes(x = position, fill = sample))+ 
	geom_histogram(position = "identity", alpha=0.5, binwidth=1000000, boundary = 1)+
	labs(x="Position in chromosome 9 (Mb)", y="The number of heterozygous SNPs (/Mb)")+
	scale_x_continuous(expand = c(0,0), breaks=c(0,5000000,10000000,15000000), labels=c( "0","5", "10","15"))+
	scale_fill_manual(name = "Genotype", labels=c(bd_bd="bd/bd",bd_plus="bd/+",p50T="WT (Dazao)"),
	values=c(bd_bd="#56B4E9",bd_plus="#E69F00",p50T="#000000"))+
	scale_y_continuous(expand = c(0,0))+
	theme_classic()

g9_info <- ggplot_build(g9)
print(g9_info$data)

#g9 was not included into the final version of the article.
pdf("Chr9_SNPs.pdf", height=3.2, width=6)
g9
dev.off()

#bd_plus only
g9_2 <- ggplot(data_chr9[data_chr9$sample %in% "bd_plus",], aes(x = position))+ 
	geom_histogram(binwidth=1000000, boundary = 1, fill="#E69F00")+
	labs(x="Position in chromosome 9 (Mb)", y="The number of heterozygous SNPs (/Mb)")+
	scale_x_continuous(expand = c(0,0), breaks=c(0,5000000,10000000,15000000), labels=c( "0","5", "10","15"))+
	scale_y_continuous(expand = c(0,0))+
	theme_classic()

pdf("Chr9_SNPs_revise.pdf", height=3.2, width=4.7)
g9_2
dev.off()
