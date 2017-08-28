options(stringsAsFactors = F)
library(tidyr)
library(dplyr)
readVcf = function(chr, pos, pheno, genome){
	# determine the source vcf file
	vcf = paste0("/mnt/HUNT/merging/", genome, "_", chr, ".vcf.gz")
	
	# extract the snp
	header = scan(pipe(paste("bcftools query -l", vcf)), what="character")
	raw = read.table(pipe(paste0("bcftools view -r ", chr, ":", pos, " ", vcf)))
	pheno = read.table(pheno, h=T)
	
	# transform the geno file
	geno = raw[,-(1:9)]
	info = raw[,(1:9)]
	geno = as.data.frame(t(geno))

	geno = separate(geno, "V1", into = strsplit(info$V9, ":")[[1]], sep=":")
	geno = separate(geno, "GP", into = c("AA", "AB", "BB"), sep=",")
	geno$DS = as.numeric(geno$DS)
	geno$AA = as.numeric(geno$AA)
	geno$AB = as.numeric(geno$AB)
	geno$BB = as.numeric(geno$BB)
	
	# attach phenotypes
	geno$IID = header
	geno = inner_join(geno, pheno, by=c("IID"="SentrixID_1"))
	
	return(geno)
}
