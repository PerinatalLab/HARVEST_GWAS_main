options(stringsAsFactors = F)
library(tidyr)
library(dplyr)
readVcf = function(chr, pos, pheno=NA, genome, project="H"){
	# determine the source vcf file
	if(project == "H"){
		vcf = paste0("/mnt/HUNT/merging/", genome, "_", chr, ".vcf.gz")
	} else if (project == "M"){
		vcf = paste0("~/data/geno/imputed/fresh/", chr, ".vcf.gz")
	} else {
		print("Please select project, either H (HARVEST) or M (MoBa08)")
		return()
	}
	
	# extract the snp
	header = scan(pipe(paste("bcftools query -l", vcf)), what="character")
	raw = read.table(pipe(paste0("bcftools view -r ", chr, ":", pos, " ", vcf)))
	
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
	
	geno$IID = header
	# attach phenotypes
	if(!is.na(pheno)){
		print("Attaching phenotype file")
		pheno = read.table(pheno, h=T)
		geno = inner_join(geno, pheno, by=c("IID"="id"))
	} else {
		print("Phenotype file not requested, returning genotypes only.")
	}
	
	return(geno)
}
