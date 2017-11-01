#!/usr/bin/Rscript

############### COMMON INIT #####################

options(stringsAsFactors = F)
library(dplyr)
library(tidyr)
library(ggplot2)

harvdir = "~/data/geno/harvest-aux/"
flags = read.table(paste0(harvdir, "harvest-flag-list.txt"), h=T)
link = read.table("~/data/mobaqs/p1724/harvest_linkage.csv", sep=";", h=T)
linkm = coremoms = corekids = corepaired = data.frame()

# Usage: attachPheno(mfr) or attachPheno(q1)
attachPheno = function(phenodf){
		linkm <<- inner_join(link, phenodf, by="PREG_ID_1724")
	print(sprintf("%i individuals with phenotype info found", nrow(linkm)))
}

# Usage: getCore(linkm)
getCore = function(pheno){
	coremoms <<- filter(flags, genotypesOK, phenotypesOK, coreOK) %>%
		semi_join(pheno, ., by=c("SentrixID_1"="IID")) %>%
		filter(Role=="Mother")
		# 7065, incl. multiple pregs
		print(sprintf("%i core mothers found", nrow(coremoms)))

	corekids <<- filter(flags, genotypesOK, phenotypesOK, coreOK) %>%
		semi_join(pheno, ., by=c("SentrixID_1"="IID")) %>%
		filter(Role=="Child")
	print(sprintf("%i core kids found", nrow(corekids)))
}

# removes repeated pregnancies,
# keeping the one included in the other core set.
# Usage: removeRepeated(coremoms, corekids) to clean mothers
# 		removeRepeated(corekids, coremoms) to clean kids
#			- obviously not needed, because kids don't have multiple MFR rows
removeRepeated = function(dftoclean, basedon){
	corepaired = group_by(dftoclean, SentrixID_1) %>%
		mutate(hasPair = PREG_ID_1724 %in% basedon$PREG_ID_1724) %>%
		top_n(1, hasPair) %>%
		filter(rank(PREG_ID_1724)==1)
	print(sprintf("after removing repeated pregnancies, %i remain", nrow(corepaired)))
	return(corepaired)
}

# Usage: attachCovariates AA87 moms 6
attachCovariates = function(corepaired, dataset, numPCs){
	## add batch covariate
	corepaired = select(flags, one_of(c("IID", "BATCH"))) %>%
		left_join(corepaired, ., by=c("SentrixID_1" = "IID")) %>%
		mutate(BATCH = as.numeric(BATCH=="M24"))

	# read in correct ibd-exclusion list and pca-covar file
	if(dataset=="moms"){
		ibd = read.table(paste0(harvdir, "ibd_pihat_exclude_mothers"))
		pcs = read.table(paste0(harvdir, "plink_covar_mothers"))
	} else if (dataset=="fets") {
		ibd = read.table(paste0(harvdir, "ibd_pihat_exclude_offspring"))
		pcs = read.table(paste0(harvdir, "plink_covar_offspring"))
	} else {
		return()
	}

	colnames(ibd) = c("FID1", "IID1", "IID2", "PIHAT")
	colnames(pcs)[1:2] = c("FID", "SentrixID_1")

	corepaired = inner_join(corepaired, pcs[,2:(numPCs+2)], by="SentrixID_1") %>%
			anti_join(ibd, by=c("SentrixID_1"="IID1"))
	print(sprintf("after attaching covariates, %i remain", nrow(corepaired)))
	return(corepaired)
}


# For SNPTEST
# Usage: makeOutputs moms_height.txt AA87 moms
makeOutputs = function(corepaired, phenofile, phenoname, dataset){
	# output directory can be flipped for linear/survival models
	outdir = "~/Documents/harvest/linreg_files/"

	if(dataset=="moms"){
		samplelist = paste0(harvdir, "allmoms.txt")
	} else if (dataset=="fets"){
		samplelist = paste0(harvdir, "allfets.txt")
	} else {
		return()
	}
	write.table(corepaired, col.names=T, row.names=F, quote=F,
		file=paste0(outdir, phenofile, ".txt"))

	# make pheno
	processCmd = paste0("awk 'FNR==NR{p[$2]=$3 FS $4 FS $5 FS $6 FS $7 FS $8 FS $9 FS $10 FS $11; next} ",
			"FNR==1{print \"ID_1 ID_2 missing ", phenoname," SEX BATCH PC1 PC2 PC3 PC4 PC5 PC6\";",
			" print \"0 0 0 P D D C C C C C C\"} ",
			"$1 in p{print $1, $1, 0, p[$1]; next} ",
			"{print $1, $1, 0, \"NA NA NA NA NA NA NA NA NA\"}' ",
			outdir, phenofile, ".txt ",
			samplelist,
			" > ", outdir, phenofile, ".pheno")
	print(processCmd)
	system(processCmd)
	print("final pheno file is:")
	nr = system(paste0("wc -l ", outdir, phenofile, ".pheno"))

	# make snptest command
	outfile = paste0("../results_lin/", dataset, "_", phenoname)
	snptestcall = paste0("../snptest_v2.5.2 ",
 		 "-data ../../merging/", dataset, "_{}.vcf.gz ", phenofile, ".pheno -pheno ", phenoname,
		 " -use_raw_phenotypes -cov_all -genotype_field GP -method expected ",
		 "-frequentist 1 -o ", outfile, "_{}.txt &> ", outfile, "_{}.log")
	cmd = paste("seq 22 -1 1 | parallel -j $1 '", snptestcall, "'")
	cmdX = paste("echo 'X' | parallel -j 1 '", snptestcall, "'")
	writeChar(c("#!/bin/bash", cmd, cmdX),
 		  paste0(outdir, "snptest_", phenoname, ".sh"), eos = "\n")
}


################# PHENOTYPE-SPECIFIC PARTS ###################

### GA
m = read.table("~/data/mobaqs/p1724/harvest_mfr.csv", sep=";", h=T)
dim(m)
final = filter(m, FLERFODSEL==0,
	   DODKAT<6 | DODKAT>10,
	   ABRUPTIOP==0,
	   PLACENTA_PREVIA==0,
	   PREEKL_EKLAMPSI==0,
	   FOSTERV_POLYHYDRAMNION==0,
	   is.na(IVF),
	   is.na(DIABETES_MELLITUS),
	   HYPERTENSIV_TILSTAND==0 & HYPERTENSJON_KRONISK==0,
	   C00_MALF_ALL==0)
final = filter(final, SVLEN_DG<308)
nrow(final) 

attachPheno(final)
getCore(linkm)

## MATERNAL
corepaired = removeRepeated(coremoms, corekids)
corepaired = corepaired[,c("PREG_ID_1724", "SentrixID_1", "SVLEN_DG", "KJONN")]

# standard linear GA
phenoGA = attachCovariates(corepaired, "moms", 6)
makeOutputs(phenoGA, "lin_moms", "SVLEN_DG", "moms") 

# linear GA split below/above 39wk
phenoGAa39 = filter(phenoGA, SVLEN_DG>=273)
makeOutputs(phenoGAa39, "lin_moms_above39", "SVLEN_DG", "moms") 

phenoGAb39 = filter(phenoGA, SVLEN_DG<273)
makeOutputs(phenoGAb39, "lin_moms_below39", "SVLEN_DG", "moms") 

# PTD defined at 37 or 39 wk
phenoGA$above37 = as.numeric(phenoGA$SVLEN_DG>=259)
phenoGA$above39 = as.numeric(phenoGA$SVLEN_DG>=273)
phenoGA = select(phenoGA, -SVLEN_DG)

makeOutputs(select(phenoGA, -above39), "log_moms_ptd", "above37", "moms") 
makeOutputs(select(phenoGA, -above37), "log_moms_ptd39", "above39", "moms") 

## FETAL
corepaired = corekids
corepaired = corepaired[,c("PREG_ID_1724", "SentrixID_1", "SVLEN_DG", "KJONN")]

# standard linear GA
phenoGA = attachCovariates(corepaired, "fets", 6)
makeOutputs(phenoGA, "lin_fets", "SVLEN_DG", "fets") 

# linear GA split below/above 39wk
phenoGAa39 = filter(phenoGA, SVLEN_DG>=273)
makeOutputs(phenoGAa39, "lin_fets_above39", "above39", "fets") 

phenoGAb39 = filter(phenoGA, SVLEN_DG<273)
makeOutputs(phenoGAb39, "lin_fets_below39", "below39", "fets") 

# PTD defined at 37 or 39 wk
phenoGA$above37 = as.numeric(phenoGA$SVLEN_DG>=259)
phenoGA$above39 = as.numeric(phenoGA$SVLEN_DG>=273)
phenoGA = select(phenoGA, -SVLEN_DG)

makeOutputs(select(phenoGA, -above39), "log_fets_ptd", "ptd", "fets") 
makeOutputs(select(phenoGA, -above37), "log_fets_ptd39", "ptd39", "fets") 
