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

# Usage: getCore()
getCore = function(){
	coremoms <<- filter(flags, genotypesOK, phenotypesOK, coreOK) %>%
		semi_join(linkm, ., by=c("SentrixID_1"="IID")) %>%
		filter(Role=="Mother")
	# 7065, incl. multiple pregs
	print(sprintf("%i core mothers found", nrow(coremoms)))
	
	corekids <<- filter(flags, genotypesOK, phenotypesOK, coreOK) %>%
		semi_join(linkm, ., by=c("SentrixID_1"="IID")) %>%
		filter(Role=="Child")
	print(sprintf("%i core kids found", nrow(corekids)))
}

# removes repeated pregnancies,
# keeping the one included in the other core set.
# Usage: removeRepeated(coremoms, corekids) to clean mothers
# 		removeRepeated(corekids, coremoms) to clean kids
removeRepeated = function(dftoclean, basedon){
	corepaired = group_by(dftoclean, SentrixID_1) %>%
		mutate(hasPair = PREG_ID_1724 %in% basedon$PREG_ID_1724) %>%
		top_n(1, hasPair) %>%
		filter(rank(PREG_ID_1724)==1)
	# 6836, no more multiple pregs
	print(sprintf("after removing repeated pregnancies, %i remain", nrow(corepaired)))
	
	## add batch covariate
	corepaired <<- select(flags, one_of(c("IID", "BATCH"))) %>%
		left_join(corepaired, ., by=c("SentrixID_1" = "IID")) %>%
		mutate(BATCH = as.numeric(BATCH=="M24"))
}

# Usage: makeOutputs moms_height.txt AA87 moms
makeOutputs = function(phenofile, phenoname, dataset){
	corepaired = corepaired[,c("PREG_ID_1724", "SentrixID_1", phenoname, "BATCH")]

	if(dataset=="moms"){
		ibd = read.table(paste0(harvdir, "ibd_pihat_exclude_mothers"))
		pcs = read.table(paste0(harvdir, "plink_covar_mothers"))
		samplelist = paste0(harvdir, "allmoms.txt")
	} else if (dataset=="kids") {
		ibd = read.table(paste0(harvdir, "ibd_pihat_exclude_offspring"))
		pcs = read.table(paste0(harvdir, "plink_covar_offspring"))
		samplelist = paste0(harvdir, "allkids.txt")
	} else {
		return()
	}
	
	colnames(ibd) = c("FID1", "IID1", "IID2", "PIHAT")
	colnames(pcs)[1:2] = c("FID", "SentrixID_1")
	
	print(head(corepaired))
	corepaired = inner_join(corepaired, pcs[,2:5], by="SentrixID_1") %>%
		anti_join(ibd, by=c("PREG_ID_1724"="FID1", "SentrixID_1"="IID1"))
	
	print(sprintf("after attaching covariates, %i remain", nrow(corepaired)))
	write.table(corepaired, col.names=T, row.names=F, quote=F,
				file=paste0("~/Documents/harvest/phenofiles/", phenofile, ".txt"))
	
	# make pheno
	processCmd = paste0("awk 'FNR==NR{p[$2]=$3 FS $4 FS $5 FS $6 FS $7; next} ",
						"FNR==1{print \"ID_1 ID_2 missing ", phenoname," BATCH PC1 PC2 PC3\";",
							" print \"0 0 0 P D C C C\"} ",
						"$1 in p{print $1, $1, 0, p[$1]; next} ",
						"{print $1, $1, 0, \"NA NA NA NA NA\"}' ",
						"~/Documents/harvest/phenofiles/", phenofile, ".txt ",
						samplelist,
						"> ~/Documents/harvest/phenofiles/", phenofile, ".pheno")
	print(processCmd)
	system(processCmd)
	print("final pheno file is:")
	nr = system(paste0("wc -l ~/Documents/harvest/phenofiles/", phenofile, ".pheno"))
	
	# make snptest command
	outfile = paste0("results/", dataset, "_", phenoname)
	snptestcall = paste0("./snptest_v2.5.2 ",
		 "-data ../merging/", dataset, "_{}.vcf.gz ", phenofile, ".pheno -pheno ", phenoname,
		 " -use_raw_phenotypes -cov_all -genotype_field GP -method expected ",
		 "-frequentist 1 -o ", outfile, "_{}.txt &> ", outfile, "_{}.log")
	cmd = paste("seq 22 -1 1 | parallel -j $1 '", snptestcall, "'")
	cmdX = paste("echo 'X' | parallel -j 1 '", snptestcall, "'")
	writeChar(c("#!/bin/bash", cmd, cmdX),
			  paste0("~/Documents/harvest/snptest_", phenoname, ".sh"), eos = "\n")
}

################# PHENOTYPE-SPECIFIC PARTS ###################

### GA
m = read.table("~/data/mobaqs/p1724/harvest_mfr.csv", sep=";", h=T)
dim(m)
final = filter(m, FLERFODSEL==0,
			   FSTART==1 & (is.na(KSNITT) | KSNITT>1),
			   is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1,
			   INDUKSJON_PROSTAGLANDIN==0 & INDUKSJON_ANNET==0 &
			   	INDUKSJON_OXYTOCIN==0 & INDUKSJON_AMNIOTOMI==0,
			   DODKAT<6 | DODKAT>10,
			   ABRUPTIOP==0,
			   PLACENTA_PREVIA==0,
			   PREEKL_EKLAMPSI==0,
			   FOSTERV_POLYHYDRAMNION==0,
			   is.na(IVF),
			   is.na(DIABETES_MELLITUS),
			   HYPERTENSIV_TILSTAND==0 & HYPERTENSJON_KRONISK==0,
			   C00_MALF_ALL==0)
nrow(final)

final = filter(final, SVLEN_UL_DG<308)
nrow(final) # 8201

corepaired = corepaired[,c("PREG_ID_1724", "SentrixID_1", "SVLEN_UL_DG")]
write.table(corepaired, "~/Documents/harvest/linreg_files/coremoms.txt",
			col.names=T, row.names=F, quote=F)

### GA DIRTY
m = read.table("~/data/mobaqs/p1724/harvest_mfr.csv", sep=";", h=T)
dim(m)
final = filter(m, FLERFODSEL==0,
			   DODKAT<6 | DODKAT>10,
			   is.na(IVF))
nrow(final)

final = filter(final, SVLEN_UL_DG<308)
nrow(final) # 10909

attachPheno(final)
getCore()
removeRepeated(coremoms, corekids)
makeOutputs("moms_dirty", "SVLEN_UL_DG", "moms") # 8879


### HEIGHT
q1 = read.table(pipe("cut -d ';' -f 1-100 ~/data/mobaqs/p1724/q1_pdb1724_v9.csv"), h=T, sep=";")
dim(q1)
q1 = filter(q1, AA87>100)
qplot(q1$AA87)
nrow(q1)

attachPheno(q1)
getCore()
removeRepeated(coremoms, corekids)
makeOutputs("moms_height", "AA87", "moms")
