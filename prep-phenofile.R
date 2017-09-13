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
attachCovariates = function(corepaired, phenoname, dataset, numPCs){
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
	if(dataset=="moms"){
		samplelist = paste0(harvdir, "allmoms.txt")
	} else if (dataset=="fets"){
		samplelist = paste0(harvdir, "allfets.txt")
	} else {
		return()
	}
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


# For ProbABEL
# Usage: makeOutputs MotherSurvPheno moms
makeOutputsSurv = function(corepaired, phenofile, dataset){
	samplelist = read.table(paste0(harvdir, "all", dataset, ".txt"), h=F)

	samplelist = left_join(samplelist, corepaired, by=c("V1"="SentrixID_1"))
	print(table(samplelist[,3], deparse.level = 2, useNA="i"))

	colnames(samplelist)[1] = "id"
	write.table(samplelist, col.names=T, row.names=F, quote=F,
				file=paste0("~/Documents/harvest/phenofiles/", phenofile, ".txt"))

	# make IDs match probABEL's fvds:
	samplelist$id = 1:nrow(samplelist)
	write.table(samplelist, col.names=T, row.names=F, quote=F,
				file=paste0("~/Documents/harvest/phenofiles/", phenofile, ".pheno"))

	if(dataset=="fets"){
		samplelistx = read.table(paste0(harvdir, "allfetsx.txt"), h=F)
		samplelistx = left_join(samplelistx, corepaired, by=c("V1"="SentrixID_1"))

		colnames(samplelistx)[1] = "id"
		samplelistx$id = 1:nrow(samplelistx)
		write.table(samplelistx, col.names=T, row.names=F, quote=F,
					file=paste0("~/Documents/harvest/phenofiles/", phenofile, "X.pheno"))
	}

	print("final pheno file is:")
	nr = system(paste0("wc -l ~/Documents/harvest/phenofiles/", phenofile, ".pheno"))
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
getCore(linkm)
removeRepeated(coremoms, corekids)
corepaired <<- corepaired[,c("PREG_ID_1724", "SentrixID_1", "SVLEN_UL_DG", "BATCH")]
attachCovariates("SVLEN_UL_DG", "moms", 3)
makeOutputs("moms_dirty", "SVLEN_UL_DG", "moms") # 8879


### HEIGHT
q1 = read.table(pipe("cut -d ';' -f 1-100 ~/data/mobaqs/p1724/q1_pdb1724_v9.csv"), h=T, sep=";")
dim(q1)
q1 = filter(q1, AA87>100)
qplot(q1$AA87)
nrow(q1)

attachPheno(q1)
getCore(linkm)
removeRepeated(coremoms, corekids)
corepaired <<- corepaired[,c("PREG_ID_1724", "SentrixID_1", "AA87", "BATCH")]
attachCovariates("AA87", "moms", 3)
makeOutputs("moms_height", "AA87", "moms")


###### GA SURVIVAL ######

m = read.table("~/data/mobaqs/p1724/harvest_mfr.csv", sep=";", h=T)
dim(m)
final = filter(m, FLERFODSEL==0,
			   DODKAT<6 | DODKAT>10,
			   !is.na(SVLEN_DG))
final = filter(final, SVLEN_DG<308) %>%
	mutate(GAcor = SVLEN_DG-154)
final = mutate(final, PARITY0 = as.numeric(PARITET_5==0))
nrow(final)


## SPONT MAIN
phenospon = final
attachPheno(phenospon)
getCore(linkm)
# maternal
corepaired = removeRepeated(coremoms, corekids)
corepaired = mutate(corepaired, Spon = as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
						(is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
						INDUKSJON_PROSTAGLANDIN==0 & INDUKSJON_ANNET==0 &
						INDUKSJON_OXYTOCIN==0 & INDUKSJON_AMNIOTOMI==0))
corepaired = corepaired[,c("SentrixID_1", "GAcor", "Spon", "PARITY0")]
corepaired = attachCovariates(corepaired, "GAcor", "moms", 6) # 9196
makeOutputsSurv(corepaired, "MotherPhenoSpon", "moms")
# fetal
corepaired = corekids
corepaired = mutate(corepaired, Spon = as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
						(is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
						INDUKSJON_PROSTAGLANDIN==0 & INDUKSJON_ANNET==0 &
						INDUKSJON_OXYTOCIN==0 & INDUKSJON_AMNIOTOMI==0))
corepaired = corepaired[,c("SentrixID_1", "GAcor", "Spon", "PARITY0")]
corepaired = attachCovariates(corepaired, "GAcor", "fets", 6) # 9553
makeOutputsSurv(corepaired, "ChildPhenoSpon", "fets")

## PROM MAIN
# maternal
corepaired = removeRepeated(coremoms, corekids)
corepaired = mutate(corepaired, Prom = as.numeric(!is.na(VANNAVGANG)))
corepaired = corepaired[,c("SentrixID_1", "GAcor", "Prom", "PARITY0")]
corepaired = attachCovariates(corepaired, "GAcor", "moms", 6) # 9196
corepaired = makeOutputsSurv(corepaired, "MotherPhenoProm", "moms")
# fetal
corepaired = corekids
corepaired = mutate(corepaired, Prom = as.numeric(!is.na(VANNAVGANG)))
corepaired = corepaired[,c("SentrixID_1", "GAcor", "Prom", "PARITY0")]
corepaired = attachCovariates(corepaired, "GAcor", "fets", 6) # 9553
corepaired = makeOutputsSurv(corepaired, "ChildPhenoProm", "fets")


## SPONT SENS
phenosens = filter(final,
			   is.na(IVF),
			   ABRUPTIOP==0,
			   PLACENTA_PREVIA==0,
			   PREEKL_EKLAMPSI==0,
			   FOSTERV_POLYHYDRAMNION==0,
			   is.na(DIABETES_MELLITUS),
			   HYPERTENSIV_TILSTAND==0 & HYPERTENSJON_KRONISK==0,
			   C00_MALF_ALL==0)
nrow(phenosens)
attachPheno(phenosens)
getCore(linkm)

# maternal
corepaired = removeRepeated(coremoms, corekids)
corepaired = mutate(corepaired, Spon = as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
						(is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
						INDUKSJON_PROSTAGLANDIN==0 & INDUKSJON_ANNET==0 &
						INDUKSJON_OXYTOCIN==0 & INDUKSJON_AMNIOTOMI==0))
corepaired = corepaired[,c("SentrixID_1", "GAcor", "Spon", "PARITY0")]
corepaired = attachCovariates(corepaired, "GAcor", "moms", 6) # 8006
corepaired = makeOutputsSurv(corepaired, "MotherSensSpon", "moms")
# fetal
corepaired = corekids
corepaired = mutate(corepaired, Spon = as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
						(is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
						INDUKSJON_PROSTAGLANDIN==0 & INDUKSJON_ANNET==0 &
						INDUKSJON_OXYTOCIN==0 & INDUKSJON_AMNIOTOMI==0))
corepaired = corepaired[,c("SentrixID_1", "GAcor", "Spon", "PARITY0")]
corepaired = attachCovariates(corepaired, "GAcor", "fets", 6) # 8304
corepaired = makeOutputsSurv(corepaired, "ChildSensSpon", "fets")

## PROM SENS
# maternal
corepaired = removeRepeated(coremoms, corekids)
corepaired = mutate(corepaired, Prom = as.numeric(!is.na(VANNAVGANG)))
corepaired = corepaired[,c("SentrixID_1", "GAcor", "Prom", "PARITY0")]
corepaired = attachCovariates(corepaired, "GAcor", "moms", 6) # 8006
corepaired = makeOutputsSurv(corepaired, "MotherSensProm", "moms")
# fetal
corepaired = corekids
corepaired = mutate(corepaired, Prom = as.numeric(!is.na(VANNAVGANG)))
corepaired = corepaired[,c("SentrixID_1", "GAcor", "Prom", "PARITY0")]
corepaired = attachCovariates(corepaired, "GAcor", "fets", 6) # 8304
corepaired = makeOutputsSurv(corepaired, "ChildSensProm", "fets")


################# CONDITIONAL ANALYSES ###################

# recursion for greedy clumping
clump = function(df){
    best = filter(df, rank(desc(CHISQ), ties.method="first")==1)
    df = filter(df, CHR!=best$CHR | abs(POS-best$POS)>thr)
    if(nrow(df)>0){
        return(bind_rows(best, clump(df)))
    } else {
        return(best)
    }
}

# clump the snps at <1e-5 to loci
loci = data.frame()
allsuggestive = data.frame()
# max locus width
thr = 10000

for(genome in c("moms", "fets")){
    for(pheno in c("spon", "prom")){
        # read and store all suggestive snps
        infile = paste0("~/common/assoc_results/filtered/", genome, "_sens", pheno, "_suggestive.txt")
        infile = read.table(infile, h=T)
        infile$ID = 1:nrow(infile) + nrow(allsuggestive)
        allsuggestive = bind_rows(allsuggestive, infile)

        # clump them into regions
        locus = clump(infile)
        locus$genome = genome
        locus$pheno = pheno
        loci = bind_rows(loci, locus)
    }
}

# if no other analysis had any suggestive SNPs within 1Mb, no need to condition
thr2 = 1e6
hasNeighbors = inner_join(loci, allsuggestive, by="CHR") %>%
    filter(ID.x!=ID.y, abs(POS.x-POS.y)<thr2)
loci = filter(loci, ID %in% c(hasNeighbors$ID.x, hasNeighbors$ID.y))
write.table(loci, "~/common/assoc_results/conditional/loci_table.txt", col.names=T, row.names=F, quote=F)

source("~/Documents/gitrep/HARVEST_GWAS_main/read_vcf.R")


# function for attaching a column with another SNP dosage
# genomeCov = source of covariate (moms/fets)
# genomeTest = which will be tested (moms/fets)
# nr = primary key to identify files
# dfprom, dfspon = two dfs with phenotypes for the genome tested
attachGenotypeCov = function(chr, pos, genomeCov, genomeTest, nr, dfprom, dfspon){
	# attach either on same ID or relative ID:
	# same-genome conditioning
	if(genomeTest == genomeCov){
	    dfprom = left_join(dfprom, geno, by=c("SentrixID_1"="IID"))
	    dfspon = left_join(dfspon, geno, by=c("SentrixID_1"="IID"))
	} else {
	# opposite-genome testing
	    dfprom = left_join(dfprom, geno, by=c("relativeID"="IID"))
	    dfspon = left_join(dfspon, geno, by=c("relativeID"="IID"))
	}

    # remove the unnecessary column from final file
    dfprom = select(dfprom, -relativeID)
    dfspon = select(dfspon, -relativeID)

	# rename the new covariate to genome_chr_pos
	colnames(dfprom)[ncol(dfprom)] = paste(genomeCov, chr, pos, sep="_")
	colnames(dfspon)[ncol(dfspon)] = paste(genomeCov, chr, pos, sep="_")

	# filter IBD, attach PCs, make BATCH, save output
	if(genomeTest=="moms"){
		makeOutputsSurv(dfprom, paste0("MotherCondProm", nr), "moms")
		makeOutputsSurv(dfspon, paste0("MotherCondSpon", nr), "moms")
	} else {
		makeOutputsSurv(dfprom, paste0("ChildCondProm", nr), "fets")
		makeOutputsSurv(dfspon, paste0("ChildCondSpon", nr), "fets")
	}
}


## now using sens as main analyses
## same exclusions for both SPONT and PROM:
phenosens = filter(final,
				   is.na(IVF),
				   ABRUPTIOP==0,
				   PLACENTA_PREVIA==0,
				   PREEKL_EKLAMPSI==0,
				   FOSTERV_POLYHYDRAMNION==0,
				   is.na(DIABETES_MELLITUS),
				   HYPERTENSIV_TILSTAND==0 & HYPERTENSJON_KRONISK==0,
				   C00_MALF_ALL==0)
nrow(phenosens)
attachPheno(phenosens)
getCore(linkm)

# make id table for conditioning mothers on fetal genome & v.v.
corepaired = removeRepeated(coremoms, corekids)
pairtable = full_join(corepaired, corekids, by="PREG_ID_1724")
pairtable = pairtable[, c("SentrixID_1.x", "SentrixID_1.y")]
colnames(pairtable) = c("momID", "kidID")

# common prep - maternal
corepaired = mutate(corepaired, Prom = as.numeric(!is.na(VANNAVGANG)),
                   Spon = as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                            (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                            INDUKSJON_PROSTAGLANDIN==0 & INDUKSJON_ANNET==0 &
                            INDUKSJON_OXYTOCIN==0 & INDUKSJON_AMNIOTOMI==0))
coreMomsSpon = corepaired[,c("SentrixID_1", "GAcor", "Spon", "PARITY0")]
coreMomsProm = corepaired[,c("SentrixID_1", "GAcor", "Prom", "PARITY0")]
coreMomsSpon = attachCovariates(coreMomsSpon, "GAcor", "moms", 6)
coreMomsProm = attachCovariates(coreMomsProm, "GAcor", "moms", 6)

# common prep - fetal
corepaired = mutate(corekids, Prom = as.numeric(!is.na(VANNAVGANG)),
                    Spon = as.numeric(FSTART==1 & (is.na(KSNITT) | KSNITT>1) &
                              (is.na(KSNITT_PLANLAGT) | KSNITT_PLANLAGT==1) &
                              INDUKSJON_PROSTAGLANDIN==0 & INDUKSJON_ANNET==0 &
                              INDUKSJON_OXYTOCIN==0 & INDUKSJON_AMNIOTOMI==0))
coreKidsSpon = corepaired[,c("SentrixID_1", "GAcor", "Spon", "PARITY0")]
coreKidsProm = corepaired[,c("SentrixID_1", "GAcor", "Prom", "PARITY0")]
coreKidsSpon = attachCovariates(coreKidsSpon, "GAcor", "fets", 6)
coreKidsProm = attachCovariates(coreKidsProm, "GAcor", "fets", 6)


# attach the other family member ID
coreMomsSpon = left_join(coreMomsSpon, pairtable, by=c("SentrixID_1"="momID"))
coreMomsProm = left_join(coreMomsProm, pairtable, by=c("SentrixID_1"="momID"))
coreKidsSpon = left_join(coreKidsSpon, pairtable, by=c("SentrixID_1"="kidID"))
coreKidsProm = left_join(coreKidsProm, pairtable, by=c("SentrixID_1"="kidID"))
ncm = ncol(coreMomsSpon)
ncf = ncol(coreKidsSpon)
colnames(coreMomsSpon)[ncm] = colnames(coreMomsProm)[ncm] = "relativeID"
colnames(coreKidsSpon)[ncf] = colnames(coreKidsProm)[ncf] = "relativeID"

## For each genotype, we automatically generate all pairs (moms/fets x prom/spon),
## so this now includes:
## Stepwise approach to identify independent SNPs within locus ("WITHIN")
## Identifying maternal vs. fetal action of SNPs by conditioning on the opposite genome ("ACROSS")

loci = read.table("~/common/assoc_results/conditional/loci_table.txt", h=T)
loci$CHR[loci$CHR==23] = "X"

# generate files
for(l in seq_along(loci$CHR)){
    print(sprintf("working on SNP %i / %i", l, nrow(loci)))

    # read in genotype covariate
    geno = readVcf(loci$CHR[l], loci$POS[l], NA, loci$genome[l])[c("IID", "DS")]
    num = loci$ID[l]

    # attach it to moms prom & spon
    attachGenotypeCov(loci$CHR[l], loci$POS[l], loci$genome[l], "moms", num, coreMomsProm, coreMomsSpon)
    # attach it to fets prom & spon
    attachGenotypeCov(loci$CHR[l], loci$POS[l], loci$genome[l], "fets", num, coreKidsProm, coreKidsSpon)
}

