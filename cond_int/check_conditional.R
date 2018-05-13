#!/usr/bin/Rscript

# This is a basic script to quickly glance the results
# of conditional analyses before plotting.

options(stringsAsFactors=F)
library(dplyr)
library(tidyr)

loci = read.table("/media/local-disk2/jjuod/probabel/phenofiles_cond/loci_table_nox.txt", h=T)
pfull = rfull = data.frame()

for(l in seq_along(loci$CHR)){
	print(sprintf("working on file %i / %i", l, nrow(loci)))
	# read in all 4 runs
	id = loci$ID[l]
	indir = "/media/local-disk2/jjuod/probabel/results/cond/res_"
	infile1 = read.table(paste0(indir, "moms_spon_", id, "m_add.out.txt"), h=T)
	infile2 = read.table(paste0(indir, "moms_prom_", id, "m_add.out.txt"), h=T)
	infile3 = read.table(paste0(indir, "fets_spon_", id, "f_add.out.txt"), h=T)
	infile4 = read.table(paste0(indir, "fets_prom_", id, "f_add.out.txt"), h=T)
	infiles = bind_rows(moms_spon=infile1, moms_prom=infile2,
			fets_spon=infile3, fets_prom=infile4, .id="analysis") %>%
		filter(Mean_predictor_allele>=0.01, Mean_predictor_allele<=0.99)
	
	# discard crap and get P-values from chi^2
	infiles = group_by(infiles, analysis) %>%
		filter(rank(desc(chi2_SNP_add), ties.method="first")==1)
	infiles$Pnew = pchisq(infiles$chi2_SNP_add, df=1, ncp=0, lower.tail=F)

	pvals = infiles[,c("analysis", "Pnew")]
	pvals$Pnew = sprintf("%0.3g", pvals$Pnew)
	rsids = infiles[,c("analysis", "name")]
	pvals$analysis = paste("P", pvals$analysis, sep="_")
	rsids$analysis = paste("SNP", rsids$analysis, sep="_")
	
	# join all conditionals together
	pvals$ID = id
	rsids$ID = id
	pfull = spread(pvals, key=analysis, value=Pnew) %>%
		bind_rows(pfull)
	rfull = spread(rsids, key=analysis, value=name) %>%
		bind_rows(rfull)
}

# merge with old results
loci2 = left_join(loci, pfull, by="ID") %>%
	left_join(rfull, by="ID")
write.table(loci2, "/media/local-disk2/jjuod/probabel/phenofiles_cond/topres_nox.txt", quote=F, col.names=T, row.names=F)
