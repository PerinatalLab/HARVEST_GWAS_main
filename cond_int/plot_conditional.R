#!/usr/bin/Rscript

# This script shall read and plot conditional analysis results.
# USAGE: ./script.R resdir

options(stringsAsFactors = F)
library(dplyr)
library(ggplot2, lib.loc = .libPaths()[1])
args = commandArgs(TRUE)
resdir = args[1]
setwd(resdir)

filelist = strsplit(list.files(pattern="snplist.*0"), "_")
genos = sapply(filelist, "[[", 2)
nums = sapply(filelist, "[[", 3)

for(i in seq_along(genos)){
	print(sprintf("working on SNP %i / %i", i, length(genos)))
	
	# read in conditioned results
	newSp0 = read.table(paste("res", genos[i], "spon", nums[i], 0, "add.out.txt", sep="_"), h=T)
	newSp1 = read.table(paste("res", genos[i], "spon", nums[i], 1, "add.out.txt", sep="_"), h=T)
	newPr0 = read.table(paste("res", genos[i], "prom", nums[i], 0, "add.out.txt", sep="_"), h=T)
	newPr1 = read.table(paste("res", genos[i], "prom", nums[i], 1, "add.out.txt", sep="_"), h=T)
	print("new results scanned")
	
	# read in mlinfos and replace rounded positions
	snplist0 = read.table(paste("snplist", genos[i], nums[i], 0, sep="_"), h=T)
	snplist1 = read.table(paste("snplist", genos[i], nums[i], 1, sep="_"), h=T)
	newSp0$POS = snplist0$POS
	newSp1$POS = snplist1$POS
	newPr0$POS = snplist0$POS
	newPr1$POS = snplist1$POS
	print("SNP positions scanned")
	
	# identify the initial snp
	hit = newSp1$POS[1]
	
	# read the r2 values for hit snp
	ldf0 = paste0("zcat ld_", genos[i], "_", nums[i], "_0.ld.gz | awk '$5==", hit, " || NR==1'")
	ldf1 = paste0("zcat ld_", genos[i], "_", nums[i], "_1.ld.gz | awk '$2==", hit, " || NR==1'")
	ld0 = read.table(pipe(ldf0), h=T)[,c("BP_A", "R2")]
	ld1 = read.table(pipe(ldf1), h=T)
	rsid = ld1$SNP_A[1]
	chr = ld1$CHR_A[1]
	ld1 = ld1[,c("BP_B", "R2")]
	print("R2 values scanned")
	
	# attach corrs to the new results
	newSp0 = left_join(newSp0, ld0, by=c("POS"="BP_A"))
	newSp1 = left_join(newSp1, ld1, by=c("POS"="BP_B"))
	newPr0 = left_join(newPr0, ld0, by=c("POS"="BP_A"))
	newPr1 = left_join(newPr1, ld1, by=c("POS"="BP_B"))
	
	# cut away overlaps, join both half-files
	newSp = bind_rows(newSp0, newSp1[-1,])
	newPr = bind_rows(newPr0, newPr1[-1,])
	
	# missing R2 are <0.2 (not printed in PLINK's output)
	newSp$R2[is.na(newSp$R2)] = 0
	newPr$R2[is.na(newPr$R2)] = 0
	
	# convert new X^2 into p-values
	newSp$Pnew = pchisq(newSp$chi2_SNP_add, df=1, ncp=0, lower.tail=F)
	newPr$Pnew = pchisq(newPr$chi2_SNP_add, df=1, ncp=0, lower.tail=F)
	
	# read in old results for same SNPs, excluding MAF<0.01
	posrange = range(newSp$POS)
	oldSpf = paste0("../", genos[i], "_sensspon_small.txt")
	oldPrf = paste0("../", genos[i], "_sensprom_small.txt")
	comm = paste("awk '$1==", chr, "&& $2>=", posrange[1], "&& $2<=", posrange[2], "&& $3>=0.01'", oldSpf)
	print(comm)
	oldSp = read.table(pipe(comm), h=F)
	
	comm = paste("awk '$1==", chr, "&& $2>=", posrange[1], "&& $2<=", posrange[2], "&& $3>=0.01'", oldPrf)
	print(comm)
	oldPr = read.table(pipe(comm), h=F)
	colnames(oldSp) = colnames(oldPr) = c("CHR", "POS", "MAF", "Pold")
	print("old results scanned")
	
	# merge raw and conditioned results
	# (not merging on position because of multiallelic loci)
	newSp = filter(newSp, Mean_predictor_allele>=0.01, Mean_predictor_allele<=0.99)
	newPr = filter(newPr, Mean_predictor_allele>=0.01, Mean_predictor_allele<=0.99)
	if(any(newSp$POS != oldSp$POS)){
		print("ERROR: old and new SPON result positions differ")	
	} else {
		newSp$Pold = oldSp$Pold
	}
	if(any(newPr$POS != oldPr$POS)){
		print("ERROR: old and new PROM result positions differ")
	} else {
		newPr$Pold = oldPr$Pold
	}
	
	
	# for the hit SNP, assign R2=1 and P=oldP
	hitpSp = hitpPr = NA
	if(sum(newSp$POS == hit) > 1){
		print("ERROR: hit SNP position appears duplicated")	
	} else {
		newSp$R2[newSp$POS == hit] = 1
		hitpSp = newSp$Pold[which(newSp$POS == hit)]
		newSp$Pnew[which(newSp$POS == hit)] = hitpSp
		
		newPr$R2[newPr$POS == hit] = 1
		hitpPr = newPr$Pold[which(newPr$POS == hit)]
		newPr$Pnew[which(newPr$POS == hit)] = hitpPr
	}
	
	# plot
	pSp = arrange(newSp, R2) %>%
		ggplot(aes(x=POS)) + 
		geom_point(aes(y=-log(Pold, 10)), shape=1, col="azure4") +
		geom_point(aes(y=-log(Pnew, 10), col=R2), size=0.7) +
		labs(title = paste("SNP", rsid, "- position", chr, ":", hit, "- SPON -", genos[i]),
			 subtitle = sprintf("previous p-value: %.3g", hitpSp),
			 y = "-log10 p-value") +
		theme_bw()
	ggsave(paste("plot", genos[i], nums[i], "spon.png", sep="_"),
		   plot = pSp, width=7, height=4)
	
	pPr = arrange(newPr, R2) %>%
		ggplot(aes(x=POS)) + 
		geom_point(aes(y=-log(Pold, 10)), shape=1, col="azure4") +
		geom_point(aes(y=-log(Pnew, 10), col=R2), size=0.7) +
		labs(title = paste("SNP", rsid, "- position ", chr, ":", hit, "- PROM -", genos[i]),
			 subtitle = sprintf("previous p-value: %.3g", hitpPr),
			 y = "-log10 p-value") +
		theme_bw()
	ggsave(paste("plot", genos[i], nums[i], "prom.png", sep="_"),
		   plot = pPr, width=7, height=4)
}
