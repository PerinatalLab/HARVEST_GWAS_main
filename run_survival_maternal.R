#!/usr/bin/Rscript

# USAGE: ./run_survival_maternal.R NUMDIR PHENODIR
# where NUMDIR is directory suffix for parallelizing

library(GenABEL)

args = commandArgs(TRUE)
#Survival data in the Survival data folder
pathpheno<-paste0("/media/local-disk2/jjuod/probabel/", args[2])
#Genotype in the imputvcf folder
pathgeno<-paste0("/media/local-disk2/jjuod/probabel/maternalfiles", args[1])
 
maternal=TRUE
 
#Phenos
if(maternal){
	#Mom
	#List of the files names in the pheno path ,starting with M finishing with .txt
	listphenofiles<-list.files(path=pathpheno,pattern = "^M.*.pheno$")
	#CHR
	#List of the files names in the pheno path finishing with .mlinfo, unzipped
	listmlinfo<-list.files(path=pathgeno,pattern = "^m.*mlinfo")
	# index files accompanying .fvd, unzipped
	listfvi<-list.files(path=pathgeno,pattern = "^m.*fvi")
	# actual genotype files, gzipped
	listfvd<-list.files(path=pathgeno,pattern = "^m.*fvd.gz")
} else {		 
	#Child
	#List of the files names in the pheno path ,starting with C finishing with .txt
	listphenofiles<-list.files(path=pathpheno,pattern = "^C.*.pheno$")
	listfvi<-list.files(path=pathpheno,pattern = "^c.*fvi$")
	listmlinfo<-list.files(path=pathpheno,pattern = "^c.*mlinfo$")
	listfvd<-list.files(path=pathpheno,pattern = "^c.*fvd.gz")
}
write(listphenofiles, file="/media/local-disk2/jjuod/probabel/results/phenolist.txt")
 
# Need to have the same order in listfvd, listfvi and listmlinfo
listmlinfo<-sort(listmlinfo)
listfvi<-sort(listfvi)
listfvd<-sort(listfvd)
 
for( i in 1:length(listmlinfo))
{
	# set file paths:
	myfvi<-paste(pathgeno, listfvi[i], sep="/")
	mymlinfo<-paste(pathgeno, listmlinfo[i], sep="/")
	genofile<-paste(pathgeno, listfvd[i], sep="/")
	myres = paste("/media/local-disk2/jjuod/probabel/results/res",
		gsub(pattern=".dose.fvd.gz",replacement="",listfvd[i]), sep ="_")
 
	# unzip:
	print(paste("Unzipping",genofile))
	cmdline<-paste("gunzip -k",genofile)
	print(cmdline)
	system(cmdline)
 
	# run analysis:
	# ProAbel Cox command, for score test
	cmdline<-paste("parallel -j ", length(listphenofiles),
			" ./bin/pacoxph --pheno ", pathpheno, "/{} --dose ", myfvi,
			" --info ", mymlinfo, " --out ", myres ,"_{.} < ../results/phenolist.txt", sep="")
	# System only work with Linux
	print(cmdline)
	system(cmdline)
 
	# remove the unzipped fvd files
	print(paste("Removing the unzipped",genofile,sep=" "))
	cmdline<-paste("rm -I", gsub(pattern=".gz",replacement="",genofile), sep = " ")
	print(cmdline)
	system(cmdline)
}
 

