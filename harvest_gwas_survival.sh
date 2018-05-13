## This script is for documenting the process of HARVEST survival analyses.
## It is not intended to be run as-is, but it should be an exact description
## of the entire pipeline.

### PREPARING FVD FILES ###
# vcf (in ../merging) -> gen (BKP_maternal)
./vcftogen.sh

# for moms, second half of chr 2 gets NAN results, and chr 16 breaks down around middle
# therefore, these two chrs are split into chunks for moms:
# (includes gen -> fvd)
rm moms_16.* moms_2.*
./vcftogen_chunk.sh

# gen (BKP_maternal) -> fvd (BKP_maternal)
# somehow, sample IDs still get replaced with 1:nrow,
# but the order in .fvd matches ../merging/all*txt.
for c in {1..22} "X"
do
	./gentofvi_moms.sh ${c}
	./gentofvi_fets.sh ${c}
done


### PREPARING PHENOTYPE ###
# using the section for survival models:
prep-phenofile.R
# copy the generated pheno files to phenofiles/


### SETTING UP FOLDERS ###
# all genotype files are then distributed to a number of folders
# for easy parallelization:
for c in {1..22}
do
	dest=$((chr % 6))
	mkdir -p maternalfiles${dest}
	mkdir -p fetalfiles${dest}
	cp moms_${c}.* maternalfiles${dest}/
	cp fets_${c}.* fetalfiles${dest}/
done
mkdir -p maternalfiles6/
mkdir -p fetalfiles6/
# fetal X will need separate pheno files:
# (not maternal, because their samples for auto and X are same)
cp fets_X.* fetalfiles6/
# move chunked chrs as well:
cp moms_X.* maternalfiles6/
cp moms_2_*.* maternalfiles6/
cp moms_16_*.* maternalfiles6/


### MAIN GWAS RUN ###
seq 0 6 | parallel -j 7 ./run_survival_maternal.R {}
seq 0 5 | parallel -j 6 ./run_survival_fetal.R {} phenofiles
./run_survival_fetal.R 6 phenofilesX

### CONDITIONAL runs within same genome ###
# select interesting SNPs at below 1e-5:
awk '$5>19.51 && $4>=0.01' fets_sensprom_all.txt > filtered/fets_sensprom_suggestive.txt # etc.

# make conditional phenos:
prep-phenofile.R
# delete X files
rm ~/Documents/harvest/phenofiles/*Cond*X.pheno
# upload to HUNT
cp ~/Documents/harvest/phenofiles/*Cond*.pheno ...

# run:
# manually delete one Xchr snp from loci_table.txt > loci_table_nox.txt
./run_conditional.sh phenofiles_cond/loci_table_nox.txt 1000000 tmp
# get results:
./check_conditional.R
