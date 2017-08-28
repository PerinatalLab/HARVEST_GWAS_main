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
prepare_main_phenos.R

### SETTING UP FOLDERS ###

### MAIN GWAS RUN ###


