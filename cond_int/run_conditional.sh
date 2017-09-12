#!/bin/bash

# This script will launch ./run_pacoxphVCF.sh on provided phenofiles,
# to allow convenient run of conditional analyses.

# USAGE: ./run_conditional.sh phenodir range

set -e

INDIR=/media/local-disk2/jjuod/merging
RESDIR=/media/local-disk2/jjuod/probabel/results
for file in $1/*.txt
do
	echo "Working on file ${file}"
	
	# trim filename to get only the number	
	num=${file%%.txt}
	num=${num##*Prom}
	num=${num##*Spon}
	if [[ $file == *"Mother"* ]]
	then
		echo "Analyzing maternal genomes"
		genome="moms"
	else
		echo "Analyzing fetal genomes"
		genome="fets"
	fi

	# check 3rd column for Prom/Spon
	pheno=$(awk 'NR==1{print tolower($3)}' $file)
	echo "Working on phenotype ${pheno}"

	# check 4th position to identify region
	pos=( $(awk 'NR==1{split($5, a, "_"); print a[2], a[3]}' $file) )
	chr=${pos[0]}
	pstart=$(( ${pos[1]} - $2 ))
	pstop=$(( ${pos[1]} + $2 ))
	echo "Analyzing chr ${chr}, region ${pstart}-${pstop}"

	# main run
	./run_pacoxphVCF.sh \
		${INDIR}/${genome}_${chr}.vcf.gz \
		${file} \
		${RESDIR}/cond/res_${genome}_${pheno}_${num} \
		${chr}:${pstart}-${pstop}
	echo "Analysis complete."
	cp tempfile_0.mlinfo ${RESDIR}/cond/snplist_${genome}_${num}_0
	cp tempfile_1.mlinfo ${RESDIR}/cond/snplist_${genome}_${num}_1

	# run PLINK to get regional LD for easier plotting later
	if [ ! -f ${RESDIR}/cond/ld_${genome}_${num} ]
	then
		echo "Generating LD matrix"
		plink \
			--vcf tempfile_0.vcf \
			--r2 inter-chr gz \
			--out ${RESDIR}/cond/ld_${genome}_${num}_0
		echo "Generating LD matrix"
		plink \
			--vcf tempfile_1.vcf \
			--r2 inter-chr gz \
			--out ${RESDIR}/cond/ld_${genome}_${num}_1
	else
		echo "LD matrix already exists."
	fi
done
