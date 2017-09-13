#!/bin/bash

# This script will launch ./run_pacoxphVCF.sh on provided phenofiles,
# to allow convenient run of conditional analyses.

# USAGE: ./run_conditional.sh phenolist range tmpdir

set -e

INDIR=/media/local-disk2/jjuod/merging
RESDIR=/media/local-disk2/jjuod/probabel/results
PHENODIR=/media/local-disk2/jjuod/probabel/phenofiles_cond
TMPDIR=$3
RANGE=$2

# 1. run ID
# 2. chromosome
# 3. position (index SNP)
function startrun {
	echo "Working on run ID ${1}"
	scriptdir=/media/local-disk2/jjuod/probabel/HARVEST_GWAS_main/cond_int

	# get lead SNP and identify region
	pstart=$(( ${3} - $RANGE ))
	pstop=$(( ${3} + $RANGE ))
	echo "Analyzing chr ${2}, region ${pstart}-${pstop}"

	# main runs
	# note the use of ${id}m as run indicator,
	# to keep mom/fet tmp files separate
	${scriptdir}/run_pacoxphVCF.sh \
		${INDIR}/moms_${2}.vcf.gz \
		${PHENODIR}/MotherCondSpon${1}.txt \
		${RESDIR}/cond/res_moms_spon ${TMPDIR} \
		${2}:${pstart}-${pstop} ${1}m
	echo "Analysis moms+spon complete."
	${scriptdir}/run_pacoxphVCF.sh \
		${INDIR}/moms_${2}.vcf.gz \
		${PHENODIR}/MotherCondProm${1}.txt \
		${RESDIR}/cond/res_moms_prom ${TMPDIR} \
		${2}:${pstart}-${pstop} ${1}m
	echo "Analysis moms+prom complete."
	${scriptdir}/run_pacoxphVCF.sh \
		${INDIR}/fets_${2}.vcf.gz \
		${PHENODIR}/ChildCondSpon${1}.txt \
		${RESDIR}/cond/res_fets_spon ${TMPDIR} \
		${2}:${pstart}-${pstop} ${1}f
	echo "Analysis fets+spon complete."
	${scriptdir}/run_pacoxphVCF.sh \
		${INDIR}/fets_${2}.vcf.gz \
		${PHENODIR}/ChildCondProm${1}.txt \
		${RESDIR}/cond/res_fets_prom ${TMPDIR} \
		${2}:${pstart}-${pstop} ${1}f
	echo "Analysis fets+prom complete."

	# save snplists to get unrounded positions
	cp ${TMPDIR}/tempfile_${1}m.mlinfo ${RESDIR}/cond/snplist_moms_${1}
	cp ${TMPDIR}/tempfile_${1}f.mlinfo ${RESDIR}/cond/snplist_fets_${1}

	# run PLINK to get regional LD for easier plotting later
	echo "Generating LD matrix for moms"
	plink \
		--vcf ${TMPDIR}/tempfile_${1}m.vcf \
		--r2 inter-chr gz \
		--out ${RESDIR}/cond/ld_moms_${1}
	echo "Generating LD matrix for fets"
	plink \
		--vcf ${TMPDIR}/tempfile_${1}f.vcf \
		--r2 inter-chr gz \
		--out ${RESDIR}/cond/ld_fets_${1}
	rm ${TMPDIR}/tempfile_${1}m.*
	rm ${TMPDIR}/tempfile_${1}f.*
}

export -f startrun
export INDIR
export PHENODIR
export TMPDIR
export RESDIR
export RANGE

# Requires a space-delimited file with header,
# listing all SNPs to be processed.
parallel -j 40 --header : --colsep ' ' startrun {ID} {CHR} {POS} < $1

