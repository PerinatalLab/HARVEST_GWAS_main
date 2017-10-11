#!/bin/bash

## This script is meant to provide an easy interface for calling
## ProbABEL's Cox PH module on VCF files.
## pacoxphVCF requires a raw .vcf file and an .mlinfo file,
## which are created by this script.
## Requires bcftools.
## For conditioned analysis, just prepare .pheno file with additional SNP column.

## runid argument allows easy parallelization and reuse of already generated VCFs.

## USAGE: ./run_pacoxphVCF.sh infile.vcf.gz phenofile outfile tmpdir chr:regionstart-regionstop runid

set -e

# paths for bcftools and pacoxphVCF
export PATH=/media/local-disk2/jjuod/erc-genotypes/bin:$PATH
export PATH=/media/local-disk2/jjuod/probabel/HARVEST_GWAS_main/bin:$PATH

echo "Checking input..."
if [ "$#" -ne 6 ]
then
	echo "Provided $# parameters instead of 6."
	echo "USAGE: ./run_pacoxphVCF.sh infile.vcf.gz phenofile outfile tmpdir chr:regionstart-regionstop runid"
	exit
fi

if [ ! -r $1 ]
then
	echo "Genotype file $1 does not exist or is not readable."
	exit
fi

if [ ! -r $2 ]
then
	echo "Phenotype file $2 does not exist or is not readable."
	exit
fi

echo "Checking .pheno file..."
touch ${4}/tempfile_${6}.diff
diff <(bcftools query -l $1) <(awk 'NR>1{print $1}' $2) > ${4}/tempfile_${6}.diff
if [ -s "tempfile.diff" ]
then
	echo "Error: found mismatches between provided .pheno file and VCF."
	echo "Mismatches listed in tempfile.diff."
	exit
fi

chr=${5%%:*}
posstart=${5#*:}
posstop=${posstart#*-}
posstart=${posstart%%-*}

echo "Checking region size..."
if [ $((posstop-posstart)) -gt 10000000 ]
then
	echo "Error: max allowed region size is 10 Mbp. You requested $((posstop-posstart))."
	exit
fi
echo "Working on chr ${chr}, region ${posstart}-${posstop}"

if [ -f ${4}/tempfile_${6}.mlinfo ]
then
	echo "Found already existing MLINFO file."
else
	echo "Creating .mlinfo file..."
	echo "rsid REF/NonEffAl ALT/EffAl POS AC INFO Rsq" > ${4}/tempfile_${6}.mlinfo
	bcftools query $1 \
		-f '%ID %REF %ALT %POS %INFO/AC %INFO/INFO 1\n' \
		-r ${chr}:${posstart}-${posstop} >> ${4}/tempfile_${6}.mlinfo
	echo "MLINFO file created."
fi

if [ -f ${4}/tempfile_${6}.vcf ]
then
	echo "Found already existing VCF file."
else
	echo "Extracting region ${chr}:${posstart}-${posstop} from the VCF file..."
	bcftools view $1 -r ${chr}:${posstart}-${posstop} -Ov -o ${4}/tempfile_${6}.vcf
	echo "Region extracted."
fi

echo "Calling ProbABEL..."
pacoxphVCF \
	-d ${4}/tempfile_${6}.vcf \
	-i ${4}/tempfile_${6}.mlinfo \
	-p $2 \
	-o ${3}_${6}
echo "Run complete."

