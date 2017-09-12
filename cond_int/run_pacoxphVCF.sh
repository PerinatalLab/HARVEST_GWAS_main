#!/bin/bash

## This script is meant to provide an easy interface for calling
## ProbABEL's Cox PH module on VCF files.
## pacoxphVCF requires a raw .vcf file and an .mlinfo file,
## which are created by this script.
## Requires bcftools.
## For conditioned analysis, just prepare .pheno file with additional SNP column.
## TODO: interaction analysis.

## USAGE: ./run_pacoxphVCF.sh infile.vcf.gz phenofile outfile chr:regionstart-regionstop [SNPtoconditionon]

set -e

export PATH=/media/local-disk2/jjuod/erc-genotypes/bin:$PATH

echo "Checking input..."
if [ "$#" -ne 5 && "$#" -ne 6 ]
then
	echo "Provided $# parameters instead of 5 or 6."
	echo "USAGE: ./run_pacoxphVCF.sh infile.vcf.gz phenofile outfile chr:regionstart-regionstop [SNPtoconditionon]"
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
touch tempfile.diff
diff <(bcftools query -l $1) <(awk 'NR>1{print $1}' $2) > tempfile.diff
if [ ! -s "tempfile.diff" ]
then
	echo "Error: found mismatches between provided .pheno file and VCF."
	echo "Mismatches listed in tempfile.diff."
	exit
fi

chr=${4%%:*}
posstart=${4#*:}
posstop=${posstart%%-*}
posstart=${posstart#*-}

echo "Checking region size..."
if [ $((posstop-posstart)) -gt 10000000 ]
then
	echo "Error: max allowed region size is 10 Mbp. You requested $((posstop-posstart))."
	exit
fi

echo "Creating .mlinfo file..."
echo "rsid REF/NonEffAl ALT/EffAl POS AC INFO Rsq" > tempfile.mlinfo
bcftools query -f '%ID %REF %ALT %POS %INFO/AC %INFO/INFO 1\n' $1 >> tempfile.mlinfo
echo "MLINFO file created."

echo "Extracting region ${chr}:${posstart}-${posstop} from the VCF file..."
bcftools view $1 -r ${chr}:${posstart}-${posstop} -Ov -o tempfile.vcf
echo "Region extracted."

echo "Calling ProbABEL..."
./bin/pacoxphVCF \
	-d tempfile.vcf \
	-i tempfile.mlinfo \
	-p $2 \
	-o $3
echo "Run complete."

# rm tempfile.vcf
rm tempfile.mlinfo
rm tempfile.diff

