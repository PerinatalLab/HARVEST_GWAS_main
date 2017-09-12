#!/bin/bash

## This script is meant to provide an easy interface for calling
## ProbABEL's Cox PH module on VCF files.
## pacoxphVCF requires a raw .vcf file and an .mlinfo file,
## which are created by this script.
## Requires bcftools.
## For conditioned analysis, just prepare .pheno file with additional SNP column.
## TODO: interaction analysis.

## USAGE: ./run_pacoxphVCF.sh infile.vcf.gz phenofile outfile tmpdir chr:regionstart-regionstop

set -e

export PATH=/media/local-disk2/jjuod/erc-genotypes/bin:$PATH

echo "Checking input..."
if [ "$#" -ne 5 ]
then
	echo "Provided $# parameters instead of 5."
	echo "USAGE: ./run_pacoxphVCF.sh infile.vcf.gz phenofile outfile tmpdir chr:regionstart-regionstop"
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

# temporary hack to split into sub-10k snp chunks:
posmid=$(( posstart + 1000000 ))
bothposstart=($posstart $posmid)
bothposstop=($posmid $posstop)

for i in {0..1}
do
	posstart=${bothposstart[$i]}
	posstop=${bothposstop[$i]}

	 echo "Checking region size..."
	 if [ $((posstop-posstart)) -gt 10000000 ]
	 then
	 	echo "Error: max allowed region size is 10 Mbp. You requested $((posstop-posstart))."
	 	exit
	 fi
	 echo "Working on chr ${chr}, region ${posstart}-${posstop}"
	 
	 echo "Creating .mlinfo file..."
	 echo "rsid REF/NonEffAl ALT/EffAl POS AC INFO Rsq" > ${4}/tempfile_${i}.mlinfo
	 bcftools query $1 \
	 	-f '%ID %REF %ALT %POS %INFO/AC %INFO/INFO 1\n' \
	 	-r ${chr}:${posstart}-${posstop} >> ${4}/tempfile_${i}.mlinfo
	 echo "MLINFO file created."
	 
	 echo "Extracting region ${chr}:${posstart}-${posstop} from the VCF file..."
	 bcftools view $1 -r ${chr}:${posstart}-${posstop} -Ov -o ${4}/tempfile_${i}.vcf
	 echo "Region extracted."
	 
	 echo "Calling ProbABEL..."
	 ../bin/pacoxphVCF \
	 	-d ${4}/tempfile_${i}.vcf \
	 	-i ${4}/tempfile_${i}.mlinfo \
	 	-p $2 \
	 	-o ${3}_${i}
	 echo "Run complete."
done

# rm tempfile.vcf
# rm tempfile_0.mlinfo
# rm tempfile_1.mlinfo
# rm tempfile.diff

