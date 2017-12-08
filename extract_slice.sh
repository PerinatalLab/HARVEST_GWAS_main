#!/bin/bash

set -e

# USAGE: ./extract_slice.sh CHROM START END

PATH=$PATH:/media/local-disk2/jjuod/erc-genotypes/bin

chr=$1
start=$2
end=$3
outfile=/media/local-disk2/jjuod/other/slices/moms_${chr}_${start}-${end}.txt.gz

echo "extracting slice $chr $start - $end"

bcftools query -l moms_${chr}.vcf.gz > /media/local-disk2/jjuod/other/slices/samplelist_moms.txt
bcftools view moms_${chr}.vcf.gz \
	-r ${chr}:${start}-${end} \
	-q 0.01:minor | \
	bcftools query \
	-f '%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n' | \
	gzip -c > ${outfile}

