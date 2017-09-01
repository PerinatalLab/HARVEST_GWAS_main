#!/bin/bash
set -e

### chr 16 chunks

for c in {1..9}
do
       posstart=$(((c-1) * 11000000+1))
       posstop=$((c * 11000000))
       echo "working on chunk ${c}/9, extracting pos ${posstart}-${posstop}";
       ../erc-genotypes/bin/bcftools view ../merging/moms_16.vcf.gz -r 16:${posstart}-${posstop} -Oz -o maternalfiles7/chunk_${c}.vcf.gz
       ./plink2 --vcf maternalfiles7/chunk_${c}.vcf.gz dosage=DS --export oxford --out maternalfiles7/moms_16_${c}

	# prepare mlinfo file
	echo "rsid REF/EffAl ALT/NonEffAl POS AC INFO Rsq" > maternalfiles7/moms_16_${c}.mlinfo
	awk -v l=${posstart} -v r=${posstop} '$4>=l && $4<=r' BKP_maternal/moms_16.mlinfo >> maternalfiles7/moms_16_${c}.mlinfo
	echo "preparations for chr ${c} complete"

	./gentofvi.R maternalfiles7/moms_16_${c}.gen BKP_maternal/moms_16_catids.sample maternalfiles7/moms_16_${c}
	echo "fvd file created for moms, chunk ${c}"

done

### 2a

posstart=117082876
posstop=260000000
echo "working on chunk 2a, extracting pos 1-${posstart}";
../erc-genotypes/bin/bcftools view ../merging/moms_2.vcf.gz -r 2:1-${posstart} -Oz -o maternalfiles7/chunk_2a.vcf.gz
./plink2 --vcf maternalfiles7/chunk_2a.vcf.gz dosage=DS --export oxford --out maternalfiles7/moms_2_a

# prepare mlinfo file
echo "rsid REF/EffAl ALT/NonEffAl POS AC INFO Rsq" > maternalfiles7/moms_2_a.mlinfo
awk -v l=1 -v r=${posstart} '$4>=l && $4<=r' BKP_maternal/moms_2.mlinfo >> maternalfiles7/moms_2_a.mlinfo
echo "preparations for chr 2a complete"

./gentofvi.R maternalfiles7/moms_2_a.gen BKP_maternal/moms_2_catids.sample maternalfiles7/moms_2_a
echo "fvd file created for moms, chunk 2a"

### 2b

echo "working on chunk 2b, extracting pos ${posstart}-${posstop}";
../erc-genotypes/bin/bcftools view ../merging/moms_2.vcf.gz -r 2:${posstart}-${posstop} -Oz -o maternalfiles7/chunk_2b.vcf.gz
./plink2 --vcf maternalfiles7/chunk_2b.vcf.gz dosage=DS --export oxford --out maternalfiles7/moms_2_b

# prepare mlinfo file
echo "rsid REF/EffAl ALT/NonEffAl POS AC INFO Rsq" > maternalfiles7/moms_2_b.mlinfo
awk -v l=${posstart} -v r=${posstop} '$4>=l && $4<=r' BKP_maternal/moms_2.mlinfo >> maternalfiles7/moms_2_b.mlinfo
echo "preparations for chr 2b complete"

./gentofvi.R maternalfiles7/moms_2_b.gen BKP_maternal/moms_2_catids.sample maternalfiles7/moms_2_b
echo "fvd file created for moms, chunk 2b"

gzip maternalfiles7/moms_*.dose.fvd
rm maternalfiles7/moms_*.gen
