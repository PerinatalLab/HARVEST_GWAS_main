for c in {1..22} "X"
do
	echo "working on chr ${c}";
	./plink2 --vcf ../merging/moms_${c}.vcf.gz dosage=DS --export oxford --out moms_${c} --threads 1
	gzip moms_${c}.gen
	./plink2 --vcf ../merging/fets_${c}.vcf.gz dosage=DS --export oxford --out fets_${c} --threads 1
	gzip fets_${c}.gen
done
