#!/bin/bash

set -e

# argument: chromosome nr
c=$1
INDIR=/media/local-disk2/jjuod

# only the second column in the sample file is used,
# so it should contain FID_IID
awk 'NR==1{print "id0", "id"; next} NR==2{next}
	{print 0, $1 "_" $2}' fets_${c}.sample > fets_${c}_catids.sample

# prepare mlinfo file:
# pacoxph reports the effect and freq of REF allele
# (because they flip during conversion to .gen)
# so we should just export it properly:
echo "rsid REF/EffAl ALT/NonEffAl POS AC INFO Rsq" > fets_${c}.mlinfo
${INDIR}/erc-genotypes/bin/bcftools \
	query -f '%ID %REF %ALT %POS %INFO/AC %INFO/INFO 1\n' \
	${INDIR}/merging/fets_${c}.vcf.gz >> fets_${c}.mlinfo
echo "preparations for chr ${c} complete"


## create the actual dose.fvd and dose.fvi files
gunzip fets_${c}.gen
./gentofvi.R fets_${c}.gen fets_${c}_catids.sample fets_${c}
echo "fvd file created for fets, chr ${c}"

# cleanup:
rm fets_${c}.gen
gzip fets_${c}.dose.fvd
echo "done with fets, chr ${c}"
