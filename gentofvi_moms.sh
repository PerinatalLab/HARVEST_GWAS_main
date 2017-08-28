#!/bin/bash

set -e

# argument: chromosome nr
c=$1
INDIR=/media/local-disk2/jjuod

# only the second column in the sample file is used,
# so it should contain FID_IID
# third isn't even necessary
awk 'NR==1{print "id0", "id", "extras"; next} NR==2{next}
	{print 0, $1 "_" $2, 0}' moms_${c}.sample > moms_${c}_catids.sample

# prepare mlinfo file:
# pacoxph reports the effect and freq of REF allele
# (because they flip during conversion to .gen)
# so we should just export it properly:
echo "rsid REF/EffAl ALT/NonEffAl POS AC INFO Rsq" > moms_${c}.mlinfo
${INDIR}/erc-genotypes/bin/bcftools \
	query -f '%ID %REF %ALT %POS %INFO/AC %INFO/INFO 1\n' \
	${INDIR}/merging/moms_${c}.vcf.gz >> moms_${c}.mlinfo
echo "preparations for chr ${c} complete"


## create the actual dose.fvd and dose.fvi files
# gunzip moms_${c}.gen
./gentofvi.R moms_${c}.gen moms_${c}_catids.sample moms_${c}
echo "fvd file created for moms, chr ${c}"

# cleanup:
rm moms_${c}.gen
gzip moms_${c}.dose.fvd
echo "done with moms, chr ${c}"

