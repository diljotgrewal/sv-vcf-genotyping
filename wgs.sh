#!/usr/bin/env bash
set -e

ONTBAM=$1
DESTRUCT=$2
LUMPY=$3
OUTDIR=$4

mkdir -p $OUTDIR

python /code/get_wgs_vcf.py run --destruct_csv $DESTRUCT --lumpy_vcf $LUMPY --outdir $OUTDIR

sniffles --input $ONTBAM --genotype-vcf ${OUTDIR}/lumpy.vcf --vcf ${OUTDIR}/lumpy_genotyped_sniffles.vcf
