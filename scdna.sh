#!/usr/bin/env bash

set -e

ONTBAM=$1
DESTRUCT=$2
LUMPY=$3
OUTDIR=$4

mkdir -p $OUTDIR


cat $LUMPY | awk -F '\t' '{if($1 != "")print $0}' > ${OUTDIR}/lumpy_removed_evidence.bed

SURVIVOR bedpetovcf ${OUTDIR}/lumpy_removed_evidence.bed ${OUTDIR}/lumpy_removed_evidence.vcf


python /code/get_wgs_vcf.py run --destruct_csv $DESTRUCT --lumpy_vcf ${OUTDIR}/lumpy_removed_evidence.vcf --outdir $OUTDIR

sniffles --input $ONTBAM --genotype-vcf ${OUTDIR}/lumpy.vcf --vcf ${OUTDIR}/lumpy_genotyped_sniffles.vcf
