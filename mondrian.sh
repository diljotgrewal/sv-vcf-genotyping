#!/usr/bin/env bash
set -e

ONTBAM=$1
DESTRUCT=$2
LUMPY=$3
SVABA=$4
GRIDSS=$5
OUTDIR=$6

mkdir -p $OUTDIR

python /code/get_mondrian_vcfs.py run --destruct_csv $DESTRUCT --lumpy_vcf $LUMPY --svaba_vcf $SVABA --gridss_vcf $GRIDSS --outdir $OUTDIR


sniffles --input $ONTBAM --genotype-vcf ${OUTDIR}/lumpy.vcf --vcf ${OUTDIR}/lumpy_genotyped_sniffles.vcf
sniffles --input $ONTBAM --genotype-vcf ${OUTDIR}/svaba.vcf --vcf ${OUTDIR}/svaba_genotyped_sniffles.vcf
sniffles --input $ONTBAM --genotype-vcf ${OUTDIR}/gridss.vcf --vcf ${OUTDIR}/gridss_genotyped_sniffles.vcf
