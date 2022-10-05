


```

module load singularity/3.7.1



singularity build image.sif docker://quay.io/diljotgrewal/sv-vcf-genotyping


# mondrian
singularity run --bind /juno image.sif /code/mondrian.sh test_sniffles/sorted.s25.bam mondrian_parsing/results/breakpoint_table.csv mondrian_parsing/results/lumpy.vcf  mondrian_parsing/results/svaba.somatic.sv.vcf.gz mondrian_parsing/results/gridss.vcf.gz mondrian_outdir

#WGS
singularity run --bind /juno image.sif /code/wgs.sh test_sniffles/sorted.s25.bam  destruct_parsing/results/breakpoint_table.csv destruct_parsing/results/lumpy.vcf wgs_outdir


# DLP
singularity run --bind /juno image.sif /code/scdna.sh test_sniffles/sorted.s25.bam dlp_parsing/results/destruct_breakpoints.csv.gz dlp_parsing/results/lumpy_breakpoints.bed scdna_outdir

```
