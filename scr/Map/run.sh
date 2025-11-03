#!/bin/bash
#PBS -N 3
#PBS -l nodes=node04:ppn=10
#PBS -l walltime=10000:00:00
#PBS -j n
#PBS -q batch
#PBS -e ${PBS_JOBNAME}.err
#PBS -o ${PBS_JOBNAME}.out

cd $PBS_O_WORKDIR

##############################################
## Tools
bwa=/data/public/software/bwa-0.7.15/bwa
samtools=/data/public/software/samtools-1.3.1/samtools
bedtools=/data/public/software/bedtools.2.25.0/bin/bedtools
MarkDuplicates=/data/public/software/picard-tools-1.107/MarkDuplicates.jar
bedGraphToBigWig=/home/yuhan/Software/others/bin/bedGraphToBigWig

## CPU
cpu_num=10

## input
Genome_Index=/data/public/refGenome/bwa_index/hg38/hg38
chrSize=/data/public/refGenome/bwa_index/hg38/hg38.chrom.sizes

## Options
EXTENSION=250

PREFIX=21_3_clean
cd /public1/xinyu/tmp/z/analysis
mkdir ${PREFIX}
cd ${PREFIX}
Fastq_R1=/public1/xinyu/tmp/z/CleanFq/${PREFIX}1.fq.gz
Fastq_R2=/public1/xinyu/tmp/z/CleanFq/${PREFIX}2.fq.gz
## mapping
${bwa} mem -t ${cpu_num} -k 18 -M ${Genome_Index} ${Fastq_R1} ${Fastq_R2} > ${PREFIX}.sam 2>${PREFIX}.bwa.log.txt
${samtools} view -b -o ${PREFIX}.bam -q 20 -@ ${cpu_num} ${PREFIX}.sam 2>${PREFIX}.sam2bam.log.txt

## sort bam by coordinate
${samtools} sort --threads ${cpu_num} -o ${PREFIX}.sorted.bam ${PREFIX}.bam

## delete unsorted bam
rm ${PREFIX}.sam

## remove PCR duplications
java -Xmx40g -jar ${MarkDuplicates} INPUT=${PREFIX}.sorted.bam OUTPUT=${PREFIX}.nodup.bam METRICS_FILE=${PREFIX}.stat.MarkDuplicates.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=SILENT

## read bam convered to bed and calculate coverage
${bedtools} bamtobed -i ${PREFIX}.nodup.bam | awk -v EXTENSIONV=$EXTENSION 'BEGIN{OFS="\t"}{if($NF=="+") print $1,$2,$2+EXTENSIONV;else if($NF=="-"){s=$3-EXTENSIONV;if(s>0) print $1,s,$3;else print $1,0,$3;}}' | grep -v [M_] | LANG=C sort -k1,1 -k2,2n > ${PREFIX}.bed
${bedtools} genomecov -bg -i ${PREFIX}.bed -g ${chrSize} > ${PREFIX}.bedgraph
