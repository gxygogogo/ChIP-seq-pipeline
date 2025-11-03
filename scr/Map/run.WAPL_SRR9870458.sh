#!/bin/bash
#PBS -N WAPL_SRR9870459
#PBS -l nodes=node02:ppn=15
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
cpu_num=15

## input
GENOME=/data/public/refGenome/bwa_index/mm10/mm10
GENOME_SIZE=/data/public/refGenome/bwa_index/mm10/mm10.chrom.sizes

#####################################################################
PREFIX=mESC_WT_ChIPseq_24h

Fastq=/public1/xinyu/CohesinProject/WAPL/mESC_WT_ChIPseq_24h/SRR9870459.fastq.gz
${bwa} mem -t ${cpu_num} -k 18 -M ${GENOME} ${Fastq} > ${PREFIX}.sam 2>${PREFIX}.bwa.log.txt
${samtools} view -b -o ${PREFIX}.bam -q 20 -@ ${cpu_num} ${PREFIX}.sam 2>${PREFIX}.sam2bam.log.txt

## sort bam by coordinate
${samtools} sort --threads ${cpu_num} -o ${PREFIX}.sorted.bam ${PREFIX}.bam

## delete unsorted bam
rm ${PREFIX}.bam
## remove PCR duplications
java -Xmx40g -jar ${MarkDuplicates} INPUT=${PREFIX}.sorted.bam OUTPUT=${PREFIX}.nodup.bam METRICS_FILE=${PREFIX}.stat.MarkDuplicates.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=SILENT
${bedtools} bamtobed -i ${PREFIX}.nodup.bam |\
     awk -v EXTENSIONV=150 'BEGIN{OFS="\t"}{if($NF=="+") print $1,$2,$2+EXTENSIONV;else if($NF=="-"){s=$3-EXTENSIONV;if(s>0) print $1,s,$3;else print $1,0,$3;}}' |\
     grep -v [M_] |\
     LANG=C sort -k1,1 -k2,2n > ${PREFIX}.ext150.bed

awk '($1!~/_/)&&($1!~/M/)&&($1!~/EBV/){print }' ${PREFIX}.ext150.bed | ${bedtools} genomecov -bg -i - -g ${GENOME_SIZE} > ${PREFIX}.ext150.bedgraph
