#!/bin/bash

bwa=/data/public/software/bwa-0.7.15/bwa
samtools=/data/public/software/samtools-1.3.1/samtools
bedtools=/data/public/software/bedtools.2.25.0/bin/bedtools
MarkDuplicates=/data/public/software/picard-tools-1.107/MarkDuplicates.jar
bedGraphToBigWig=/home/yuhan/Software/others/bin/bedGraphToBigWig

## CPU
cpu_num=10

## input
Genome_Index=/data/public/refGenome/bwa_index/hg38/hg38
GENOME_SIZE=/data/public/refGenome/bwa_index/hg38/hg38.chrom.sizes

for PREFIX in GM12878_SA1KOC12_ChIPseq_SA2_202_20231005_B2 GM12878_SA1KOC12_ChIPseq_SA2_202_20231005_B3 GM12878_SA1KOC12_ChIPseq_SA2_202_20231005_B4
do
cd /public2/TangLabData/ProcessedData/ChIP-seq
mkdir ${PREFIX}
cd ${PREFIX}

Fastq_R1=/public2/TangLabData/CleanData/ChIP-seq/${PREFIX}/${PREFIX}_R1.fq.gz
Fastq_R2=/public2/TangLabData/CleanData/ChIP-seq/${PREFIX}/${PREFIX}_R2.fq.gz

echo "
      ${bwa} mem -t ${cpu_num} -k 18 -M ${Genome_Index} ${Fastq_R1} ${Fastq_R2} > ${PREFIX}.sam 2>${PREFIX}.bwa.log.txt
      ${samtools} view -b -o ${PREFIX}.bam -q 20 -@ ${cpu_num} ${PREFIX}.sam 2>${PREFIX}.sam2bam.log.txt

      ${samtools} sort --threads ${cpu_num} -o ${PREFIX}.sorted.bam ${PREFIX}.bam

      java -Xmx40g -jar ${MarkDuplicates} INPUT=${PREFIX}.sorted.bam OUTPUT=${PREFIX}.nodup.bam METRICS_FILE=${PREFIX}.stat.MarkDuplicates.txt REMOVE_DUPLICATES=TRUE ASSUME_SORTED=TRUE VALIDATION_STRINGENCY=SILENT

      ${samtools} sort -n --threads ${cpu_num} -o ${PREFIX}.sortedByName.bam ${PREFIX}.nodup.bam 2>${PREFIX}.sortedByName.log.txt

      ${bedtools} bamtobed -bedpe -i ${PREFIX}.sortedByName.bam > ${PREFIX}.nodup.bedpe 2>${PREFIX}.bam2bedpe.log.txt

      cat ${PREFIX}.nodup.bedpe | \
          perl -lane 'if($F[0] eq $F[3] && $F[0] =~ /^chr[0-9XYAB]+$/){print}' | \
          perl -lane 'unless($F[-2] eq $F[-1]){print join("\t", (@F[0..2], $F[-2])); print join("\t", (@F[3..5], $F[-1]));}' | \
          perl -lane 'if($F[-1] eq "+"){print join("\t", ($F[0], $F[1], $F[1]+250));} else{print join("\t", ($F[0], $F[2]-250, $F[2]));}' | \
          perl -lane 'unless($F[1]<0){print}' | \
          LANG=C sort -k1,1 -k2,2n | \
          ${bedtools} genomecov -bg -i - -g ${GENOME_SIZE} \
          > ${PREFIX}.ext250.cov.bedgraph 2>${PREFIX}.bedpe2cov.log.txt" | qsub -d ./ -N ${PREFIX} -l nodes=node02:ppn=5
done

for PREFIX in GM12878_SA1KOC12_ChIPseq_SA2_202_20231005_B2 GM12878_SA1KOC12_ChIPseq_SA2_202_20231005_B3 GM12878_SA1KOC12_ChIPseq_SA2_202_20231005_B4
do
cd /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}
cat ${PREFIX}.nodup.bedpe | \
          perl -lane 'if($F[0] eq $F[3] && $F[0] =~ /^chr[0-9XYAB]+$/){print}' | \
          perl -lane 'unless($F[-2] eq $F[-1]){print join("\t", (@F[0..2], $F[-2])); print join("\t", (@F[3..5], $F[-1]));}' | \
          perl -lane 'if($F[-1] eq "+"){print join("\t", ($F[0], $F[1], $F[1]+250));} else{print join("\t", ($F[0], $F[2]-250, $F[2]));}' | \
          perl -lane 'unless($F[1]<0){print}' | \
          LANG=C sort -k1,1 -k2,2n | \
          bedtools genomecov -bg -i - -g ${GENOME_SIZE} \
          > ${PREFIX}.ext250.cov.bedgraph 2>${PREFIX}.bedpe2cov.log.txt
done

for PREFIX in GM12878_ChIPseq_SA2_202_20230826_B1 GM12878_ChIPseq_SA2_202_20230826_B2 GM12878_ChIPseq_SA2_202_20230826_B3 GM12878_ChIPseq_SA2_202_20230826_B4 GM12878_ChIPseq_SA2_202_20230826_B5 GM12878_ChIPseq_SA2_202_20230826_B6
do
cd /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}
rm align.stat
done
