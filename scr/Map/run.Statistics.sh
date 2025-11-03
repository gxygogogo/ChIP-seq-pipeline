#!/bin/bash

#######################################################################################################################################
######################## Read Pairs ########################
cat list.txt | while read PREFIX
do
echo "${PREFIX}"
# zcat /public2/TangLabData/RawData/ChIP-seq/${PREFIX}/*_R1.fq.gz | grep "@" | wc -l
# zcat /public2/TangLabData/CleanData/ChIP-seq/${PREFIX}/*_R1.fq.gz | grep "@" | wc -l
zcat /public2/TangLabData/CleanData/ChIP-seq/${PREFIX}/*_R1.fq.gz | grep "^@" | wc -l
done


######################## Mapped Read Pairs (*.sorted.bam) ########################
samtools="/data/public/software/samtools-1.3.1/samtools"
cat list.txt | while read PREFIX
do
cd /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/
${samtools} flagstat ${PREFIX}.sorted.bam > ${PREFIX}.sorted.flagstat
#cat flagstat.txt | grep "properly paired" | awk '{print $1/2;}'
done

cat list.txt | while read PREFIX
do
cd /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/
${samtools} flagstat ${PREFIX}.nodup.bam > ${PREFIX}.nodup.flagstat
#cat flagstat.txt | grep "properly paired" | awk '{print $1/2;}'
done

######################## Mapped Read Pairs (*.nodup.bam;Fragments) ########################
cat list.txt | while read PREFIX
do
cd /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/
echo -e "sample\tsorted\tnodup" >align.stat; for f in *.sorted.bam;do file=${f/.sorted.bam/};echo $file;for type in sorted nodup;do  grep "QC-passed reads" $file.$type.flagstat|sed "s/ .*//";done|sed ':a;N;s/\n/\t/;ta;' ;done|sed 'N;s/\n/\t/g'>>align.stat
done

######################## PerC Duplication ########################
cat list.txt | while read PREFIX
do
echo ${PREFIX}
cat /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/${PREFIX}.stat.MarkDuplicates.txt | grep -v "^#" | awk 'NR==3{print $(NF-1);}'
done


#######################################################################################################################################
cat list.txt | while read PREFIX
do
cd /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}
#gzip ${PREFIX}.ext150.cov.bedgraph
#gzip ${PREFIX}.nodup.bedpe
rm ${PREFIX}.sam
#rm ${PREFIX}.bam
#rm ${PREFIX}.sortedByName.bam
rm ${PREFIX}.sorted.bam
mkdir log
mv *.log.txt log/
mkdir statistic
mv *.stat statistic/
mv *.flagstat statistic/
mv *.stat.MarkDuplicates.txt statistic/
done

bedtools bamtobed -bedpe -i GM12878_SA1KOC10_ChIPseq_SA2202_abcamAb_20250414_B5.sorted.bam > GM12878_SA1KOC10_ChIPseq_SA2202_abcamAb_20250414_B5.nodup.bedpe 2>GM12878_SA1KOC10_ChIPseq_SA2202_abcamAb_20250414_B5.bam2bedpe.log.txt

