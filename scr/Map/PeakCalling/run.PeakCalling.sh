#!/bin/bash
# ssh node02
cat list.txt | while read PREFIX
do
cd /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}
macs14 -t ${PREFIX}.bed -f BED -g hs -s 150 --bw 3500 --slocal 1500 --llocal 50000 -n ${PREFIX}.CallNarrow --keep-dup 1
done

cat list.txt | while read PREFIX
do
cd /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}
echo ${PREFIX}
awk 'BEGIN{OFS="\t";}$1~"chr"&&$2!="start"{print $1,$2,$3;}' ${PREFIX}.CallNarrow_peaks.xls | wc -l
awk 'BEGIN{OFS="\t";}$1~"chr"&&$2!="start"{print $1,$2,$3;}' ${PREFIX}.CallNarrow_peaks.xls > ${PREFIX}.peaks
done







for PREFIX in 21_2_clean 21_3_clean
do
cd /public1/xinyu/tmp/z/analysis/${PREFIX}
macs14 -t ${PREFIX}.bed -f BED -g hs -s 150 --bw 3500 --slocal 1500 --llocal 50000 -n ${PREFIX}.CallNarrow --keep-dup 1
done

for PREFIX in 21_2_clean 21_3_clean
do
cd /public1/xinyu/tmp/z/analysis/${PREFIX}
echo ${PREFIX}
awk 'BEGIN{OFS="\t";}$1~"chr"&&$2!="start"{print $1,$2,$3;}' ${PREFIX}.CallNarrow_peaks.xls | wc -l
awk 'BEGIN{OFS="\t";}$1~"chr"&&$2!="start"{print $1,$2,$3;}' ${PREFIX}.CallNarrow_peaks.xls > ${PREFIX}.peaks
done