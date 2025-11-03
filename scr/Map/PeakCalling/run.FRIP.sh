#!/bin/bash

echo -e "sample\treadsInPeak\treads\tpeaks">frip.stat
cat list.txt | while read PREFIX
do
a=`bedtools intersect -a /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/Peak/${PREFIX}.peaks -b /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/${PREFIX}.bed | wc -l`
b=`wc -l /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/${PREFIX}.bed | awk '{print $1}'`
c=`wc -l /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/Peak/${PREFIX}.peaks | awk '{print $1}'`
echo -e ${PREFIX}"\t"$a"\t"$b"\t"$c >> frip.stat
done
cat frip.stat | awk 'BEGIN{OFS="\t"}{if(NR==1) print $0,"FRIP"; else print $0,($2/$3);}' > FRIP_stat.txt
