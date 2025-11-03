# ChIP-seq流程
## 染色质免疫沉淀-测序
染色质免疫沉淀-测序（英语：ChIP-sequencing，简称为ChIP-seq）被用于分析蛋白质与DNA的交互作用。该技术将染色质免疫沉淀（ChIP）与大规模并行DNA测序结合起来以鉴定与DNA相关蛋白的结合部位。其可被用于精确绘制任意目的蛋白在全基因组上的结合位点。在此之前，ChIP-on-chip是研究这些蛋白-DNA联系的最常用的技术。
![概况](https://github.com/gxygogogo/ChIP-seq-pipeline/blob/main/img/overview.png)
## ChIP-seq应用
ChIP-seq主要用于确定转录因子和其他染色质相关蛋白质如何影响表型影响机制。确定蛋白质如何与DNA相互作用来调控基因表达对于充分了解许多生物过程和疾病状态至关重要。这种表观遗传信息与基因型和表达分析是互补的。 ChIP-seq技术目前主要被看作是需要杂交阵列的ChIP-on-Chip的替代品。 这必然会引入一些偏差，因为一个阵列仅限于有固定数量的探针。相反，尽管不同测序技术的测序偏差尚未完全了解，测序被认为偏差较小。
与转录因子和其他蛋白质直接物理相互作用的特定DNA位点可通过染色质免疫沉淀被分离出来。 ChIP产生在活体内In vivo兴趣蛋白质结合的目标DNA位点文库。 大规模平行序列分析与全基因组序列数据库结合使用来分析任何蛋白质与DNA的相互作用模式或任何表观遗传染色质修饰的模式。 这可以应用于一系列ChIP蛋白和修饰，如转录因子，聚合酶和转录机制，结构蛋白质，蛋白修饰和DNA修饰。 作为对特异性抗体依赖性的替代方法，已经开发了不同的方法来发现基因组中所有核小体耗尽的(nucleosome-depleted)或核小体扰动(nucleosome-disrupted)的活性调控区的超集，如DNase-Seq和FAIRE-Seq。
## 分析流程
### Step1. Raw data质控
```{shell}
## 安诺过滤(CutTag, ChIP-seq)
# -q 19: 质量值低于此值认为是低质量碱基，默认15，表示质量为Q15
# -u 50: 限制低质量碱基的百分比，默认40，表示40%
# -n 5: 限制N的数量
alias fastp=/home/xinyu/Software/fastp

for PREFIX in GM12878_ndualCUTTag_MCM3RAD21_20251021_B1 GM12878_ndualCUTTag_MCM3RAD21_20251021_B2 GM12878_ndualCUTTag_MCM3RAD21_20251021_B3 GM12878_ndualCUTTag_MCM3RAD21_20251021_B4 GM12878_SA1KOC10_ndualCUTTag_MCM3RAD21_20251021_B1 GM12878_SA1KOC10_ndualCUTTag_MCM3RAD21_20251021_B2 GM12878_SA1KOC10_ndualCUTTag_MCM3RAD21_20251021_B3 GM12878_SA1KOC10_ndualCUTTag_MCM3RAD21_20251021_B4 GM12878_SA2KOA13_ndualCUTTag_MCM3RAD21_20251021_B1 GM12878_SA2KOA13_ndualCUTTag_MCM3RAD21_20251021_B2 GM12878_SA2KOA13_ndualCUTTag_MCM3RAD21_20251021_B3 GM12878_SA2KOA13_ndualCUTTag_MCM3RAD21_20251021_B4 hela_ndualCUTTag_MCM3RAD21_20251021_B1 hela_ndualCUTTag_MCM3RAD21_20251021_B2 hela_ndualCUTTag_MCM3RAD21_20251021_B3 hela_ndualCUTTag_MCM3RAD21_20251021_B4 hela_SA1KOC10_ndualCUTTag_MCM3RAD21_20251021_B1 hela_SA1KOC10_ndualCUTTag_MCM3RAD21_20251021_B2 hela_SA1KOC10_ndualCUTTag_MCM3RAD21_20251021_B3 hela_SA1KOC10_ndualCUTTag_MCM3RAD21_20251021_B4 hela_SA2KOA13_ndualCUTTag_MCM3RAD21_20251021_B1 hela_SA2KOA13_ndualCUTTag_MCM3RAD21_20251021_B2 hela_SA2KOA13_ndualCUTTag_MCM3RAD21_20251021_B3 hela_SA2KOA13_ndualCUTTag_MCM3RAD21_20251021_B4
do
echo ${PREFIX}
mkdir /public2/TangLabData/CleanData/DualCut-Tag/${PREFIX}
cd /public2/TangLabData/CleanData/DualCut-Tag/${PREFIX}
fastp \
   -i /public2/TangLabData/RawData/DualCut-Tag/${PREFIX}/${PREFIX}_R1.fq.gz \
   -o /public2/TangLabData/CleanData/DualCut-Tag/${PREFIX}/${PREFIX}_R1.fq.gz \
   -I /public2/TangLabData/RawData/DualCut-Tag/${PREFIX}/${PREFIX}_R2.fq.gz \
   -O /public2/TangLabData/CleanData/DualCut-Tag/${PREFIX}/${PREFIX}_R2.fq.gz \
   -q 19 \
   -u 50 \
   -n 5
mv *.html ${PREFIX}.html
mv *.json ${PREFIX}.json
done
```
### Step2. 比对分析
```{shell}
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
```
### Step3. 比对统计
```{shell}
#!/bin/bash

# 指定 samtools 路径
samtools="/data/public/software/samtools-1.3.1/samtools"

# 输出文件
output_file="Statistic_output.tsv"
# 添加表头
# echo -e "Library\tCleanReads\tMappedReads\tUniqFragments\tDuplications\n" > ${output_file}
# 读取样本列表
cat list.txt | while read PREFIX
do
    echo "Processing ${PREFIX}"

    # 读取文件路径
    read1="/public2/TangLabData/CleanData/ChIP-seq/${PREFIX}/*_R1.fq.gz"
    # read1="/public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/*_R1.fq.gz"
    sorted_bam="/public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/${PREFIX}.sorted.bam"
    nodup_bedpe="/public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/*.nodup.bedpe"
    markdup_txt="/public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/${PREFIX}.stat.MarkDuplicates.txt"

    # 1. 统计 Read Pairs 数量
    read_pairs_count=$(zcat ${read1} | grep "^@" | wc -l)

    # 2. 统计 Mapped Read Pairs(*.sorted.bam) 数量
    ${samtools} flagstat ${sorted_bam} > /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/flagstat.txt
    mapped_reads_count=$(cat /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/flagstat.txt | grep "properly paired" | awk '{print $1/2;}')

    # 3. 统计可用的 uniq Fragments(*.nodup.bam) 数目
    fragments_count=$(wc -l ${nodup_bedpe} | awk '{print $1;}')

    # 4. 统计 PerC Duplication
    dup_perc=$(grep -v "^#" ${markdup_txt} | awk 'NR==3{print $(NF-1);}')

    # 将所有统计输出到output_file
    echo -e "${PREFIX}\t${read_pairs_count}\t${mapped_reads_count}\t${fragments_count}\t${dup_perc}" >> ${output_file}
done

echo "Script completed. Results saved in ${output_file}"
```
![概况](https://github.com/gxygogogo/ChIP-seq-pipeline/blob/main/img/stat.png)
### Step4. peak calling以及观察
#### 1. Call peak
```{shell}
cat list.txt | while read PREFIX
do
cd /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}
macs14 -t ${PREFIX}.bed -f BED -g hs -s 150 --bw 3500 --slocal 1500 --llocal 50000 -n ${PREFIX}.CallNarrow --keep-dup 1
done
```
#### 2. Peak质量评估
使用run.FRIP.sh进行质量评估
```{shell}
echo -e "sample\treadsInPeak\treads\tpeaks">frip.stat
cat list.txt | while read PREFIX
do
a=`bedtools intersect -a /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/Peak/${PREFIX}.peaks -b /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/${PREFIX}.bed | wc -l`
b=`wc -l /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/${PREFIX}.bed | awk '{print $1}'`
c=`wc -l /public2/TangLabData/ProcessedData/ChIP-seq/${PREFIX}/Peak/${PREFIX}.peaks | awk '{print $1}'`
echo -e ${PREFIX}"\t"$a"\t"$b"\t"$c >> frip.stat
done
cat frip.stat | awk 'BEGIN{OFS="\t"}{if(NR==1) print $0,"FRIP"; else print $0,($2/$3);}' > FRIP_stat.txt
```
#### 3. peak 观察
使用将peak文件转换为bed.gz文件，并且原始的bedgraph文件转换为bigwig文件。
### Step5. peak基因组注释
```{shell}
cat 001_GM12878_SA2KO_CutTag_NIPBL.bed | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > SA2KO_only_sorted.bed
cat 010_GM12878_SA1KO_CutTag_NIPBL.bed | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > SA1KO_only_sorted.bed
cat 011_GM12878_SA1KO_CutTag_NIPBL_GM12878_SA2KO_CutTag_NIPBL.bed | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > SA1KO_SA2KO_sorted.bed
cat 100_GM12878_CutTag_NIPBL.bed | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > WT_only_sorted.bed
cat 101_GM12878_CutTag_NIPBL_GM12878_SA2KO_CutTag_NIPBL.bed | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > WT_SA2KO_sorted.bed
cat 110_GM12878_CutTag_NIPBL_GM12878_SA1KO_CutTag_NIPBL.bed | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > WT_SA1KO_sorted.bed
cat 111_GM12878_CutTag_NIPBL_GM12878_SA1KO_CutTag_NIPBL_GM12878_SA2KO_CutTag_NIPBL.bed | sort -k1,1 -k2,2n | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1"_"$2"_"$3}' > common_sorted.bed

bedtools intersect -a WT_only_sorted.bed -b /data1/yuhan/BasicBrowser/BASIC_ChromhmmGm12878/ChromHMM.bed -sorted -wa -wb > tmp1
cut -f4,9 tmp1 | sort | uniq > tmp2
awk 'BEGIN{ORS="\t";TMP="TMP"}{if($1==TMP) print $2;else print "\n"$0;TMP=$1}' tmp2 | awk 'BEGIN{OFS="\t";}NR!=1{if($0~"DistalNIPBL/CandidateInsulator") print $1,"DistalNIPBL/CandidateInsulator";else if($0~"ActivePromoter") print $1,"ActivePromoter";else if($0~"CandidateStrongEnhancer") print $1,"CandidateStrongEnhancer";else if($0~"CandidateWeakEnhancer/DNase") print $1,"CandidateWeakEnhancer/DNase";else if($0~"TranscriptionAssociated") print $1,"TranscriptionAssociated";else if($0~"Heterochromatin/Repetitive/CopyNumberVariation") print $1,"Heterochromatin/Repetitive/CopyNumberVariation";else if($0~"InactivePromoter") print $1,"InactivePromoter";else if($0~"LowActivityProximalToActive") print $1,"LowActivityProximalToActive";else if($0~"PolycombRepressed") print $1,"PolycombRepressed";}' > WT_only_chromHMM.txt
```
### Step6. peak下游分析
依据peak文件和bam文件，计算每个peak区域的覆盖强度，作为count值。
```{shell}
WT_H3K27me3=/public1/xinyu/CohesinProject/FinalizedCutTag/HDAC1/peaks/GM12878_WT_CutTag_HDAC1.peaks
SA1KO_HDAC1=/public1/xinyu/CohesinProject/FinalizedCutTag/HDAC1/peaks/GM12878_SA1KO_CutTag_HDAC1.peaks
SA2KO_HDAC1=/public1/xinyu/CohesinProject/FinalizedCutTag/HDAC1/peaks/GM12878_SA2KO_CutTag_HDAC1.peaks

cat ${WT_BRD4} ${SA1KO_BRD4} ${SA2KO_BRD4} \
    ${WT_H3K27ac} ${SA1KO_H3K27ac} ${SA2KO_H3K27ac} \
    ${WT_NIPBL} ${SA1KO_NIPBL} ${SA2KO_NIPBL} \
    ${WT_P300} ${SA1KO_P300} ${SA2KO_P300} \
    ${WT_HDAC1} ${SA1KO_HDAC1} ${SA2KO_HDAC1} \
    ${WT_CTCF} ${SA1KO_CTCF} ${SA2KO_CTCF} | \
   sort -k1,1 -k2,2n | uniq | bedtools merge -i - | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$1":"$2"-"$3}' > GM12878_CutTag_peakRegion.bed



factor=HDAC1
for PREFIX in GM12878_CutTag_${factor} GM12878_SA1KO_CutTag_${factor} GM12878_SA2KO_CutTag_${factor}
do
bedtools coverage -a GM12878_CutTag_peakRegion.bed -b /public1/xinyu/CohesinProject/FinalizedCutTag/${factor}/${PREFIX}.bam | cut -f4-5 > ${PREFIX}.counts
done
```
使用DESeq2计算差异peak，然后进行下游分析。
