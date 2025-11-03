TABLE="/opt/basic/_py/bin/python /opt/basic/console/table_util.py"
TRACK="/opt/basic/_py/bin/python /opt/basic/console/track_util.py"
File=/public2/TangLabData/ProcessedData/ChIP-seq/GM12878_SA2KOA13_SA2201_ChIPseq_SA2Clone61Ab_20230322_B1S2/GM12878_SA2KOA13_SA2201_ChIPseq_SA2Clone61Ab_20230322_B1S2.bedgraph
${TABLE} create hg38 -l "GM12878 ChIP-seq" "GM12878_SA2KOA13_SA2201_ChIPseq_SA2Clone61Ab_20230322_B1S2.ext150_cov"
TABLE_ID=6292
${TRACK} gen_cov max ${TABLE_ID} ${File}


## WAPL

TABLE="/opt/basic/_py/bin/python /opt/basic/console/table_util.py"
TRACK="/opt/basic/_py/bin/python /opt/basic/console/track_util.py"
File=/public2/TangLabData/ProcessedData/ChIP-seq/GM12878_WT_ChIPseq_WAPL_Bethyl_20230425_B1/GM12878_WT_ChIPseq_WAPL_Bethyl_20230425_B1.bedgraph
${TABLE} create hg38 -l "GM12878 ChIP-seq" "GM12878_WT_ChIPseq_WAPL_Bethyl_20230425_B1.cov"
TABLE_ID=6361
${TRACK} gen_cov max ${TABLE_ID} ${File}

File=/public2/TangLabData/ProcessedData/ChIP-seq/GM12878_WT_ChIPseq_WAPL_Bethyl_20230425_B2/GM12878_WT_ChIPseq_WAPL_Bethyl_20230425_B2.bedgraph
${TABLE} create hg38 -l "GM12878 ChIP-seq" "GM12878_WT_ChIPseq_WAPL_Bethyl_20230425_B2.cov"
TABLE_ID=6362
${TRACK} gen_cov max ${TABLE_ID} ${File}

File=/public2/TangLabData/ProcessedData/ChIP-seq/GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B1/GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B1.bedgraph
${TABLE} create hg38 -l "GM12878 ChIP-seq" "GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B1.cov"
TABLE_ID=6363
${TRACK} gen_cov max ${TABLE_ID} ${File}

File=/public2/TangLabData/ProcessedData/ChIP-seq/GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B2/GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B2.bedgraph
${TABLE} create hg38 -l "GM12878 ChIP-seq" "GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B2.cov"
TABLE_ID=6364
${TRACK} gen_cov max ${TABLE_ID} ${File}


File=/public2/TangLabData/ProcessedData/ChIP-seq/GM12878_ChIPseq_SA2202_abcamAb_20250414_B6/GM12878_ChIPseq_SA2202_abcamAb_20250414_B6.bedgraph
${TABLE} create hg38 -l "GM12878 ChIP-seq" "GM12878_ChIPseq_SA2202_abcamAb_20250414_B6.ext250_cov"
TABLE_ID=8126
${TRACK} gen_cov max ${TABLE_ID} ${File}


