TABLE="/opt/basic/_py/bin/python /opt/basic/console/table_util.py"
TRACK="/opt/basic/_py/bin/python /opt/basic/console/track_util.py"
File=/public1/xinyu/CohesinProject/WAPL/mESC_WT_ChIPseq_0h/mESC_WT_ChIPseq_0h.ext150.bedgraph
${TABLE} create mm10 -l "WAPL ChIP-seq" "mESC_WT_ChIPseq_0h.ext150_cov"
TABLE_ID=6303
${TRACK} gen_cov max ${TABLE_ID} ${File}


File=/public1/xinyu/CohesinProject/WAPL/mESC_WT_ChIPseq_24h/mESC_WT_ChIPseq_24h.ext150.bedgraph
${TABLE} create mm10 -l "WAPL ChIP-seq" "mESC_WT_ChIPseq_24h.ext150_cov"
TABLE_ID=6304
${TRACK} gen_cov max ${TABLE_ID} ${File}
