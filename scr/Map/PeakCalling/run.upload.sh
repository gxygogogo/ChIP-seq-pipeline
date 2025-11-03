#!/bin/bash

alias TABLE="/opt/basic/_py/bin/python /opt/basic/console/table_util.py"
alias TRACK="/opt/basic/_py/bin/python /opt/basic/console/track_util.py"

File=/public2/TangLabData/ProcessedData/ChIP-seq/GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B1S2/Peak/GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B1S2.peaks
TABLE create hg38 -l "GM12878 ChIP-seq" "GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B1S2.peaks"
TABLE_ID=6399
TABLE load ${TABLE_ID} 1:chrom 2:start 3:end 4:name 5:score -i ${File}
TRACK new ${TABLE_ID} scls


File=/public2/TangLabData/ProcessedData/ChIP-seq/GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B2S2/Peak/GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B2S2.peaks
TABLE create hg38 -l "GM12878 ChIP-seq" "GM12878_WT_ChIPseq_WAPL_Proteintech_20230425_B2S2.peaks"
TABLE_ID=6400
TABLE load ${TABLE_ID} 1:chrom 2:start 3:end 4:name 5:score -i ${File}
TRACK new ${TABLE_ID} scls


## Configure
options
{
  "showLabel": false,
  "glyph": {
    "hpad": -1,
    "height": 16.0
},
  "colors": "blue",
  label: "(${chrom}:${start}-${end})",
 tooltip: "(${chrom}:${start}-${end})"
}
series
{}
