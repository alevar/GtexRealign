#!/usr/bin/env bash

# ./readlen_stats.sh /ccb/salz4-4/gpertea/GTEX_align/by_tissue/ ./realen.stats

inputDir=${1}
outFP=${2}

touch ${outFP}
echo fp,readlen > ${outFP}

for file in ${inputDir}/*/*cram ; do
    
    readlen="$(samtools view ${file} | cut -f10 | awk '{print length}' | sort | uniq -c | awk '{$1=$1;print}' | tr '\n' ';' | tr ' ' ':')"
    readlen=${readlen%?}
    
    echo ${file},${readlen} >> ${outFP}

done;
